### This script uses the rentrez package to scrape NCBI metadata from the NCBI IDs output by blastx.
### Specifically, it collects the TAIR IDs and gene names of the homologous sequences detected by
### blastx.

library(qdapRegex)
library(rentrez)

### Read in the BLAST output of 16S sequences
bac <-  read.table("16S_blast_results_final.fmt6")
bac <- bac[!duplicated(bac),]
colnames(bac) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

### Collect all unique NCBI IDs in the blast output
bac_accessions <- as.vector(sapply(bac$NCBI_ID,function(x) strsplit(x,".",fixed = T)[[1]][1]))
bac.full <- unique(bac$NCBI_ID)


### The below two functions together collect the necessary information. This code requires a lot of
### failsafes and delays because the NCBI webserver gives an error when receiving too many requests.
fetch_with_retry <- function(id_string, max_retries = 6, delay_base = 2) {
  attempt <- 1
  repeat {
    res <- tryCatch(
      rentrez::entrez_fetch(db = "protein", rettype = "txt", id = id_string),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
    message("Fetch error: ", conditionMessage(res), " (attempt ", attempt, " of ", max_retries, ")")
    if (attempt >= max_retries) stop(res)
    Sys.sleep(delay_base * 2^(attempt - 1))   # exponential backoff
    attempt <- attempt + 1
  }
}

parse_record <- function(fetch_txt, accession) {
  recs <- strsplit(fetch_txt, "\n//\n", fixed = TRUE)[[1]]
  out <- vector("list", length(recs))
  for (i in seq_along(recs)) {
    record <- strsplit(recs[[i]], "\n")[[1]]
    lt_line <- record[grep('/locus_tag=', record, fixed = TRUE)]
    if (length(lt_line) == 0) {
      ara_id <- ""
    } else {
      ara_id <- gsub("[^A-Za-z0-9.]", "", lt_line)
      ara_id <- strsplit(ara_id,"locustag",fixed = T)[[1]][2]
    }
    
    if(any(grepl("gene=",record))){
      gene.names <- record[grep("gene=",record)]
      gene <- strsplit(gene.names,"\"")[[1]][2]
    } else { 
      gs_lines <- record[grep('/gene_synonym=', record, fixed = TRUE)]
      if (length(gs_lines) == 0) {
        gene <- ""
      } else {
        gene.name <- gs_lines[!grepl(".",gs_lines,fixed = T)]
        gen.length <- nchar(gene.name)
        gene.name <- gene.name[which.min(gen.length)]
        gene <- if (length(gene.name)) gene.name[which.min(nchar(gene.name))] else ""
        gene <- gsub("[^A-Za-z0-9.]", "", gene.name)
        if(length(gene)==0){
          gene <- ""
        }  else{
          gene <- strsplit(gene,"genesynonym",fixed = T)[[1]][2]
        }
        
      }
    }
    
    out[[i]] <- data.frame(Accession_ID = if (length(accession) > 1) accession[i] else accession,
                           AT_ID = ara_id, Gene_name = gene, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

ids <- bac.full
batch_size <- 200
delay <- 2   # seconds between batches; increase if you still get rate limits
results <- vector("list", ceiling(length(ids)/batch_size))
k <- 1

### Call these functions in a for loop
for (start in seq(1, length(ids), by = batch_size)) {
  b <- ids[start:min(start + batch_size - 1, length(ids))]
  id_string <- paste(b, collapse = ",")
  res_txt <- tryCatch(fetch_with_retry(id_string), error = function(e) e)
  if (inherits(res_txt, "error")) {
    message("Batch failed for IDs: ", paste(b, collapse = ","), " — retrying individually")
    # try individuals
    for (single in b) {
      res1 <- tryCatch(fetch_with_retry(single), error = function(e) e)
      if (inherits(res1, "error")) {
        warning("Failed to fetch ", single, ": ", conditionMessage(res1))
        next
      }
      results[[k]] <- parse_record(res1, single)
      k <- k + 1
      Sys.sleep(delay)
    }
  } else {
    results[[k]] <- parse_record(res_txt, b)
    k <- k + 1
    Sys.sleep(delay)
  }
}

### Put output in a dataframe
bac.hom <- do.call(rbind, results[seq_len(k-1)])

save(bac.hom,file = "ara.hom.16s.v2.out")

### Read in the BLAST output of ITS sequences and repeat
its <- read.table("ITS_blast_results_final.fmt6")
its <- its[!duplicated(its),]

colnames(its) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

its_accessions <- as.vector(sapply(its$NCBI_ID,function(x) strsplit(x,".",fixed = T)[[1]][1]))

its.full <- unique(its$NCBI_ID)

ids <- its.full
batch_size <- 200
delay <- 2  
results <- vector("list", ceiling(length(ids)/batch_size))
k <- 1

for (start in seq(1, length(ids), by = batch_size)) {
  b <- ids[start:min(start + batch_size - 1, length(ids))]
  id_string <- paste(b, collapse = ",")
  res_txt <- tryCatch(fetch_with_retry(id_string), error = function(e) e)
  if (inherits(res_txt, "error")) {
    message("Batch failed for IDs: ", paste(b, collapse = ","), " — retrying individually")
    # try individuals
    for (single in b) {
      res1 <- tryCatch(fetch_with_retry(single), error = function(e) e)
      if (inherits(res1, "error")) {
        warning("Failed to fetch ", single, ": ", conditionMessage(res1))
        next
      }
      results[[k]] <- parse_record(res1, single)
      k <- k + 1
      Sys.sleep(delay)
    }
  } else {
    results[[k]] <- parse_record(res_txt, b)
    k <- k + 1
    Sys.sleep(delay)
  }
}

its.hom <- do.call(rbind, results[seq_len(k-1)])

save(its.hom,file = "ara.hom.its.v2.out")

### And finally, we did the same for all sequences in the GFF

seq <-  read.table("All_sequences_blast_results.fmt6")
seq <- seq[!duplicated(seq),]

colnames(seq) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

seq_accessions <- as.vector(sapply(seq$NCBI_ID,function(x) strsplit(x,".",fixed = T)[[1]][1]))

seq.full <- unique(seq$NCBI_ID)

ids <- seq.full
batch_size <- 200
delay <- 2   # seconds between batches; increase if you still get rate limits
results <- vector("list", ceiling(length(ids)/batch_size))
k <- 1

for (start in seq(1, length(ids), by = batch_size)) {
  b <- ids[start:min(start + batch_size - 1, length(ids))]
  id_string <- paste(b, collapse = ",")
  res_txt <- tryCatch(fetch_with_retry(id_string), error = function(e) e)
  if (inherits(res_txt, "error")) {
    message("Batch failed for IDs: ", paste(b, collapse = ","), " — retrying individually")
    # try individuals
    for (single in b) {
      res1 <- tryCatch(fetch_with_retry(single), error = function(e) e)
      if (inherits(res1, "error")) {
        warning("Failed to fetch ", single, ": ", conditionMessage(res1))
        next
      }
      results[[k]] <- parse_record(res1, single)
      k <- k + 1
      Sys.sleep(delay)
    }
  } else {
    results[[k]] <- parse_record(res_txt, b)
    k <- k + 1
    Sys.sleep(delay)
  }
}

seq.hom <- do.call(rbind, results[seq_len(k-1)])

save(seq.hom,file = "ara.hom.all.out")



