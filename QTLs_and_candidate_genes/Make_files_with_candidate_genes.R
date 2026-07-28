### This script processes the output of Collect_AT_IDs_from_NCBI.R to make the file with candidate
### genes found in the supplemental data and also the supplemental file with Arabidopsis thaliana
### homologs for all sequences in the gff.

library(dplyr)

### Here we make the file with the candidate genes for the 16S traits.

### Read in A. thaliana IDs
load("ara.hom.16s.v2.out")

nrs <- sapply(gregexpr("\\d", bac.hom$Gene_name), function(m) if (m[1] == -1) 0 else length(m))
bac.hom$Gene_name[nrs>3] <- ""

### Read in blast output
bac <-  read.table("16S_blast_results_final.fmt6")
colnames(bac) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

### Add A. thaliana homolog information to the BLAST output
bac$AT_ID <- bac.hom$AT_ID[match(bac$NCBI_ID,bac.hom$Accession_ID)]
bac$AT_synonym <- bac.hom$Gene_name[match(bac$NCBI_ID,bac.hom$Accession_ID)]

### Make a dataframe with one row per sequence
bac.frame <- data.frame(bac %>%
  mutate(bitscore = as.numeric(as.character(bitscore))) %>%
  group_by(Gene) %>%
  arrange(desc(bitscore)) %>%
  slice_head(n = 10) %>%
  summarise(
    NCBI_ID    = paste(NCBI_ID, collapse = "|"),
    E_value    = paste(E_value, collapse = "|"),
    bitscore   = paste(bitscore, collapse = "|"),
    AT_ID      = paste(AT_ID, collapse = "|"),
    AT_synonym = paste(AT_synonym, collapse = "|"),
    .groups = "drop"
  ))

meta.bac <- strsplit(bac.frame$Gene,";")
meta.bac  <-  data.frame(do.call(rbind, meta.bac))
colnames(meta.bac) <- c("Chromosome","Gene_start","Gene_end","Phenotype","Compartment","Peak_SNP_location","P-value")

final.bac <- data.frame(meta.bac,bac.frame[-1])

### Save data
write.csv(final.bac,"Arabidopsis_homologs_near_QTLs_16S.csv")

### Subset the dataframe to only those rows where the peak SNP is inside the gene sequence
bac.in <- final.bac[final.bac$Peak_SNP_location > final.bac$Gene_start & final.bac$Peak_SNP_location < final.bac$Gene_end,]

write.csv(bac.in,file = "Arabidopsis_homologs_near_QTLs_16S_in_genes.csv")

### Now do the same for the ITS data

load("ara.hom.its.v2.out")

nrs.its <- sapply(gregexpr("\\d", its.hom$Gene_name), function(m) if (m[1] == -1) 0 else length(m))
its.hom$Gene_name[nrs.its>3] <- ""

its <-  read.table("ITS_blast_results_final.fmt6")
colnames(its) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

its$AT_ID <- its.hom$AT_ID[match(its$NCBI_ID,its.hom$Accession_ID)]
its$AT_synonym <- its.hom$Gene_name[match(its$NCBI_ID,its.hom$Accession_ID)]

head(its)

its.frame <- data.frame(its %>%
                          mutate(bitscore = as.numeric(as.character(bitscore))) %>%
                          group_by(Gene) %>%
                          arrange(desc(bitscore)) %>%
                          slice_head(n = 10) %>%
                          summarise(
                            NCBI_ID    = paste(NCBI_ID, collapse = "|"),
                            E_value    = paste(E_value, collapse = "|"),
                            bitscore   = paste(bitscore, collapse = "|"),
                            AT_ID      = paste(AT_ID, collapse = "|"),
                            AT_synonym = paste(AT_synonym, collapse = "|"),
                            .groups = "drop"
                          ))

meta.its <- strsplit(its.frame$Gene,";")
meta.its  <-  data.frame(do.call(rbind, meta.its))
colnames(meta.its) <- c("Chromosome","Gene_start","Gene_end","Phenotype","Compartment","Peak_SNP_location","P-value")

final.its <- data.frame(meta.its,its.frame[-1])

write.csv(final.its,"Arabidopsis_homologs_near_QTLs_ITS.csv")

its.in <- final.its[final.its$Peak_SNP_location > final.its$Gene_start & final.its$Peak_SNP_location < final.its$Gene_end,]

write.csv(its.in,file = "Arabidopsis_homologs_near_QTLs_ITS_in_genes.csv")

### Finally, make the file containing the homologs of all sequences in the genome

load("ara.hom.all.out")

seq <-  read.table("All_sequences_blast_results.fmt6")
seq <- seq[!duplicated(seq),]

colnames(seq) <- c("Gene","NCBI_ID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","E_value","bitscore")

seq_accessions <- as.vector(sapply(seq$NCBI_ID,function(x) strsplit(x,".",fixed = T)[[1]][1]))

seq.full <- unique(seq$NCBI_ID)

nrs <- sapply(gregexpr("\\d", seq.hom$Gene_name), function(m) if (m[1] == -1) 0 else length(m))
seq.hom$Gene_name[nrs>3] <- ""

seq$AT_ID <- seq.hom$AT_ID[match(seq$NCBI_ID,seq.hom$Accession_ID)]
seq$AT_synonym <- seq.hom$Gene_name[match(seq$NCBI_ID,seq.hom$Accession_ID)]

frame <- data.frame(seq %>%
                      mutate(bitscore = as.numeric(as.character(bitscore))) %>%
                      group_by(Gene) %>%
                      arrange(desc(bitscore)) %>%
                      slice_head(n = 10) %>%
                      summarise(
                        NCBI_ID    = paste(NCBI_ID, collapse = "|"),
                        E_value    = paste(E_value, collapse = "|"),
                        bitscore   = paste(bitscore, collapse = "|"),
                        AT_ID      = paste(AT_ID, collapse = "|"),
                        AT_synonym = paste(AT_synonym, collapse = "|"),
                        .groups = "drop"
                      ))

meta <- strsplit(frame$Gene,"|",fixed = T)
meta  <-  data.frame(do.call(rbind, meta))
colnames(meta) <- c("Chromosome","Gene_start","Gene_end","Type","ID")

final <- data.frame(meta,frame[-1])

write.xlsx(final,"homologs_all_seqs.xlsx")










