### This script makes fasta's to be used by blastx with the aim of identifying Arabidopsis homologs.
### First, we make fastas for the genes in a 0.2Mbp interval around each peak SNP. At the end of the
### script, we also made a fasta including all sequences in the gff labelled as "gene" or "non-coding RNA".
### In the annotation of this script we also show the command line code we used to blast the sequences

library(Biostrings)
library(stringr)
library(seqinr)

### Load in the genome assembly (from Byrne et al., 2015 https://doi.org/10.1111/tpj.13037)
assembly <- readDNAStringSet("Lolium_2.6.1_Chromosomes_Scaffolds.fasta")
gff <- read.table("Lolium_2.6.1_V3_gene_annotation.gff3")
colnames(gff)[c(1,3,4,5)] <- c("Chromosome","Type","Start","Stop")
gff$Chromosome[gff$Chromosome=="Lp_chr7_0.1"] <- "Lp_chr7_0_1"

#### Read in the genus-rank relative abundances (to select the samples with available data)

load("agglomerated_counts/genus_rel_filt_mitchlo.out")
genus.rel<-genus.rel[,-which(colnames(genus.rel)%in%c("GWr_170","GWs_179","GWr_065","GWs_018","GWs_074","GWs_001","GWs_002","GWs_120","GWr_014","GWr_189","GWs_000"))]
genus.rel.root <- genus.rel[grep("GWr",colnames(genus.rel))]
genus.rel.soil <- genus.rel[,grep("GWs",colnames(genus.rel))]
colnames(genus.rel.soil)[which(colnames(genus.rel.soil)=="GWs_215a")] <- "GWs_215"
in.root <- as.numeric(gsub("GWr_","",colnames(genus.rel.root)))
in.soil <- as.numeric(gsub("GWs_","",colnames(genus.rel.soil)))
genus.rel.soil <- genus.rel.soil[,!is.na(in.soil)]
in.soil <- in.soil[!is.na(in.soil)]

### Read in all the peak SNPs from all 16S endosphere traits 
load("Submission_data/16S_endo_submission.out")
pheno <- plot.frame$Phenotype
max.snps <- plot.frame$max_snps
chr <- plot.frame$Chromosome
pos <- plot.frame$Position
pvalues <- plot.frame$Pval

### Read in the SNP metadata
load("~/Documents/Lolium/LoliumSNPs/Snp_info_submission/snp.info.full.out")

pheno <- str_replace(pheno," ","_")

length_chr <- c(length(assembly[[1]]),length(assembly[[2]]),length(assembly[[3]]),length(assembly[[4]]),length(assembly[[5]]),length(assembly[[6]]),
                length(assembly[[7]]))

chr_names <- c("Lp_chr1_0","Lp_chr2_0","Lp_chr3_0","Lp_chr4_0","Lp_chr5_0","Lp_chr6_0","Lp_chr7_0_1")

### This function takes all sequences in a 0.2 interval around the peak SNPs out of the genome 
### assembly (using the gff coordinates)
get_interval_sequence <- function(pos,chro,span,phenotype,habitat,pvalue){
  sequences <- c()
  interval <- span*10^6
  loc <- pos*10^6
  min.int <- max(0,loc-interval)
  max.int <- min(loc+interval,length_chr[as.numeric(chro)])
  seq <- gff[((gff$Start<min.int&gff$Stop>min.int)|(gff$Start>min.int & gff$Stop<max.int )|(gff$Start<max.int&gff$Stop>max.int))&gff$Chromosome==chr_names[as.numeric(chro)],c(1,3,4,5)]
  seq <- seq[seq$Type %in% c("gene","ncRNA"),]
  if(nrow(seq>0)){
    for(i in 1:nrow(seq)){
      sequences <- c(sequences, as.character(assembly[[as.numeric(chro)]][seq$Start[i]:seq$Stop[i]]))
    }
    names(sequences) <- paste(chro,seq$Start,seq$Stop,phenotype,habitat,loc,pvalue,sep=";")
    return(sequences)
  }
}

### Call function
seq_16S_root <- sapply(1:length(chr),function(x) get_interval_sequence(pos[x],chr[x],0.1,pheno[x],"endosphere",pvalues[x]))
seq_16S_root <- unlist(seq_16S_root)
seq.names.16S.root <- names(seq_16S_root)
seq_16S_root <- as.list(seq_16S_root)
seq_16S_root <- unname(seq_16S_root)

### Now the same thing for the 16S rhizosphere traits

load("Submission_data/16S_rhizo_submission.out")
pheno <- plot.frame$Phenotype
max.snps <- plot.frame$max_snps
chr <- plot.frame$Chromosome
pos <- plot.frame$Position
pvalues <- plot.frame$Pval

seq_16S_soil <- sapply(1:length(chr),function(x) get_interval_sequence(pos[x],chr[x],0.1,pheno[x],"rhizosphere",pvalues[x]))
seq_16S_soil <- unlist(seq_16S_soil)
seq.names.16S.soil <- names(seq_16S_soil)
seq_16S_soil <- as.list(seq_16S_soil)
seq_16S_soil <- unname(seq_16S_soil)

### Now put the endosphere and rhizosphere sequences in the same vector

seq_16S <- c(seq_16S_root,seq_16S_soil)
seq.names.16S <- c(seq.names.16S.root,seq.names.16S.soil)

bp.tot <- lapply(seq_16S,function(x) nchar(x))

### Write to a fasta!
write.fasta(sequences = seq_16S,names = seq.names.16S,file.out = "~/Documents/Lolium/LoliumSNPs/Submission_data/16S_seqs_final.fasta",open = "w")

### Repeat the above for the ITS data

OTU.gwas<-read.table("otutable_fungigwasfinal.txt")
read.counts <- apply(OTU.gwas,2,sum)

load("agglomerate_counts_fungi/genus_rel_filt_mitchlo.out")
genus.rel<-genus.rel[,-c(which(read.counts<1000)+6)]
genus.rel <- genus.rel[,-which(grepl("clesn",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("215b",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("r240",colnames(genus.rel)))]
colnames(genus.rel)[which(colnames(genus.rel)%in%"s215a")] <- "s215"

genus.rel.root <- genus.rel[,grepl("r\\d",colnames(genus.rel))]
in.root <- as.numeric(gsub("r","",colnames(genus.rel.root)))

genus.rel.soil <- genus.rel[,grepl("s\\d",colnames(genus.rel))]
in.soil <- as.numeric(gsub("s","",colnames(genus.rel.soil)))

### Endosphere sequences

load("Submission_data/ITS_endo_submission.out")
pheno <- plot.frame$Phenotype
max.snps <- plot.frame$max_snps
chr <- plot.frame$Chromosome
pos <- plot.frame$Position
pvalues <- plot.frame$Pval

pheno <- str_replace(pheno," ","_")

seq_ITS_root <- sapply(1:length(chr),function(x) get_interval_sequence(pos[x],chr[x],0.1,pheno[x],"endosphere",pvalues[x]))
seq_ITS_root <- unlist(seq_ITS_root)
seq.names.ITS.root <- names(seq_ITS_root)
seq_ITS_root <- as.list(seq_ITS_root)
seq_ITS_root <- unname(seq_ITS_root)

### Rhizosphere sequences

load("Submission_data/ITS_rhizo_submission.out")
pheno <- plot.frame$Phenotype
max.snps <- plot.frame$max_snps
chr <- plot.frame$Chromosome
pos <- plot.frame$Position
pvalues <- plot.frame$Pval

pheno <- str_replace(pheno," ","_")

seq_ITS_soil <- sapply(1:length(chr),function(x) get_interval_sequence(pos[x],chr[x],0.1,pheno[x],"rhizosphere",pvalues[x]))
seq_ITS_soil <- unlist(seq_ITS_soil)
seq.names.ITS.soil <- names(seq_ITS_soil)
seq_ITS_soil <- as.list(seq_ITS_soil)
seq_ITS_soil <- unname(seq_ITS_soil)

seq_ITS <- c(seq_ITS_root,seq_ITS_soil)
seq.names.ITS <- c(seq.names.ITS.root,seq.names.ITS.soil)

bp.tot <- lapply(seq_ITS,function(x) nchar(x))

write.fasta(sequences = seq_ITS,names = seq.names.ITS,file.out = "~/Documents/Lolium/LoliumSNPs/Submission_data/ITS_seqs_final.fasta",open = "w")

##### 24-06-2026
##### conda create -n blast-env blast python=3
##### source activate blast-env
##### conda install -c conda-forge ncbi-datasets-cli
##### datasets download genome taxon "Arabidopsis thaliana" --include protein --filename arabidopsis_proteome.zip
##### unzip arabidopsis_proteome.zip -d arabidopsis_data
##### makeblastdb -in arabidopsis_data/ncbi_dataset/data/GCA_000001735.2/protein.faa -dbtype prot -out ath_protein_db
##### blastx -query ITS_seqs.fasta -db ath_protein_db -out ITS_blast_results.fmt6 -evalue 1e-5 -num_threads 15 -outfmt 6

############## Here, we make a fasta containing all sequences ##########################

length_chr <- c(length(assembly[[1]]),length(assembly[[2]]),length(assembly[[3]]),length(assembly[[4]]),length(assembly[[5]]),length(assembly[[6]]),
                length(assembly[[7]]))

chr_names <- c("Lp_chr1_0","Lp_chr2_0","Lp_chr3_0","Lp_chr4_0","Lp_chr5_0","Lp_chr6_0","Lp_chr7_0_1")
chr <- c(1,2,3,4,5,6,7)
names(chr) <- c("Lp_chr1_0","Lp_chr2_0","Lp_chr3_0","Lp_chr4_0","Lp_chr5_0","Lp_chr6_0","Lp_chr7_0_1")

gff.gene <- gff[gff$Type%in%c("gene","ncRNA_gene","ncRNA"),]

seq <- list()

for(i in 1:nrow(gff.gene)){
  chro <- gff$Chromosome[i]
  start <- gff$Start[i]
  stop <- gff$Stop[i]
  seq[[i]] <- as.character(assembly[[chro]][start:stop])
}

seq.names <- paste(gff.gene$Chromosome,gff.gene$Start,gff.gene$Stop,gff.gene$Type,gff.gene$V9,sep = "|")

write.fasta(sequences = seq,names = seq.names,file.out = "All_seqs.fasta",open = "w")

##### blastx -query All_seqs.fasta -db ath_protein_db -out all_sequences.fmt6 -evalue 1e-5 -num_threads 15 -outfmt 6



