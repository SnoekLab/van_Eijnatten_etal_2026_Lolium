### This script calls QTLs for the 16S GWAS results (including genus-rank relative abundances,
### beta diversity, and alpha diversity)

library(ggplot2)
library(stringr)
library(withr)
library(viridis)
library(pals)

### Read in the genus-rank relative abundances and select all valid samples/genotypes

load("genus_rel_filt_mitchlo.out")

genus.rel <- genus.rel[,-which(colnames(genus.rel)%in%c("GWr_170","GWs_179","GWr_065","GWs_018","GWs_074","GWs_001","GWs_002","GWs_120","GWr_014","GWr_189","GWs_181_clesn","GWs_220_clesn",
                                                        "GWs_64_clesn","GWs_215b","GWs_000"))]
colnames(genus.rel)[which(colnames(genus.rel)=="GWs_215a")] <- "GWs_215"
genus.rel.root <- genus.rel[grep("GWr",colnames(genus.rel))]
in.root <- as.numeric(gsub("GWr_","",colnames(genus.rel.root)))
genus.rel.soil <- genus.rel[,grep("GWs",colnames(genus.rel))]
in.soil <- as.numeric(gsub("GWs_","",colnames(genus.rel.soil)))

### Process SNP matrix
load("snp.info.full.out")

max.snps <- c()
pheno <- c()
pvalues <- c()

### This loop calls collects QTLs with -log10(P) > 7 for the rhizosphere traits. It does so by cutting
### the genome into ~1Mbp intervals and collecting the SNP with the highest p-value from each interval.

for(i in c(list.files("GWAS_lolium_genus_16s_corgwas_corrected_v2/rhizo/frames",full.names = T),list.files("GWAS_lolium_shannon_16S_corrected_v2/soil_depth/frames",full.names = T),list.files("GWAS_lolium_pcoa_16S_corrected_v3/soil/frames",full.names = T))){
  load(i)
  phe <- str_extract(i,"(?<=/frames/)(.*)(?=[.])")
  
  pval <- -log10(plot_frame$pval)
  chr <- as.numeric(plot_frame$chr)
  pos <- as.numeric(plot_frame$pos)
  
  to.pl <- data.frame(chr,pos,pval)
  
  maxpmil <- as.matrix(aggregate(to.pl$pval,list(to.pl$chr,as.numeric(cut(to.pl$pos,400))),max))
  
  
  if(sum(maxpmil[,3]>7)>1){
    mm.snps <- maxpmil[maxpmil[,3]>7,]
    
    snp.nums <- NULL
    
    for (j in 1:nrow(mm.snps)){
      snp.cor <- which(snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),1]==mm.snps[j,1])[1]
      snp.nums[j] <- which(snp.info$Chr == as.character(mm.snps[j,1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),2][snp.cor])
      print(which(snp.info$Chr == as.character(mm.snps[j,1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),2][snp.cor]))
    }
    
    max.snps <- c(max.snps,snp.nums)
    pheno <- c(pheno,rep(phe,sum(maxpmil[,3]>7)))
    pvalues <- c(pvalues,mm.snps[,3])
  }
  
  if(sum(maxpmil[,3]>7)==1){
    mm.snps <- maxpmil[maxpmil[,3]>7,]
    
    snp.nums <- NULL
    
    for (j in c(1)){
      snp.cor <- which(snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[3]]*10^7))),1]==mm.snps[1])[1]
      snp.nums[j] <- which(snp.info$Chr == as.character(mm.snps[1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[3]]*10^7))),2][snp.cor])
    }
    
    max.snps <- c(max.snps,snp.nums)
    pheno <- c(pheno,phe)
    pvalues <- c(pvalues,mm.snps[3])
  }
}

pvalues <- as.vector(pvalues)

pos <- snp.info[max.snps,2]
chr <- snp.info[max.snps,1]

### Save QTLs
plot.frame <- data.frame(Phenotype = pheno,Chromosome = chr,Position = pos/10^6,Pval=pvalues,max_snps = max.snps)
save(plot.frame,file = "16S_rhizo_submission.out")

### Now do the same for the endosphere phenotypes

max.snps <- c()
pheno <- c()
pvalues <- c()

for(i in c(list.files("GWAS_lolium_genus_16s_corgwas_corrected_v2/endo/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_shannon_16S_corrected_v2/root_depth/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_pcoa_16S_corrected_v3/root/frames",full.names = T))){
  load(i)
  phe <- str_extract(i,"(?<=/frames/)(.*)(?=[.])")
  
  pval <- -log10(plot_frame$pval)
  chr <- as.numeric(plot_frame$chr)
  pos <- as.numeric(plot_frame$pos)
  
  to.pl <- data.frame(chr,pos,pval)
  
  maxpmil <- as.matrix(aggregate(to.pl$pval,list(to.pl$chr,as.numeric(cut(to.pl$pos,400))),max))
  
  
  if(sum(maxpmil[,3]>7)>1){
    mm.snps <- maxpmil[maxpmil[,3]>7,]
    
    snp.nums <- NULL
    
    for (j in 1:nrow(mm.snps)){
      snp.cor <- which(snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),1]==mm.snps[j,1])[1]
      snp.nums[j] <- which(snp.info$Chr == as.character(mm.snps[j,1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),2][snp.cor])
      print(which(snp.info$Chr == as.character(mm.snps[j,1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[j,3]]*10^7))),2][snp.cor]))
    }
    
    max.snps <- c(max.snps,snp.nums)
    pheno <- c(pheno,rep(phe,sum(maxpmil[,3]>7)))
    pvalues <- c(pvalues,mm.snps[,3])
  }
  
  if(sum(maxpmil[,3]>7)==1){
    mm.snps <- maxpmil[maxpmil[,3]>7,]
    
    snp.nums <- NULL
    
    for (j in c(1)){
      snp.cor <- which(snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[3]]*10^7))),1]==mm.snps[1])[1]
      snp.nums[j] <- which(snp.info$Chr == as.character(mm.snps[1]) & snp.info$Pos == snp.info[snp.info$Pos%in%c(as.integer(as.character(plot_frame$pos[pval == mm.snps[3]]*10^7))),2][snp.cor])
    }
    
    max.snps <- c(max.snps,snp.nums)
    pheno <- c(pheno,phe)
    pvalues <- c(pvalues,mm.snps[3])
  }
}

pvalues <- as.vector(pvalues)

pos <- snp.info[max.snps,2]
chr <- snp.info[max.snps,1]

plot.frame <- data.frame(Phenotype = pheno,Chromosome = chr,Position = pos/10^6,Pval=pvalues,max_snps = max.snps)
save(plot.frame,file = "16S_endo_submission.out")




