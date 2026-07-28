#### this script makes genetic correlation matrices from the SNP matrix to be used in the mantel
#### tests

library(ggplot2)

#### Read in the IDs of the L. perenne accessions to be used.
load("~/Documents/Lolium/LoliumSNPs/acc.16s.endo.out")
load("~/Documents/Lolium/LoliumSNPs/acc.16s.rhizo.out")

acc.16s <- union(in.root,in.soil)

load("~/Documents/Lolium/LoliumSNPs/acc.its.endo.out")
load("~/Documents/Lolium/LoliumSNPs/acc.its.rhizo.out")

acc.its <- union(in.root,in.soil)

accs <- union(acc.16s,acc.its)
accs <- accs[order(accs)]

### Read in the SNP matrix
load("~/Documents/Lolium/LoliumSNPs/ubec.out")

snp.info <- ubec[,1:10]

ubec <- data.matrix(ubec[,11:249])
ubec.sub <- ubec[,accs]
rm(ubec)
gc()

### Filter SNPs on MAF of 0.1
threshold <- round(ncol(ubec.sub)*0.1, digits=0)

selc <- apply(ubec.sub==1,1,sum) >= threshold & apply(ubec.sub==3,1,sum) >= threshold

ubec.sub <- ubec.sub[selc,]
snp.info <-snp.info[selc,]

### Code missing genotypes as alternative
ubec.sub[ubec.sub==9]<-3

### Calculate correlations
cor.ubec <- cor(ubec.sub)

save(cor.ubec,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.mat.out")

chr1 <- ubec.sub[snp.info$Chr==1,]
chr1.cor <- cor(chr1)
save(chr1.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr1.out")

chr2 <- ubec.sub[snp.info$Chr==2,]
chr2.cor <- cor(chr2)
save(chr2.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr2.out")

chr3 <- ubec.sub[snp.info$Chr==3,]
chr3.cor <- cor(chr3)
save(chr3.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr3.out")

chr4 <- ubec.sub[snp.info$Chr==4,]
chr4.cor <- cor(chr4)
save(chr4.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr4.out")

chr5 <- ubec.sub[snp.info$Chr==5,]
chr5.cor <- cor(chr5)
save(chr5.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr5.out")

chr6 <- ubec.sub[snp.info$Chr==6,]
chr6.cor <- cor(chr6)
save(chr6.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr6.out")

chr7 <- ubec.sub[snp.info$Chr==7,]
chr7.cor <- cor(chr7)
save(chr7.cor,file = "~/Documents/Lolium/LoliumSNPs/Submission_data/cor.chr7.out")

