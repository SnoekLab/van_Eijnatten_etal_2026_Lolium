### This script calls QTLs for the 16S GWAS results (including genus-rank relative abundances,
### beta diversity, and alpha diversity)

library(ggplot2)
library(stringr)
library(withr)
library(viridis)
library(pals)

### Read in the genus-rank relative abundances and select all valid samples/genotypes

load("genus_rel_filt_mitchlo.out")

OTU.gwas<-read.table("otutable_fungigwasfinal.txt")
read.counts <- apply(OTU.gwas,2,sum)

genus.rel <- genus.rel[genus.rel$genus!= "unidentified",]
genus.rel<-genus.rel[,-c(which(read.counts<1000)+6)]
genus.rel <- genus.rel[,-which(grepl("clesn",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("215b",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("r240",colnames(genus.rel)))]
colnames(genus.rel)[which(colnames(genus.rel)%in%"s215a")] <- "s215"

genus.rel.root <- genus.rel[,grepl("r\\d",colnames(genus.rel))]
in.root <- as.numeric(gsub("r","",colnames(genus.rel.root)))

genus.rel.soil <- genus.rel[,grepl("s\\d",colnames(genus.rel))]
in.soil <- as.numeric(gsub("s","",colnames(genus.rel.soil)))

### Process SNP matrix

load("snp.info.full.out")

max.snps <- c()
pheno <- c()
pvalues <- c()

### This loop calls collects QTLs with -log10(P) > 7 for the rhizosphere traits. It does so by cutting
### the genome into ~1Mbp intervals and collecting the SNP with the highest p-value from each interval.
for(i in c(list.files("GWAS_lolium_genus_its_corgwas_corrected_v2/rhizo/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_shannon_fungi_corrected_v2/soil_depth/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_pcoa_ITS_corrected_v3/soil/frames",full.names = T))){
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
save(plot.frame,file = "ITS_rhizo_submission.out")


### Now do the same for the endosphere phenotypes

max.snps <- c()
pheno <- c()
pvalues <- c()

for(i in c(list.files("GWAS_lolium_genus_its_corgwas_corrected_v2/endo/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_shannon_fungi_corrected_v2/root_depth/frames",full.names = T),list.files("~/NOBINFBACKUP/Lolium/rerun_submission/GWAS_lolium_pcoa_ITS_corrected_v3/root/frames",full.names = T))){
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

plot.frame <- data.frame(Phenotype = pheno,Chromosome = factor(chr,levels = c("1","2","3","4","5","6","7")),Position = pos/10^6,Pval=pvalues,max_snps = max.snps)
save(plot.frame,file = "ITS_endo_submission.out")

facet_bounds <- data.frame(chr = factor(c(1,2,3,4,5,6,7),levels = c("1","2","3","4","5","6","7")),
                           xmin = c(0,0,0,0,0,0,0),
                           xmax= c(max(snp.info$Pos[snp.info$Chr==1]),max(snp.info$Pos[snp.info$Chr==2]),max(snp.info$Pos[snp.info$Chr==3]),
                                   max(snp.info$Pos[snp.info$Chr==4]),max(snp.info$Pos[snp.info$Chr==5]),max(snp.info$Pos[snp.info$Chr==6]),
                                   max(snp.info$Pos[snp.info$Chr==7]))/10^6)

facet_hist <- facet_bounds

ff.hist <- with(facet_hist,data.frame(Position=c(xmin,xmax),Chromosome=c(chr,chr)))

peak.plot <- ggplot(plot.frame,aes(x = Position,y = as.factor(Phenotype),fill = as.factor(Chromosome)))+
  geom_point(aes(size = Pval),shape = 24) + 
  facet_grid(~as.factor(Chromosome),scales = "free",space = "free",drop = F)+
  scale_y_discrete(breaks = NULL) +
  geom_point(data=ff.hist,y=NA)+
  scale_x_continuous(expand = c(0.05,0.05),breaks = c(0,100,200,300)) +
  scale_fill_discrete(type = as.vector(polychrome(7+1)[-2])) +
  theme_minimal() +
  xlab("Position in Mb")+
  ylab("")+
  guides(fill=guide_legend(title="Chromosome",override.aes = list(shape=24))) +
  theme(axis.title.x =element_text(size = 25),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 20),panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size = 25),legend.title = element_text(size = 25),strip.text.x = element_blank(),panel.border = element_rect(fill = NA,color = "black"))

png(filename = "manhat_manhat_Endosphere_point_all_its.png",width = 1400,height = 1000)
peak.plot
dev.off()

peak.plot <- ggplot(plot.frame,aes(x = Position,fill = as.factor(Chromosome)))+
  geom_histogram(data=ff.hist,y=NA)+
  geom_histogram(size = 1,bins = 50) + 
  facet_grid(~as.factor(Chromosome),scales = "free_x",space = "free")+
  geom_histogram(data=ff.hist,y=NA)+
  scale_x_continuous(expand = c(0.05,0.05),breaks = c(0,100,200,300)) +
  scale_fill_discrete(type = as.vector(polychrome(7+1)[-2])) +
  theme_minimal() +
  xlab("Position in Mb")+
  guides(fill=guide_legend(title="Chromosome")) +
  theme(axis.title.x =element_text(size = 25),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 20),panel.spacing.x = unit(1, "lines"),
        legend.text = element_text(size = 25),legend.title = element_text(size = 25),strip.text.x = element_text(size = 25),panel.border = element_rect(fill = NA,color = "black"),
        strip.text.y = element_blank())

png(filename = "manhat_manhat_Endosphere_hist_all_its.png",width = 1400,height = 1000)
peak.plot
dev.off()




