### This script read in the GWAS output of the pleiotropic locus and performs Benjamini-Hochberg correction

### Read in the GWAS results of all traits mapping to the pleiotropic locus
gwas.files <- c(list.files("GWAS_full_pleitropic_locus_corrected_v2/frames",full.names = T))
gwas.files <- gwas.files[-which(grepl("PC2",gwas.files)|grepl("PC3",gwas.files))]

########### Acidothermus relative abundance, P = 0.014
load(gwas.files[1])
acido <- plot_frame
acido$pval_bh <- p.adjust(acido$pval,method = "BH")
acido <- acido[acido$chr==7&acido$pos>28&acido$pos<30,]
ac.bh <- min(acido$pval_bh)
acido[which.min(acido$pval_bh),]
ac.log.pval <- max(-log10(acido$pval))
ac.pval <- min(acido$pval)

########### Gaiella relative abundance, P = 0.080
load(gwas.files[2])
ga <- plot_frame
ga$pval_bh <- p.adjust(ga$pval,method = "BH")
ga <- ga[ga$chr==7&ga$pos>28&ga$pos<30,]
ga.bh <- min(ga$pval_bh)
ga[which.min(ga$pval_bh),]
ga.log.pval <- max(-log10(ga$pval))
ga.pval <- min(ga$pval)

########### mle1-7 relative abundance, P = 0.054
load(gwas.files[3])
ml <- plot_frame
ml$pval_bh <- p.adjust(ml$pval,method = "BH")
ml <- ml[ml$chr==7&ml$pos>28&ml$pos<30,]
ml.bh <- min(ml$pval_bh)
ml[which.min(ml$pval_bh),]
ml.log.pval <- max(-log10(ml$pval))
ml.pval <- min(ml$pval)

########### PC1, P = 0.075
load(gwas.files[4])
pc1 <- plot_frame
pc1$pval_bh <- p.adjust(pc1$pval,method = "BH")
pc1 <- pc1[pc1$chr==7&pc1$pos>28&pc1$pos<30,]
pc1.bh <- min(pc1$pval_bh)
pc1[order(pc1$pval_bh)[1:10],]
pc1.log.pval <- max(-log10(pc1$pval))
pc1.pval <- min(pc1$pval)

########### Pseudolabrys relative abundance, P = 0.695
load(gwas.files[5])
ps <- plot_frame
ps$pval_bh <- p.adjust(ps$pval,method = "BH")
ps <- ps[ps$chr==7&ps$pos>28&ps$pos<30,]
ps.bh <- min(ps$pval_bh)
ps[which.min(ps$pval_bh),]
ps.log.pval <- max(-log10(ps$pval))
ps.pval <- min(ps$pval)

########### Pseudonocardia relative abundance, P = 0.120
load(gwas.files[6])
pc <- plot_frame
pc$pval_bh <- p.adjust(pc$pval,method = "BH")
pc <- pc[pc$chr==7&pc$pos>28&pc$pos<30,]
pc.bh <- min(pc$pval_bh)
pc[which.min(pc$pval_bh),]
pc.log.pval <- max(-log10(pc$pval))
pc.pval <- min(pc$pval)

########### RB41 relative abundance, P = 0.072
load(gwas.files[7])
RB41 <- plot_frame
RB41$pval_bh <- p.adjust(RB41$pval,method = "BH")
RB41 <- RB41[RB41$chr==7&RB41$pos>28&RB41$pos<30,]
RB41.bh <- min(RB41$pval_bh)
RB41[order(RB41$pval_bh)[1:10],]
RB41.log.pval <- max(-log10(RB41$pval))
RB41.pval <- min(RB41$pval)

########### Shannon Depth, P = 0.038
load(gwas.files[8])
shannon <- plot_frame
shannon$pval_bh <- p.adjust(shannon$pval,method = "BH")
shannon <- shannon[shannon$chr==7&shannon$pos>28&shannon$pos<30,]
shannon.bh <- min(shannon$pval_bh)
shannon[order(shannon$pval_bh)[1:3],]
shannon.log.pval <- max(-log10(shannon$pval))
shannon.pval <- min(shannon$pval)

########### Tumebacillus relative abundance, P = 0.016
load(gwas.files[9])
tumebacillus <- plot_frame
tumebacillus$pval_bh <- p.adjust(tumebacillus$pval,method = "BH")
tumebacillus <- tumebacillus[tumebacillus$chr==7&tumebacillus$pos>28&tumebacillus$pos<30,]
tumebacillus.bh <- min(tumebacillus$pval_bh)
tumebacillus[order(tumebacillus$pval_bh)[1:10],]
tumebacillus.log.pval <- max(-log10(tumebacillus$pval))
tumebacillus.pval <- min(tumebacillus$pval)

### Put all results in a dataframe
dat <- data.frame(Phenotype = c("Acidothermus","Gaiella","mle1-7","PC1","Pseudolabrys","Pseudonocardia","RB41",
                         "Shannon depth","Tumebacillus"),
           Pval = c(ac.pval,ga.pval,ml.pval,pc1.pval,ps.pval,pc.pval,RB41.pval,shannon.pval,tumebacillus.pval),
           Log10_pval = c(ac.log.pval,ga.log.pval,ml.log.pval,pc1.log.pval,ps.log.pval,pc.log.pval,RB41.log.pval,shannon.log.pval,tumebacillus.log.pval),
           BH_pval = c(ac.bh,ga.bh,ml.bh,pc1.bh,ps.bh,pc.bh,RB41.bh,shannon.bh,tumebacillus.bh))

write.csv(dat,"pvals_pleiotropic_locus.csv")





