#### This script does mantel tests for correlation between plant genetic distance and plant
#### Microbiome distance

library(vegan)

### Read in bray curtis distance based on 16S microbiome and correlation matrix based on 
### the SNP matrix
load("bray.16s.out")
load("cor.mat.out")

### Separate endosphere and rhizosphere
endo <- bray.mat[grepl("r",rownames(bray.mat),fixed = T),grepl("r",colnames(bray.mat),fixed = T)]
rhizo <- bray.mat[grepl("s",rownames(bray.mat),fixed = T),grepl("s",colnames(bray.mat),fixed = T)]
colnames(rhizo)[which(colnames(rhizo)=="GWs_215a")] <- "GWs_215"

### Make sure the genotypes match and are in the same order
endo.acc <- as.numeric(as.vector(sapply(colnames(endo),function(x)strsplit(x,"_",fixed = T)[[1]][2])))
rhizo.acc <- as.numeric(as.vector(sapply(colnames(rhizo),function(x)strsplit(x,"_",fixed = T)[[1]][2])))

endo.cor <- cor.ubec[rownames(cor.ubec)%in%endo.acc,colnames(cor.ubec)%in%endo.acc]
colnames(endo.cor) <- colnames(endo)
rownames(endo.cor) <- rownames(endo)

rhizo.cor <- cor.ubec[rownames(cor.ubec)%in%rhizo.acc,colnames(cor.ubec)%in%rhizo.acc]
colnames(rhizo.cor) <- colnames(rhizo)
rownames(rhizo.cor) <- rownames(rhizo)

### Turn into distance matrix
endo <- as.dist(as.matrix(endo))
rhizo <- as.dist(as.matrix(rhizo))

endo.dist <- as.dist(as.matrix(1 - endo.cor))
rhizo.dist <- as.dist(as.matrix(1 - rhizo.cor))

### Do mantel test
mantel.16s.endo <- mantel(endo, endo.dist, method = "spearman", permutations = 9999)
mantel.16s.rhizo <- mantel(rhizo, rhizo.dist, method = "spearman", permutations = 9999)

### Do the same tests, but now per chromosome
load("cor.chr1.out")
load("cor.chr2.out")
load("cor.chr3.out")
load("cor.chr4.out")
load("cor.chr5.out")
load("cor.chr6.out")
load("cor.chr7.out")

endo.chr1 <- chr1.cor[rownames(chr1.cor)%in%endo.acc,colnames(chr1.cor)%in%endo.acc]
colnames(endo.chr1) <- colnames(endo)
rownames(endo.chr1) <- rownames(endo)

rhizo.chr1 <- chr1.cor[rownames(chr1.cor)%in%rhizo.acc,colnames(chr1.cor)%in%rhizo.acc]
colnames(rhizo.chr1) <- colnames(rhizo)
rownames(rhizo.chr1) <- rownames(rhizo)

endo.chr1<- as.dist(as.matrix(1 - endo.chr1))
rhizo.chr1 <- as.dist(as.matrix(1 - rhizo.chr1))

endo.chr2 <- chr2.cor[rownames(chr2.cor)%in%endo.acc,colnames(chr2.cor)%in%endo.acc]
colnames(endo.chr2) <- colnames(endo)
rownames(endo.chr2) <- rownames(endo)

rhizo.chr2 <- chr2.cor[rownames(chr2.cor)%in%rhizo.acc,colnames(chr2.cor)%in%rhizo.acc]
colnames(rhizo.chr2) <- colnames(rhizo)
rownames(rhizo.chr2) <- rownames(rhizo)

endo.chr2<- as.dist(as.matrix(1 - endo.chr2))
rhizo.chr2 <- as.dist(as.matrix(1 - rhizo.chr2))

endo.chr3 <- chr3.cor[rownames(chr3.cor)%in%endo.acc,colnames(chr3.cor)%in%endo.acc]
colnames(endo.chr3) <- colnames(endo)
rownames(endo.chr3) <- rownames(endo)

rhizo.chr3 <- chr3.cor[rownames(chr3.cor)%in%rhizo.acc,colnames(chr3.cor)%in%rhizo.acc]
colnames(rhizo.chr3) <- colnames(rhizo)
rownames(rhizo.chr3) <- rownames(rhizo)

endo.chr3<- as.dist(as.matrix(1 - endo.chr3))
rhizo.chr3 <- as.dist(as.matrix(1 - rhizo.chr3))

endo.chr4 <- chr4.cor[rownames(chr4.cor)%in%endo.acc,colnames(chr4.cor)%in%endo.acc]
colnames(endo.chr4) <- colnames(endo)
rownames(endo.chr4) <- rownames(endo)

rhizo.chr4 <- chr4.cor[rownames(chr4.cor)%in%rhizo.acc,colnames(chr4.cor)%in%rhizo.acc]
colnames(rhizo.chr4) <- colnames(rhizo)
rownames(rhizo.chr4) <- rownames(rhizo)

endo.chr4<- as.dist(as.matrix(1 - endo.chr4))
rhizo.chr4 <- as.dist(as.matrix(1 - rhizo.chr4))

endo.chr5 <- chr5.cor[rownames(chr5.cor)%in%endo.acc,colnames(chr5.cor)%in%endo.acc]
colnames(endo.chr5) <- colnames(endo)
rownames(endo.chr5) <- rownames(endo)

rhizo.chr5 <- chr5.cor[rownames(chr5.cor)%in%rhizo.acc,colnames(chr5.cor)%in%rhizo.acc]
colnames(rhizo.chr5) <- colnames(rhizo)
rownames(rhizo.chr5) <- rownames(rhizo)

endo.chr5<- as.dist(as.matrix(1 - endo.chr5))
rhizo.chr5 <- as.dist(as.matrix(1 - rhizo.chr5))

endo.chr6 <- chr6.cor[rownames(chr6.cor)%in%endo.acc,colnames(chr6.cor)%in%endo.acc]
colnames(endo.chr6) <- colnames(endo)
rownames(endo.chr6) <- rownames(endo)

rhizo.chr6 <- chr6.cor[rownames(chr6.cor)%in%rhizo.acc,colnames(chr6.cor)%in%rhizo.acc]
colnames(rhizo.chr6) <- colnames(rhizo)
rownames(rhizo.chr6) <- rownames(rhizo)

endo.chr6<- as.dist(as.matrix(1 - endo.chr6))
rhizo.chr6 <- as.dist(as.matrix(1 - rhizo.chr6))

endo.chr7 <- chr7.cor[rownames(chr7.cor)%in%endo.acc,colnames(chr7.cor)%in%endo.acc]
colnames(endo.chr7) <- colnames(endo)
rownames(endo.chr7) <- rownames(endo)

rhizo.chr7 <- chr7.cor[rownames(chr7.cor)%in%rhizo.acc,colnames(chr7.cor)%in%rhizo.acc]
colnames(rhizo.chr7) <- colnames(rhizo)
rownames(rhizo.chr7) <- rownames(rhizo)

endo.chr7<- as.dist(as.matrix(1 - endo.chr7))
rhizo.chr7 <- as.dist(as.matrix(1 - rhizo.chr7))

mantel.16s.endo.chr1 <- mantel(endo, endo.chr1, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr1 <- mantel(rhizo, rhizo.chr1, method = "spearman", permutations = 9999)

mantel.16s.endo.chr2 <- mantel(endo, endo.chr2, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr2 <- mantel(rhizo, rhizo.chr2, method = "spearman", permutations = 9999)

mantel.16s.endo.chr3 <- mantel(endo, endo.chr3, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr3 <- mantel(rhizo, rhizo.chr3, method = "spearman", permutations = 9999)

mantel.16s.endo.chr4 <- mantel(endo, endo.chr4, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr4 <- mantel(rhizo, rhizo.chr4, method = "spearman", permutations = 9999)

mantel.16s.endo.chr5 <- mantel(endo, endo.chr5, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr5 <- mantel(rhizo, rhizo.chr5, method = "spearman", permutations = 9999)

mantel.16s.endo.chr6 <- mantel(endo, endo.chr6, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr6 <- mantel(rhizo, rhizo.chr6, method = "spearman", permutations = 9999)

mantel.16s.endo.chr7 <- mantel(endo, endo.chr7, method = "spearman", permutations = 9999)
mantel.16s.rhizo.chr7 <- mantel(rhizo, rhizo.chr7, method = "spearman", permutations = 9999)


### Repeat for the ITS data
load("bray.its.out")
load("cor.mat.out")

endo <- bray.mat[grepl("r",rownames(bray.mat),fixed = T),grepl("r",colnames(bray.mat),fixed = T)]
rhizo <- bray.mat[grepl("s",rownames(bray.mat),fixed = T),grepl("s",colnames(bray.mat),fixed = T)]

endo.acc <- as.numeric(as.vector(sapply(colnames(endo),function(x)strsplit(x,"r",fixed = T)[[1]][2])))
rhizo.acc <- as.numeric(as.vector(sapply(colnames(rhizo),function(x)strsplit(x,"s",fixed = T)[[1]][2])))

endo <- endo[order(endo.acc),order(endo.acc)]
rhizo <- rhizo[order(rhizo.acc),order(rhizo.acc)]

endo.cor <- cor.ubec[rownames(cor.ubec)%in%endo.acc,colnames(cor.ubec)%in%endo.acc]
colnames(endo.cor) <- colnames(endo)
rownames(endo.cor) <- rownames(endo)

rhizo.cor <- cor.ubec[rownames(cor.ubec)%in%rhizo.acc,colnames(cor.ubec)%in%rhizo.acc]
colnames(rhizo.cor) <- colnames(rhizo)
rownames(rhizo.cor) <- rownames(rhizo)

endo <- as.dist(as.matrix(endo))
rhizo <- as.dist(as.matrix(rhizo))

endo.dist <- as.dist(as.matrix(1 - endo.cor))
rhizo.dist <- as.dist(as.matrix(1 - rhizo.cor))

mantel.its.endo <- mantel(endo, endo.dist, method = "spearman", permutations = 9999)
mantel.its.rhizo <- mantel(rhizo, rhizo.dist, method = "spearman", permutations = 9999)

####################################### ITS per chromosome #######################################

load("cor.chr1.out")
load("cor.chr2.out")
load("cor.chr3.out")
load("cor.chr4.out")
load("cor.chr5.out")
load("cor.chr6.out")
load("cor.chr7.out")

endo.chr1 <- chr1.cor[rownames(chr1.cor)%in%endo.acc,colnames(chr1.cor)%in%endo.acc]
colnames(endo.chr1) <- colnames(endo)
rownames(endo.chr1) <- rownames(endo)

rhizo.chr1 <- chr1.cor[rownames(chr1.cor)%in%rhizo.acc,colnames(chr1.cor)%in%rhizo.acc]
colnames(rhizo.chr1) <- colnames(rhizo)
rownames(rhizo.chr1) <- rownames(rhizo)

endo.chr1<- as.dist(as.matrix(1 - endo.chr1))
rhizo.chr1 <- as.dist(as.matrix(1 - rhizo.chr1))

endo.chr2 <- chr2.cor[rownames(chr2.cor)%in%endo.acc,colnames(chr2.cor)%in%endo.acc]
colnames(endo.chr2) <- colnames(endo)
rownames(endo.chr2) <- rownames(endo)

rhizo.chr2 <- chr2.cor[rownames(chr2.cor)%in%rhizo.acc,colnames(chr2.cor)%in%rhizo.acc]
colnames(rhizo.chr2) <- colnames(rhizo)
rownames(rhizo.chr2) <- rownames(rhizo)

endo.chr2<- as.dist(as.matrix(1 - endo.chr2))
rhizo.chr2 <- as.dist(as.matrix(1 - rhizo.chr2))

endo.chr3 <- chr3.cor[rownames(chr3.cor)%in%endo.acc,colnames(chr3.cor)%in%endo.acc]
colnames(endo.chr3) <- colnames(endo)
rownames(endo.chr3) <- rownames(endo)

rhizo.chr3 <- chr3.cor[rownames(chr3.cor)%in%rhizo.acc,colnames(chr3.cor)%in%rhizo.acc]
colnames(rhizo.chr3) <- colnames(rhizo)
rownames(rhizo.chr3) <- rownames(rhizo)

endo.chr3<- as.dist(as.matrix(1 - endo.chr3))
rhizo.chr3 <- as.dist(as.matrix(1 - rhizo.chr3))

endo.chr4 <- chr4.cor[rownames(chr4.cor)%in%endo.acc,colnames(chr4.cor)%in%endo.acc]
colnames(endo.chr4) <- colnames(endo)
rownames(endo.chr4) <- rownames(endo)

rhizo.chr4 <- chr4.cor[rownames(chr4.cor)%in%rhizo.acc,colnames(chr4.cor)%in%rhizo.acc]
colnames(rhizo.chr4) <- colnames(rhizo)
rownames(rhizo.chr4) <- rownames(rhizo)

endo.chr4<- as.dist(as.matrix(1 - endo.chr4))
rhizo.chr4 <- as.dist(as.matrix(1 - rhizo.chr4))

endo.chr5 <- chr5.cor[rownames(chr5.cor)%in%endo.acc,colnames(chr5.cor)%in%endo.acc]
colnames(endo.chr5) <- colnames(endo)
rownames(endo.chr5) <- rownames(endo)

rhizo.chr5 <- chr5.cor[rownames(chr5.cor)%in%rhizo.acc,colnames(chr5.cor)%in%rhizo.acc]
colnames(rhizo.chr5) <- colnames(rhizo)
rownames(rhizo.chr5) <- rownames(rhizo)

endo.chr5<- as.dist(as.matrix(1 - endo.chr5))
rhizo.chr5 <- as.dist(as.matrix(1 - rhizo.chr5))

endo.chr6 <- chr6.cor[rownames(chr6.cor)%in%endo.acc,colnames(chr6.cor)%in%endo.acc]
colnames(endo.chr6) <- colnames(endo)
rownames(endo.chr6) <- rownames(endo)

rhizo.chr6 <- chr6.cor[rownames(chr6.cor)%in%rhizo.acc,colnames(chr6.cor)%in%rhizo.acc]
colnames(rhizo.chr6) <- colnames(rhizo)
rownames(rhizo.chr6) <- rownames(rhizo)

endo.chr6<- as.dist(as.matrix(1 - endo.chr6))
rhizo.chr6 <- as.dist(as.matrix(1 - rhizo.chr6))

endo.chr7 <- chr7.cor[rownames(chr7.cor)%in%endo.acc,colnames(chr7.cor)%in%endo.acc]
colnames(endo.chr7) <- colnames(endo)
rownames(endo.chr7) <- rownames(endo)

rhizo.chr7 <- chr7.cor[rownames(chr7.cor)%in%rhizo.acc,colnames(chr7.cor)%in%rhizo.acc]
colnames(rhizo.chr7) <- colnames(rhizo)
rownames(rhizo.chr7) <- rownames(rhizo)

endo.chr7<- as.dist(as.matrix(1 - endo.chr7))
rhizo.chr7 <- as.dist(as.matrix(1 - rhizo.chr7))

mantel.its.endo.chr1 <- mantel(endo, endo.chr1, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr1 <- mantel(rhizo, rhizo.chr1, method = "spearman", permutations = 9999)

mantel.its.endo.chr2 <- mantel(endo, endo.chr2, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr2 <- mantel(rhizo, rhizo.chr2, method = "spearman", permutations = 9999)

mantel.its.endo.chr3 <- mantel(endo, endo.chr3, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr3 <- mantel(rhizo, rhizo.chr3, method = "spearman", permutations = 9999)

mantel.its.endo.chr4 <- mantel(endo, endo.chr4, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr4 <- mantel(rhizo, rhizo.chr4, method = "spearman", permutations = 9999)

mantel.its.endo.chr5 <- mantel(endo, endo.chr5, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr5 <- mantel(rhizo, rhizo.chr5, method = "spearman", permutations = 9999)### significance 0.0553

mantel.its.endo.chr6 <- mantel(endo, endo.chr6, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr6 <- mantel(rhizo, rhizo.chr6, method = "spearman", permutations = 9999)

mantel.its.endo.chr7 <- mantel(endo, endo.chr7, method = "spearman", permutations = 9999)
mantel.its.rhizo.chr7 <- mantel(rhizo, rhizo.chr7, method = "spearman", permutations = 9999)


#### Store all results in a dataframe
mantel.frame <- data.frame(P_value = c(mantel.16s.endo$signif,mantel.16s.rhizo$signif,
                           mantel.its.endo$signif,mantel.its.rhizo$signif,
                           mantel.16s.endo.chr1$signif,mantel.16s.rhizo.chr1$signif,
                           mantel.its.endo.chr1$signif,mantel.its.rhizo.chr1$signif,
                           mantel.16s.endo.chr2$signif,mantel.16s.rhizo.chr2$signif,
                           mantel.its.endo.chr2$signif,mantel.its.rhizo.chr2$signif,
                           mantel.16s.endo.chr3$signif,mantel.16s.rhizo.chr3$signif,
                           mantel.its.endo.chr3$signif,mantel.its.rhizo.chr3$signif,
                           mantel.16s.endo.chr4$signif,mantel.16s.rhizo.chr4$signif,
                           mantel.its.endo.chr4$signif,mantel.its.rhizo.chr4$signif,
                           mantel.16s.endo.chr5$signif,mantel.16s.rhizo.chr5$signif,
                           mantel.its.endo.chr5$signif,mantel.its.rhizo.chr5$signif,
                           mantel.16s.endo.chr6$signif,mantel.16s.rhizo.chr6$signif,
                           mantel.its.endo.chr6$signif,mantel.its.rhizo.chr6$signif,
                           mantel.16s.endo.chr7$signif,mantel.16s.rhizo.chr7$signif,
                           mantel.its.endo.chr7$signif,mantel.its.rhizo.chr7$signif),
                           R_statistic = c(mantel.16s.endo$statistic,mantel.16s.rhizo$statistic,
                                           mantel.its.endo$statistic,mantel.its.rhizo$statistic,
                                           mantel.16s.endo.chr1$statistic,mantel.16s.rhizo.chr1$statistic,
                                           mantel.its.endo.chr1$statistic,mantel.its.rhizo.chr1$statistic,
                                           mantel.16s.endo.chr2$statistic,mantel.16s.rhizo.chr2$statistic,
                                           mantel.its.endo.chr2$statistic,mantel.its.rhizo.chr2$statistic,
                                           mantel.16s.endo.chr3$statistic,mantel.16s.rhizo.chr3$statistic,
                                           mantel.its.endo.chr3$statistic,mantel.its.rhizo.chr3$statistic,
                                           mantel.16s.endo.chr4$statistic,mantel.16s.rhizo.chr4$statistic,
                                           mantel.its.endo.chr4$statistic,mantel.its.rhizo.chr4$statistic,
                                           mantel.16s.endo.chr5$statistic,mantel.16s.rhizo.chr5$statistic,
                                           mantel.its.endo.chr5$statistic,mantel.its.rhizo.chr5$statistic,
                                           mantel.16s.endo.chr6$statistic,mantel.16s.rhizo.chr6$statistic,
                                           mantel.its.endo.chr6$statistic,mantel.its.rhizo.chr6$statistic,
                                           mantel.16s.endo.chr7$statistic,mantel.16s.rhizo.chr7$statistic,
                                           mantel.its.endo.chr7$statistic,mantel.its.rhizo.chr7$statistic),
                           Compartment = rep(c("Endosphere","Rhizosphere"),16),
                           Community = rep(c("Bacterial","Bacterial","Fungal","Fungal"),8),
                           Chromosome = c(rep("Full genome",4),rep("chr1",4),rep("chr2",4),rep("chr3",4),
                                          rep("chr4",4),rep("chr5",4),rep("chr6",4),rep("chr7",4)))


write.csv(mantel.frame,file = "mantel_results.csv")

