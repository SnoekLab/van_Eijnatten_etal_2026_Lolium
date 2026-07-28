### This script does a PCoA (separately for each kingdom/compartment) to be used for mapping with
### GWAS. Figure S2 is also made in this figure

library(ggplot2)
library(cowplot)

########################################### 16S ################################################

### Read in the bray-curtis distances of the 16S samples
load("bray.16s.out")

### Separate the samples from each compartment
endo <- bray.mat[grepl("GWr",rownames(bray.mat)),grepl("GWr",colnames(bray.mat))]
rhizo <- bray.mat[grepl("GWs",rownames(bray.mat)),grepl("GWs",colnames(bray.mat))]

### PCA on endosphere distances
pca <- prcomp(t(endo))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
pc3 <- pca$rotation[,3]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]
var.pc3 <- pc.sum$importance[2,3]

pca.frame.endo <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pc3,Genotype = rownames(endo))

pca.plot.endo <- ggplot(pca.frame.endo, aes(x = PC1, y = PC2)) +
  geom_point(size = 4.5, alpha = 0.9, stroke = 0.4, shape = 21, col = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.05, 0.05), 
        legend.justification = c(0, 0),         
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)) +
  xlab(paste0("PC1 (", round(var.pc1 * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))

png(filename = "Submission_figures/pca_fig_bac_endo.png",width =400,height = 400)
pca.plot.endo
dev.off()

save(pca.frame.endo, file = "pca.bac.endo.out")

### PCA on rhizosphere distances

pca <- prcomp(t(rhizo))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
pc3 <- pca$rotation[,3]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]
var.pc3 <- pc.sum$importance[2,3]

pca.frame.rhizo <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pc3,Genotype = rownames(rhizo))

pca.plot.rhizo <- ggplot(pca.frame.rhizo, aes(x = PC1, y = PC2)) +
  geom_point(size = 4.5, alpha = 0.9, stroke = 0.4, shape = 21, col = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.05, 0.05), 
        legend.justification = c(0, 0),         
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)) +
  xlab(paste0("PC1 (", round(var.pc1 * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))

png(filename = "pca_fig_bac_rhizo.png",width =400,height = 400)
pca.plot.rhizo
dev.off()

### Load in the SNP matrix

load("ubec.out")

snp.info <- ubec[,1:10]
ubec <- ubec[,11:249]

### Select the peak SNP (all traits mapping to the pleiotropic locus have the same peak SNP).
pleio.pos <- 290459937
pleio.chr <- 7
  
pleio.geno <- ubec[snp.info$Pos==pleio.pos&snp.info$Chr==7,]  

save(pca.frame.rhizo, file = "pca.bac.rhizo.out")
save(pleio.geno,file = "pleiotropic_locus_peak_snp_genotypes.out")

pca.frame.rhizo$Genotype[which(pca.frame.rhizo$Genotype=="GWs_215a")] <- "GWs_215"
gen.nr <- as.numeric(as.vector(sapply(pca.frame.rhizo$Genotype,function(x)strsplit(x,"_")[[1]][2])))

pca.frame.rhizo$Pleio_geno <- as.numeric(pleio.geno[match(gen.nr,names(pleio.geno))])

### Plot the SNP genotype with color on the projection of each accession on the first two 
### principal components (Figure S2)

pca.frame.rhizo$Pleio_geno[pca.frame.rhizo$Pleio_geno==9] <- 3
pca.frame.rhizo$Pleio_geno[pca.frame.rhizo$Pleio_geno==1] <- "REF"
pca.frame.rhizo$Pleio_geno[pca.frame.rhizo$Pleio_geno==2] <- "HET"
pca.frame.rhizo$Pleio_geno[pca.frame.rhizo$Pleio_geno==3] <- "ALT"

pca.plot.rhizo.pleio.alt <- ggplot(pca.frame.rhizo, aes(x = PC1, y = PC2,fill = factor(Pleio_geno,levels = c("REF","HET","ALT")))) +
  geom_point(size = 4.5, alpha = 0.9, stroke = 0.4, shape = 21, col = "black") +
  scale_fill_manual(values = c("REF" = "#D55E00","HET" = "#0072B2","ALT" = "#009E73")) + 
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.02, 0.02), 
        legend.justification = c(0, 0),         
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  xlab(paste0("PC1 (", round(var.pc1 * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))

png(filename = "pca_fig_bac_rhizo_pleio_alt.png",width =425,height = 400)
pca.plot.rhizo.pleio.alt
dev.off()

########################################### ITS ################################################

load("bray.its.out")

endo <- bray.mat[grepl("r",rownames(bray.mat)),grepl("r",colnames(bray.mat))]
rhizo <- bray.mat[grepl("s",rownames(bray.mat)),grepl("s",colnames(bray.mat))]

pca <- prcomp(t(endo))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
pc3 <- pca$rotation[,3]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]
var.pc3 <- pc.sum$importance[2,3]

pca.frame.endo <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pc3,Genotype = rownames(endo))

pca.plot.endo <- ggplot(pca.frame.endo, aes(x = PC1, y = PC2)) +
  geom_point(size = 4.5, alpha = 0.9, stroke = 0.4, shape = 21, col = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.05, 0.05), 
        legend.justification = c(0, 0),         
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)) +
  xlab(paste0("PC1 (", round(var.pc1 * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))

png(filename = "pca_fig_fun_endo.png",width =400,height = 400)
pca.plot.endo
dev.off()

save(pca.frame.endo, file = "pca.fun.endo.out")

pca <- prcomp(t(rhizo))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
pc3 <- pca$rotation[,3]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]
var.pc3 <- pc.sum$importance[2,3]

pca.frame.rhizo <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pc3,Genotype = rownames(rhizo))

pca.plot.rhizo <- ggplot(pca.frame.rhizo, aes(x = PC1, y = PC2)) +
  geom_point(size = 4.5, alpha = 0.9, stroke = 0.4, shape = 21, col = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  #scale_fill_manual(values = c("#081D58","#D73027")) + 
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.05, 0.05), 
        legend.justification = c(0, 0),         
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)) +
  xlab(paste0("PC1 (", round(var.pc1 * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))

png(filename = "pca_fig_fun_rhizo.png",width =400,height = 400)
pca.plot.rhizo
dev.off()

save(pca.frame.rhizo, file = "pca.fun.rhizo.out")






