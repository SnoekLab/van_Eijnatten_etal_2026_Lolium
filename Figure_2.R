### This script makes figure 2.

library(ggplot2)

### Read in the 16S genus rank relative abundances.
load("genus_rel_filt_mitchlo.out")

### Remove bad and duplicate samples
genus.rel <- genus.rel[,-which(colnames(genus.rel)%in%c("GWr_170","GWs_179","GWr_065","GWs_018","GWs_074","GWs_001","GWs_002","GWs_120","GWr_014","GWr_189","GWs_181_clesn","GWs_220_clesn",
                                                        "GWs_64_clesn","GWs_215b","GWs_000"))]
tax <- genus.rel[,1:6]
genus.rel <- genus.rel[,-c(1:6)]

### Calculate bray curtis distance
sim<-function(data, index="bc"){
  nplot=dim(data)[2]
  similar=matrix(1,ncol=nplot,nrow=nplot)
  rownames(similar)<-colnames(data)
  colnames(similar)<-colnames(data)
  for(i in 1:(nplot-1)){
    m=i+1
    for(m in m:nplot){
      if(index=="bc"){
        bc.sim=2*sum(apply(data[,c(i,m)], 1, min))/(sum (data[,c(i,m)]))
        similar[i,m]=bc.sim
        similar[m,i]=bc.sim
      }
    }
  }
  return(round(similar,3))
}

bray.mat <- sim(genus.rel)

save(bray.mat,file = "bray.16s.out")

### PCA on bray curtis matrix
pca <- prcomp(t(bray.mat))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]

pca.frame.bac <- data.frame(PC1 = pc1, PC2 = pc2,Habitat = ifelse(grepl("GWr",colnames(bray.mat),fixed = T), "Endosphere", "Rhizosphere"))

pca.plot.bac <- ggplot(pca.frame.bac, aes(x = PC1, y = PC2,fill = factor(Habitat))) +
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
  ylab("PC2 (4.0%)")+ 
  annotate("text", x = -0.07, y = 0.12, label = "Bacteria", 
           hjust = -0.1, vjust = 1.1, size = 6)

### Read in the ASV-level data for calculating Shannon diversity
OTU.gwas<-read.delim("ASVs_counts_gwas.txt")
TAX.gwas<-read.delim("ASVs_taxonomy_gwas.txt")

### Remove all chloroplast and mitochondrial data.
rselc <- TAX.gwas$Order != "Chloroplast"&TAX.gwas$Family!= "Mitochondria"
rselc[is.na(rselc)] <- T

TAX.gwas <- TAX.gwas[rselc,]
OTU.gwas <- OTU.gwas[rselc,]

### Remove bad and duplicate samples
OTU.gwas<-OTU.gwas[,-which(colnames(OTU.gwas)%in%c("GWr_170","GWs_179","GWr_065","GWs_018","GWs_074","GWs_001","GWs_002","GWs_120","GWr_014","GWr_189","GWs_215b",
                                                   "GWs_220_clesn","GWs_181_clesn","GWs_64_clesn","GWs_000"))]

### Calculate read depth of each sample
depth <- apply(OTU.gwas,2,sum)

### Change abundances into proportions.
OTU.gwas <- apply(OTU.gwas,2,function(x) x/sum(x))

### Calculate in how many samples each clade occurs.
occ <- apply(OTU.gwas,1,function(x) sum(x>0))

OTU.gwas <- OTU.gwas[occ>=1,]

### Calculate Shannon diversity
shannon <- apply(OTU.gwas,2,function(x) -sum(x*log(x+0.0001)))

shannon.frame <- data.frame(Shannon = shannon, Habitat = ifelse(grepl("GWr",names(shannon),fixed = T), "Endosphere", "Rhizosphere"))

bac.shannon.box <- ggplot(shannon.frame, aes(x = factor(Habitat), y = Shannon)) +
  geom_violin(aes(fill = Habitat), alpha = 0.6, show.legend = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, fill = NA, color = "black", show.legend = FALSE) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  labs(fill = "Environment", y = "Bacterial diversity")+
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        legend.position = "none")


### Now do the same steps for the fungi

OTU.gwas<-read.table("otutable_fungigwasfinal.txt")
read.counts <- apply(OTU.gwas,2,sum)

load("genus_rel_filt_mitchlo.out")

genus.rel<-genus.rel[,-c(which(read.counts<1000)+6)]
genus.rel <- genus.rel[,-which(grepl("clesn",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("215b",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("r240",colnames(genus.rel)))]
colnames(genus.rel)[which(colnames(genus.rel)%in%"s215a")] <- "s215"

tax <- genus.rel[,1:6]
genus.rel <- genus.rel[,-c(1:6)]

sim<-function(data, index="bc"){
  nplot=dim(data)[2]
  similar=matrix(1,ncol=nplot,nrow=nplot)
  rownames(similar)<-colnames(data)
  colnames(similar)<-colnames(data)
  for(i in 1:(nplot-1)){
    m=i+1
    for(m in m:nplot){
      if(index=="bc"){
        bc.sim=2*sum(apply(data[,c(i,m)], 1, min))/(sum (data[,c(i,m)]))
        similar[i,m]=bc.sim
        similar[m,i]=bc.sim
      }
    }
  }
  return(round(similar,3))
}

bray.mat <- sim(genus.rel)

save(bray.mat,file = "bray.its.out")

pca <- prcomp(t(bray.mat))
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]
plot(pc1,pc2)

pc.sum <- summary(pca)
var.pc1 <- pc.sum$importance[2,1]
var.pc2 <- pc.sum$importance[2,2]

pca.frame <- data.frame(PC1 = pc1, PC2 = pc2,Habitat = ifelse(grepl("r",colnames(bray.mat),fixed = T), "Endosphere", "Rhizosphere"))

pca.plot.fun <- ggplot(pca.frame, aes(x = PC1, y = PC2,fill = factor(Habitat))) +
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
  ylab(paste0("PC2 (", round(var.pc2 * 100, 1), "%)"))+ 
  annotate("text", x = -0.11, y = 0.13, label = "Fungi", 
           hjust = -0.1, vjust = 1.1, size = 6)

png(filename = "pca_fig_fun.png",width =500,height = 400)
pca.plot.fun
dev.off()

OTU.gwas<-read.table("otutable_fungigwasfinal.txt")
read.counts <- apply(OTU.gwas,2,sum)

OTU.gwas<-OTU.gwas[,-c(which(read.counts<1000))]
OTU.gwas <- OTU.gwas[,-which(grepl("clesn",colnames(OTU.gwas)))]
OTU.gwas <- OTU.gwas[,-which(grepl("215b",colnames(OTU.gwas)))]
OTU.gwas <- OTU.gwas[,-which(grepl("r240",colnames(OTU.gwas)))]
colnames(OTU.gwas)[which(colnames(OTU.gwas)%in%"s215a")] <- "s215"
OTU.gwas <- apply(OTU.gwas,2,function(x) x/sum(x))

### Calculate in how many samples each clade occurs.
occ <- apply(OTU.gwas,1,function(x) sum(x>0))

OTU.gwas <- OTU.gwas[occ>=1,]

### Calculate shannon without filtering
shannon <- apply(OTU.gwas,2,function(x) -sum(x*log(x+0.0001)))

shannon.fun <- data.frame(Shannon = shannon, Habitat = ifelse(grepl("r",colnames(bray.mat),fixed = T), "Endosphere", "Rhizosphere"))

fun.shannon.box <-  ggplot(shannon.fun, aes(x = factor(Habitat), y = Shannon)) +
  geom_violin(aes(fill = Habitat), alpha = 0.6, show.legend = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, fill = NA, color = "black", show.legend = FALSE) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  labs(fill = "Environment", y = "Fungal diversity")+
  theme_minimal() +
  theme(text = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        legend.position = "none")

library(cowplot)

legend <-  get_legend(ggplot(shannon.fun, aes(x = factor(Habitat), y = Shannon)) +
  geom_violin(aes(fill = Habitat), alpha = 0.6, show.legend = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, fill = NA, color = "black", show.legend = FALSE) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  labs(fill = "Environment", y = "Fungal diversity")+
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        legend.position = "top",
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22)))


library(tidyr)


### Make the barplots
load("genus_rel_filt_mitchlo.out")

genus.rel <- genus.rel[,-which(colnames(genus.rel)%in%c("GWr_170","GWs_179","GWr_065","GWs_018","GWs_074","GWs_001","GWs_002","GWs_120","GWr_014","GWr_189","GWs_181_clesn","GWs_220_clesn",
                                                        "GWs_64_clesn","GWs_215b","GWs_000"))]

tax <- genus.rel[,1:6]
genus.rel <- genus.rel[,-c(1:6)]
tax[tax == "unknown"] <- "Unknown"

Habitat = ifelse(grepl("GWr",colnames(genus.rel)),"Endosphere","Rhizosphere")

genus.rhizo <- genus.rel[,Habitat == "Rhizosphere"]
genus.endo <- genus.rel[,Habitat == "Endosphere"]

mean.rhizo <- apply(genus.rhizo,1,mean)
mean.endo <- apply(genus.endo,1,mean)

class.frame <- data.frame(class = tax$class, Rhizosphere = mean.rhizo, Endosphere = mean.endo)
class.frame <- aggregate(class.frame[,-which(colnames(class.frame)=="class")],by = list(class.frame$class),FUN = sum)
class.frame <- pivot_longer(class.frame,cols = c(2,3),names_to = "Habitat")
colnames(class.frame)[1]<- "class"
class.frame$class <- ifelse(class.frame$value <0.01,"Other",class.frame$class)
class.frame$class[class.frame$class=="Blastocatellia_(Subgroup_4)"] <- "Blastocellia_4"
class.frame$class<-factor(class.frame$class, levels=c(sort(unique(class.frame$class[!(class.frame$class%in%c("Other","Unknown"))])),"Other","Unknown"))

cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#000000",
          "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#54278F",
          "#0072B2","#D55E00","#CC79A7","#009E73","#F0E442","#56B4E9","#F4A582","#FF7F00",
          "#B2DF8A","#1F78B4","#FB9A99","#984EA3","darkgrey","lightgrey")

cols.short <- c(cols[c(1:18)],"grey30","lightgrey")
cols.short <- cols.short[-8]

class.plot.bac <- ggplot(class.frame, aes(x = Habitat, y = value, fill = class)) +
  geom_col(width = 0.8, colour = "grey30", size = 0.15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cols.short) +                    # or scale_fill_viridis_d(option="D")
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 20),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "right",
    legend.justification = c(0,1),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "white", colour = NA)
  ) +
  labs(y = "Proportion", fill = "Class")+
  guides(fill = guide_legend(ncol = 1))

png(filename = "barplot_class_bac.png",width = 500,height = 600)
class.plot.bac
dev.off()

###################### Barplot class rank fungi ########################

OTU.gwas<-read.table("otutable_fungigwasfinal.txt")
read.counts <- apply(OTU.gwas,2,sum)

load("agglomerate_counts_fungi/genus_rel_filt_mitchlo.out")

genus.rel<-genus.rel[,-c(which(read.counts<1000)+6)]
genus.rel <- genus.rel[,-which(grepl("clesn",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("215b",colnames(genus.rel)))]
genus.rel <- genus.rel[,-which(grepl("r240",colnames(genus.rel)))]
colnames(genus.rel)[which(colnames(genus.rel)%in%"s215a")] <- "s215"

tax <- genus.rel[,1:6]
tax[tax=="unidentified"] <- "Unknown"
genus.rel <- genus.rel[,-c(1:6)]

Habitat = ifelse(grepl("r",colnames(genus.rel)),"Endosphere","Rhizosphere")

genus.rhizo <- genus.rel[,Habitat == "Rhizosphere"]
genus.endo <- genus.rel[,Habitat == "Endosphere"]

mean.rhizo <- apply(genus.rhizo,1,mean)
mean.endo <- apply(genus.endo,1,mean)

library(tidyr)
library(ggplot2)

class.frame <- data.frame(class = tax$class, Rhizosphere = mean.rhizo, Endosphere = mean.endo)
class.frame <- aggregate(class.frame[,-which(colnames(class.frame)=="class")],by = list(class.frame$class),FUN = sum)
class.frame <- pivot_longer(class.frame,cols = c(2,3),names_to = "Habitat")
colnames(class.frame)[1]<- "class"
class.frame$class <- ifelse(class.frame$value <0.01,"Other",class.frame$class)
class.frame$class<-factor(class.frame$class, levels=c(sort(unique(class.frame$class[!(class.frame$class%in%c("Other","Unknown"))])),"Other","Unknown"))

cols.fun <- cols.short[c(1:13,length(cols.short)-1,length(cols.short))]

class.plot.fun <- ggplot(class.frame, aes(x = Habitat, y = value, fill = class)) +
  geom_col(width = 0.8, colour = "grey30", size = 0.15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cols.fun) +                    # or scale_fill_viridis_d(option="D")
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 20),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "right",
    legend.justification = c(0,1),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "white", colour = NA)
  ) +
  labs(y = "Proportion", fill = "Order")+
  guides(fill = guide_legend(ncol = 1))

png(filename = "barplot_class_fun.png",width = 800,height = 750)
class.plot.fun
dev.off()

fig.comb <- plot_grid(pca.plot.bac,bac.shannon.box,class.plot.bac,pca.plot.fun,fun.shannon.box,class.plot.fun,nrow = 2, rel_widths = c(0.6,0.4,0.65,0.6,0.65,0.4),labels = c("A","B","C","D","E","F"),
                    label_size = 25,label_y = 1.02,align = "h",axis = "b")

png(filename = "fig2.png",width =1100,height = 900)
fig.comb
dev.off()











