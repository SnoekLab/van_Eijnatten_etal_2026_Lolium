### This script makes figure 3.

library(ggplot2)
library(pals)
library(scales)

### Read in the 16S QTLs
load("16S_rhizo_submission.out")

### Read in the GWAS results of the traits mapping to the pleiotropic locus
gwas.files <- c(list.files("GWAS_full_pleitropic_locus_corrected_v2/frames",full.names = T))
gwas.files <- gwas.files[!grepl("PC2",gwas.files)&!grepl("PC3",gwas.files)]

random.frame <- data.frame()

for(i in gwas.files){
  load(i)
  plot_frame <- plot_frame[!is.na(plot_frame$pval),]
  plot_frame$pval <- -log10(plot_frame$pval)
  if(nrow(plot_frame)!=0){
    trait <- strsplit(i,"frames/",fixed = T)[[1]][2]
    trait <- strsplit(trait,".out",fixed = T)[[1]][1]
    plot_frame$trait <- trait
    random.frame <- rbind(random.frame,plot_frame)
  }
}

random.frame$chr_pos <- paste(random.frame$chr,random.frame$pos)

### Read in the GFF
gff <- read.table("Lolium_2.6.1_V3_gene_annotation.gff3")
colnames(gff)[c(1,3,4,5)] <- c("Chromosome","Type","Start","Stop")

### Select genes in a 0.1 Mbp interval around the peak SNP
genes <- gff[gff$Chromosome=="Lp_chr7_0.1"&gff$Start>29.03*10^7&gff$Stop<29.053*10^7&gff$Type=="gene","V9"]

### Read in candidate gene data
hom.frame <- read.csv("Arabidopsis_homologs_near_QTLs_16S.csv")
hom.frame[hom.frame$Gene_start%in%gff[gff$V9%in%genes,"Start"],]

### Plot the locus (0.2Mbp interval) of the pleiotropic locus with a separate manhattan plot 
### for each trait (figure not in paper).
cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#E41A1C","#0072B2")

loc.plot.full <- ggplot(random.frame[random.frame$pos>29.04440&random.frame$pos<29.0475&random.frame$chr == 7 ,], aes(x = pos*10, y = pval, fill = trait)) +
  geom_point(shape = 24, colour = "black", size = 4.5, alpha = 0.6) +
  scale_fill_manual(values = cols) +
  facet_wrap(~trait)+
  scale_x_continuous(expand = expansion(add = 0.00002),breaks = c(290.445, 290.455,290.465,290.475)) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", colour = "black", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    legend.key = element_rect(fill = "white", colour = NA)
  ) +
  labs(x = "Position (Mbp)", y = "P-value", fill = "") 


png(filename = "Pleiotropic_locus_narrow.png",width = 800,height = 500)
loc.plot.full
dev.off()

library(scales)
library(cowplot)

### Select the locus
sub <- random.frame[random.frame$pos>(29.0459937-0.03)&random.frame$pos<29.0459937+0.03&random.frame$chr == 7&random.frame$pval>2 ,]
sub <- sub[order(sub$pval),]

### Select a color ramp for the p-values
pal <- c("grey90","grey70","cyan2","green3","orange","hotpink")
leg_breaks <- c(3,4,5,6,7,8)
col_fun <- col_numeric(pal, domain = c(2,8.34))

### Select the candidate genes in the interval
genes_0.1_MBp <- gff[gff$Start>(290459937-0.1*10^6)&gff$Stop<(290459937+0.1*10^6)&gff$Chromosome=="Lp_chr7_0.1"&gff$Type=="gene",]
hom.frame[hom.frame$Gene_start%in%genes_0.3_MBp$Start,] 
genes_0.1_MBp$Gene_name <- c("TBCC1","FTIP3","FTIP1","NAMT1")

### Define interval around peak SNP.
center <- 290.459937
lims <- c(center - 0.1, center + 0.1)

### Plot the locus
sub$trait[sub$trait == "Shannon depth"] <- "Shannon diversity"
sub$trait <- factor(sub$trait,levels = rev(c("PC1","Shannon diversity","Tumebacillus","Acidothermus","Gaiella","mle1-7","RB41","Pseudonocardia","Pseudolabrys")))                    
                    
loc.plot.full <- ggplot(sub, aes(x = pos*10, y = trait, fill = pval,size = pval)) +
  geom_point(shape = 24, alpha = 0.6,color = "black") +
  scale_fill_gradientn(colors = pal,guide = "none")+
  scale_size_continuous(
    range = c(2, 8),
    breaks = leg_breaks,
    guide = guide_legend(
      title = "Pval",
      override.aes = list(
        shape = 24,
        colour = "black",
        fill = col_fun(leg_breaks),   # set fill colors for legend keys
        alpha = 0.9
      )
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +  
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", colour = "black", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    legend.key = element_rect(fill = "white", colour = NA)
  ) +
  labs(x = "Position on chr 7 (Mbp)", y = "Phenotype", fill = "Pval") +
  coord_cartesian(xlim = lims)


png(filename = "~/Documents/Lolium/LoliumSNPs/Submission_figures/Pleiotropic_locus_narrow_different.png",width = 800,height = 350)
loc.plot.full
dev.off()

### Plot the gene track below the figure
genes <- genes_0.1_MBp
genes$startM <- genes$Start / 1e7 * 10   
genes$endM   <- genes$Stop  / 1e7 * 10
genes$midM   <- (genes$startM + genes$endM)/2

gene.track <- ggplot(genes) +
  geom_rect(aes(xmin = startM, xmax = endM, ymin = 0.8, ymax = 1.2),
            fill = "grey70", colour = "black") +
  geom_text(aes(x = startM, y = 0.6, label = Gene_name),
            size = 4, angle = 0, hjust = 0.2) +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 5, b = 10, l = 5)) +
  coord_cartesian(xlim = lims,ylim = c(0.5, 1.2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))

### Paste the gene track below the locus plot
cow.plt <- plot_grid(loc.plot.full, gene.track, ncol = 1, rel_heights = c(0.85, 0.15),
          align = "v")+panel_border(color = "black",size = 0.2)

### Save figure
png(filename = "Pleiotropic_locus_with_track.png",width = 800,height = 350)
cow.plt
dev.off()



