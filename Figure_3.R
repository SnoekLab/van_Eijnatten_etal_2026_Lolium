### This script outputs Figure 3. It includes some different versions of the figure, including the
### one in the final manuscript

library(ggplot2)
library(cowplot)

### Read in all GWAS data
endo.16s.gen.files <- list.files("GWAS_lolium_genus_16s_corgwas_corrected_v2/endo/frames",full.names = T)
endo.16s.pca.files <- list.files("GWAS_lolium_pcoa_16S_corrected_v3/root/frames",full.names = T)
endo.16s.div.files <- c(list.files("GWAS_lolium_shannon_16S_corrected_v2/root/frames",full.names = T),
                        list.files("GWAS_lolium_shannon_16S_corrected_v2/root_depth/frames",full.names = T))

rhizo.16s.gen.files <- list.files("GWAS_lolium_genus_16s_corgwas_corrected_v2/rhizo/frames",full.names = T)
rhizo.16s.pca.files <- list.files("GWAS_lolium_pcoa_16S_corrected_v3/soil/frames",full.names = T)
rhizo.16s.div.files <- c(list.files("GWAS_lolium_shannon_16S_corrected_v2/soil/frames",full.names = T),
                         list.files("GWAS_lolium_shannon_16S_corrected_v2/soil_depth/frames",full.names = T))

endo.its.gen.files <- list.files("GWAS_lolium_genus_its_corgwas_corrected_v2/endo/frames",full.names = T)
endo.its.pca.files <- list.files("GWAS_lolium_pcoa_ITS_corrected_v3/root/frames",full.names = T)
endo.its.div.files <- c(list.files("GWAS_lolium_shannon_fungi_corrected_v2/root/frames",full.names = T),
                        list.files("GWAS_lolium_shannon_fungi_corrected_v2/root_depth/frames",full.names = T))

rhizo.its.gen.files <- list.files("GWAS_lolium_genus_its_corgwas_corrected_v2/rhizo/frames",full.names = T)
rhizo.its.pca.files <- list.files("GWAS_lolium_pcoa_ITS_corrected_v3/soil/frames",full.names = T)
rhizo.its.div.files <- c(list.files("GWAS_lolium_shannon_fungi_corrected_v2/soil/frames",full.names = T),
                         list.files("GWAS_lolium_shannon_fungi_corrected_v2/soil_depth/frames",full.names = T))

lol.files <- c(endo.16s.gen.files,endo.16s.pca.files,endo.16s.div.files,
               rhizo.16s.gen.files,rhizo.16s.pca.files,rhizo.16s.div.files,
               endo.its.gen.files,endo.its.pca.files,endo.its.div.files,
               rhizo.its.gen.files,rhizo.its.pca.files,rhizo.its.div.files)

traits <- as.vector(sapply(lol.files,function(x)strsplit(x,"/",fixed = T)[[1]][10]))
traits <- as.vector(sapply(traits,function(x)strsplit(x,".",fixed = T)[[1]][1]))

bac <- ifelse(grepl("16S",lol.files,ignore.case = T),"Bacteria","Fungi")
comp <- ifelse(grepl("root",lol.files,fixed = T)|grepl("endo",lol.files,fixed = T),"Endosphere","Rhizosphere")

all.gwas <- data.frame()

for(i in 1:length(lol.files)){
  load(lol.files[i])
  plot_frame$Trait <- traits[i]
  plot_frame$Kingdom <- bac[i]
  plot_frame$Compartment <- comp[i]
  all.gwas <- rbind(all.gwas,plot_frame)
}

all.gwas$pval <- -log10(all.gwas$pval)
all.gwas <- all.gwas[all.gwas$pval>2,]
all.gwas$pos <- as.integer(as.character(all.gwas$pos*10^7))

### Add type so we can make the shape of the point show whether the QTL relates to a taxon,
### alpha diversity, or beta diversity
Type <- rep("Genus",nrow(all.gwas))
Type[grepl("Shannon",all.gwas$Trait)] <- "Shannon"
Type[grepl("PC",all.gwas$Trait)] <- "PC"

all.gwas$Type <- Type

all.man <- ggplot(all.gwas,aes(pos/10^6,pval,col = Compartment))+
  geom_point(size = 3)+
  geom_hline(yintercept = 7,col="red",linewidth=0.1)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350),expand = c(0.05,0.05))+
  scale_color_manual(values = c("#D55E00", "#0072B2"))+
  xlab("Position (Mbp)") + ylab(bquote(-log[10] (p)))+
  facet_grid(Kingdomchr,scales = "free_x",space = "free_x")+
  theme_cowplot()+
  theme(panel.background = element_rect(linewidth = 0.2,color="black"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.2),
        axis.text.x = element_text(angle = 90,size=25,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(lineheight =7, size = 25),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        panel.spacing.x = unit(3.5,"mm"),
        strip.text = element_text(face="bold", size=25),
        panel.grid.major = element_line(color = "grey2",linewidth =  0.1,linetype = 3),
        legend.position = "top",
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25),
        legend.justification = "center",
        panel.border = element_rect(linewidth = 0.2,color="black"))

png(filename = "man_all_v2.png",width = 1500,height = 1000)
all.man
dev.off()

### Load in the Bonferroni thresholds
load("bf.bac.out")
load("bf.fun.out")

all.max <- aggregate(all.gwas$pval,list(all.gwas$Trait,all.gwas$Compartment,all.gwas$Kingdom),max)
colnames(all.max) <- c("Trait","Compartment","Kingdom","P-value")
all.max$Type <- ifelse(grepl("PC",all.max$Trait),"Composition","Genus")
all.max$Type[grepl("Shannon",all.max$Trait)] <- "Diversity"
all.max <- all.max[-which(all.max$Trait=="Shannon"),]

all.gwas$sig_7 <- ifelse(all.gwas$pval<7,"not_sig",all.gwas$Kingdom)
all.gwas$sig_bf <- ifelse(all.gwas$pval< -log10(bf.fun[2]),"not_sig",all.gwas$Kingdom)

all.man.rev <- ggplot(all.gwas[all.gwas$pval>5,],aes(pos/10^6,pval,col = sig_bf))+
  geom_point(size = 3)+
  geom_hline(yintercept = -log10(bf.fun[2]),col="red",linewidth=0.3)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350),expand = c(0.05,0.05))+
  scale_color_manual(values = c("not_sig" = "lightgrey", "Bacteria" = "#D55E00", "Fungi" = "#0072B2"),
                     breaks = c("Bacteria", "Fungi"),
                     name = "Kingdom") +
  xlab("Position (Mbp)") + ylab(bquote(-log[10] (p)))+
  facet_grid(Compartmentchr,scales = "free_x",space = "free_x")+
  theme_cowplot()+
  theme(panel.background = element_rect(linewidth = 0.2,color="black"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.2),
        axis.text.x = element_text(angle = 90,size=25,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(lineheight =7, size = 25),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        panel.spacing.x = unit(3.5,"mm"),
        strip.text = element_text(face="bold", size=25),
        panel.grid.major = element_line(color = "grey2",linewidth =  0.1,linetype = 3),
        legend.position = "top",
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25),
        legend.justification = "center",
        panel.border = element_rect(linewidth = 0.2,color="black"))

png(filename = "man_all_rev_v2.png",width = 1500,height = 700)
all.man.rev
dev.off()

all.man.rev.type <- ggplot(all.gwas[all.gwas$pval>5,],aes(pos/10^6,pval,col = sig_bf,shape = Type))+
  geom_point(size = 5,alpha = 0.9)+
  geom_hline(yintercept = -log10(bf.fun[2]),col="red",linewidth=0.3)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400),expand = c(0.05,0.05))+
  scale_color_manual(values = c("not_sig" = "lightgrey", "Bacteria" = "#D55E00", "Fungi" = "#0072B2"),
                     breaks = c("Bacteria", "Fungi"),
                     name = "Kingdom") +
  xlab("Position (Mbp)") + ylab(bquote(-log[10] (p)))+
  facet_grid(Compartmentchr,scales = "free_x",space = "free_x")+
  theme_cowplot()+
  theme(panel.background = element_rect(linewidth = 0.2,color="black"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.2),
        axis.text.x = element_text(angle = 90,size=25,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(lineheight =7, size = 25),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        panel.spacing.x = unit(3.5,"mm"),
        strip.text = element_text(face="bold", size=25),
        panel.grid.major = element_line(color = "grey2",linewidth =  0.1,linetype = 3),
        legend.position = "top",
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25),
        legend.justification = "center",
        panel.border = element_rect(linewidth = 0.2,color="black"))

png(filename = "man_all_rev_type_v2.png",width = 1500,height = 700)
all.man.rev.type
dev.off()

save(all.gwas,file = "all.gwas.out")







