### This script makes figure S1.

library(ggplot2)
library(openxlsx)
library(reshape2)
library(scico)
library(cowplot)

################## Load in correlation matrix.
load("cor.mat.out")

hm.gplots <- heatmap.2(cor.ubec,scale = "none",trace = "none",col = rainbow(256),Rowv = T,Colv = T,
                       density.info = "none",dendrogram = "both",margins = c(0,0))

hr <- as.hclust(hm.gplots$rowDendrogram)
hr$labels

diag(cor.ubec) <- NA

df <- melt(cor.ubec)   

hm.plot.t <- ggplot(df,aes(x= factor(Var1,levels = hr$labels[hr$order]),y=factor(Var2,levels = hr$labels[hr$order]),fill=value)) +
  geom_tile() +
  labs(x="",y="",fill = "Correlation")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_scico(palette = 'vik') +
  xlab("Environmental group")+
  theme(legend.position="top", 
        legend.key.width=unit(.1,"npc"),legend.key.height=unit(.05,"npc"),axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

png(filename = "hm_cor.png",width = 600, height = 600)  
plot_grid(hm.plot.t) + panel_border(color = "black",size = 1)
dev.off()

range(cor.ubec,na.rm =T)
quantile(cor.ubec[upper.tri(cor.ubec)],probs = 0.95,na.rm = T)
