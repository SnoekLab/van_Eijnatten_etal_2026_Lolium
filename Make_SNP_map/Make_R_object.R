### This script reads in the filtered vcf files (obtained with the script Process_vcf.py) and saves these as Robjects. This includes separate
### R objects per chromosome as well as an object covering the whole genome.

library(viridis)
library(openxlsx)

ubec1 <- read.delim("Lolium_chr1.vcf", header = F)
ubec1[ubec1=="11"] <- 3 ; ubec1[ubec1=="00"] <- 1 ; ubec1[ubec1=="01"] <- 2 ; ubec1[ubec1==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec1) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec1,file="obj_ubec1.out")
load(file="obj_ubec1.out")
object.size(ubec1) ## 3.5 GB !!!!!


rm(ubec1)
ubec2 <- read.delim("Lolium_chr2.vcf", header = F)
ubec2[ubec2=="11"] <- 3 ; ubec2[ubec2=="00"] <- 1 ; ubec2[ubec2=="01"] <- 2 ; ubec2[ubec2==".."] <- 9


colnames(ubec2) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)
save(ubec2,file="obj_ubec2.out")

rm(ubec2)
ubec3 <- read.delim("Lolium_chr3.vcf", header = F)
ubec3[ubec3=="11"] <- 3 ; ubec3[ubec3=="00"] <- 1 ; ubec3[ubec3=="01"] <- 2 ; ubec3[ubec3==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec3) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec3,file="obj_ubec3.out")
object.size(ubec3) ## 4.3 GB !!!!!

unique(ubec3[,18])

rm(ubec3)
ubec4 <- read.delim("Lolium_chr4.vcf", header = F)
ubec4[ubec4=="11"] <- 3 ; ubec4[ubec4=="00"] <- 1 ; ubec4[ubec4=="01"] <- 2 ; ubec4[ubec4==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec4) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec4,file="obj_ubec4.out")
object.size(ubec4) ## 5.4 GB !!!!!

rm(ubec4)
ubec5 <- read.delim("Lolium_chr5.vcf", header = F)
ubec5[ubec5=="11"] <- 3 ; ubec5[ubec5=="00"] <- 1 ; ubec5[ubec5=="01"] <- 2 ; ubec5[ubec5==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec5) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec5,file="obj_ubec5.out")

rm(ubec5)
ubec6 <- read.delim("Lolium_chr6.vcf", header = F)
ubec6[ubec6=="11"] <- 3 ; ubec6[ubec6=="00"] <- 1 ; ubec6[ubec6=="01"] <- 2 ; ubec6[ubec6==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec6) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec6,file="obj_ubec6.out")
object.size(ubec6) ## 3.5 GB !!!!!

rm(ubec6)
ubec7 <- read.delim("Lolium_chr7.vcf", header = F)
ubec7[ubec7=="11"] <- 3 ; ubec7[ubec7=="00"] <- 1 ; ubec7[ubec7=="01"] <- 2 ; ubec7[ubec7==".."] <- 9

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
colnames(ubec7) <- c("Chr","Pos","REF","ALT","QUAL","FILTER","AC","AF","AN","DP",1:239)

save(ubec7,file="obj_ubec7.out")
rm(ubec7)

load("obj_ubec1")
load("obj_ubec2")
load("obj_ubec3")
load("obj_ubec4")
load("obj_ubec5")
load("obj_ubec6")
load("obj_ubec7")

ubec <- rbind(ubec1,ubec2,ubec3,ubec4,ubec5,ubec6,ubec7)

save(ubec,file = "ubec.out")




