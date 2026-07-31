### This script reads in the  vcf from Fois et al., 2021 (https://doi.org/10.3390/genes13010020) and filters the SNP on quality and allele
### frequency. The filtered SNPs are then saved to 7 VCFs, one for each chromosome.

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#### Test and inspection run
"""
lol = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r") 
lol_smallsample = open("X:/Work/Lolium_SNPS/Lolium_Small.txt","w")
for i in range(10000):
    ok = lol.readline()
    lol_smallsample.writelines(ok)
    
lol.close()
lol_smallsample.close()


"""
"""
QUAL_thr = 200
AC_thr = 2
AF_thr = 0

tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")
redtki = open("X:/Work/Lolium_SNPS/Lolium_red2.vcf","w")
### loop to skip over the 9187 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()

## we need chr, pos, 
for j in range(20000):
    myline = tki.readline()
    # print(myline)
    myline = myline.rstrip()
    myline = myline.split("\t")
    # print(myline)
    newline = []
    newline += [myline[0][-3]]  # chromosome
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]
           
    # myline[8]  ## order info can be removed
    for k in range(9,len(myline)):
        newline += [myline[k][0:3:2]]
    newline = "\t".join(newline)
    newline += "\n"
    # print(newline)
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr:
        redtki.write(newline)
    
tki.close()
redtki.close()   
"""

# readlines()

##### Make while loop to process chr 1
#"""
QUAL_thr = 200
AC_thr = 2
AF_thr = 0.00
AF_upper = 1
use_chr = "1"

tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")
redtki = open("X:/Work/Lolium_SNPS/Lolium_chr1.vcf","w")

### loop to skip over the 48 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()

## we need chr, pos, 
while use_chr == "1" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    #print(myline)
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    #print (use_chr)
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]
           
    # myline[8]  ## order info can be removed
    
    # print(newline)
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
    
tki.close()
redtki.close()   

print(use_chr)



##### Make while loop to process chr 2

QUAL_thr = 200
AC_thr = 2
AF_thr = 0.00
AF_upper = 1
use_chr = "1"

tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")
redtki = open("X:/Work/Lolium_SNPS/Lolium_chr2.vcf","w")

### loop to skip over the 9187 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()
while use_chr == "1" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    use_chr = myline[0][-3]

## we need chr, pos, 
while use_chr == "2" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]
           
    # myline[8]  ## order info can be removed
    
    # print(newline)
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
    
tki.close()
redtki.close()   

print(use_chr)






##### Make while loop to process chr 3 ------------------------------

QUAL_thr = 200
AC_thr = 2
AF_thr = 0.00
AF_upper = 1
use_chr = "1"

tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")
redtki = open("X:/Work/Lolium_SNPS/Lolium_chr3.vcf","w")

### loop to skip over the 9187 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()
while use_chr == "1" or use_chr == "2" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    use_chr = myline[0][-3]

## we need chr, pos, 
while use_chr == "3" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]
           
    # myline[8]  ## order info can be removed
    
    # print(newline)
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
    
tki.close()
redtki.close()   

print(use_chr)
print(colheader)





###################################################################################################################################
##
##
###################################################################################################################################
##### Make while loop to process rest of the chomosomes ------------------------------

QUAL_thr = 200
AC_thr = 2
AF_thr = 0.00
AF_upper = 1
use_chr = "1"


tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")

### loop to skip over the 9187 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()
while use_chr == "1" or use_chr == "2" or use_chr == "3":
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    use_chr = myline[0][-3]

## we need chr, pos, 
redtki = open("X:/Work/Lolium_SNPS/Lolium_chr4.vcf","w")

# redtki = open("X:/Work/TKI100/redtki.chr4.vcf","w")
print(use_chr)

while use_chr == "4" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

redtki = open("X:/Work/Lolium_SNPS/Lolium_chr5.vcf","w")

# redtki = open("X:/Work/TKI100/redtki.chr5.vcf","w")
print(use_chr)

while use_chr == "5" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

redtki = open("X:/Work/Lolium_SNPS/Lolium_chr6.vcf","w")

# redtki = open("X:/Work/TKI100/redtki.chr6.vcf","w")
print(use_chr)

while use_chr == "6" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

redtki = open("X:/Work/Lolium_SNPS/Lolium_chr7.vcf","w")

# redtki = open("X:/Work/TKI100/redtki.chr7.vcf","w")
print(use_chr)

while use_chr == "7" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-3]]  # chromosome
    use_chr = myline[0][-3]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

### last chr == 7 redo###########################################################################

QUAL_thr = 200
AC_thr = 2
AF_thr = 0.00
AF_upper = 1
use_chr = "1"


tki = open("X:/Work/Lolium_SNPS/Lolium_resequencing_v1.vcf","r")

### loop to skip over the 9187 lines of annotation
for i in range(9187):
    tki.readline()

colheader = tki.readline()
while use_chr == "1" or use_chr == "2" or use_chr == "3" or use_chr == "4" or use_chr == "5" or use_chr == "6":
    myline = tki.readline()
    #myline = myline.rstrip()
    #myline = myline.split("\t")
    use_chr = myline[6]

## we need chr, pos, 
redtki = open("X:/Work/Lolium_SNPS/Lolium_chr7.vcf","w")

# redtki = open("X:/Work/TKI100/redtki.chr4.vcf","w")
print(use_chr)
while use_chr == "7" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-5]]  # chromosome
    use_chr = myline[0][-5]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

"""
redtki = open("X:/Work/01_test-variant-calling_V1_20220414/redwei.chr8.vcf","w")
# redtki = open("X:/Work/TKI100/redtki.chr8.vcf","w")
print(use_chr)

while use_chr == "8" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-1]]  # chromosome
    use_chr = myline[0][-1]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

redtki = open("X:/Work/01_test-variant-calling_V1_20220414/redwei.chr9.vcf","w")
# redtki = open("X:/Work/TKI100/redtki.chr9.vcf","w")
print(use_chr)

while use_chr == "9" :
    myline = tki.readline()
    myline = myline.rstrip()
    myline = myline.split("\t")
    newline = []
    newline += [myline[0][-1]]  # chromosome
    use_chr = myline[0][-1]
    newline += [myline[1]] # position
    # myline[2]  ## can be removed
    newline += [myline[3]] ## REF allele
    newline += [myline[4]] ## ALT allele
    newline += [myline[5]] ## QUAL, we could filter on this 
    QUAL = float(myline[5]) ## for use as selction threshold
    newline += [myline[6]] ## FILTER, these names can be shortened
    tmpline7 = myline[7].split(";") ## we can loop over these to extract AC, AF, AN and DP
    for l in tmpline7:
       if l[0:2] == "AC":
           newline += [l[3:len(l)]]
           AC = int(l[3:len(l)].split(",")[0])
       if l[0:2] == "AF":
           newline += [l[3:len(l)]]
           AF = float(l[3:len(l)].split(",")[0])
       if l[0:2] == "AN":
           newline += [l[3:len(l)]]
       if l[0:2] == "DP":
           newline += [l[3:len(l)]]           
    # myline[8]  ## order info can be removed
    if QUAL > QUAL_thr and AC > AC_thr and AF > AF_thr and AF < AF_upper:
       for k in range(9,len(myline)):
           newline += [myline[k][0:3:2]]
       newline = "\t".join(newline)
       newline += "\n"
        
       redtki.write(newline)
redtki.close()   

tki.close()

"""


##### Count snps per chr
tki = open("X:/Work/TKI100/wageningen-reseq-iap.filtered_variants.vcf","r")
### loop to skip over the 123 lines of annotation
for i in range(123):
    tki.readline()
colheader = tki.readline()
snp_cnt = {}
for myline in tki:
    myline = myline.rstrip()
    myline = myline.split("\t")
    use_chr = myline[0][-1]
    if use_chr not in snp_cnt:
        snp_cnt[use_chr] = 1
    else:
        snp_cnt[use_chr] +=1
tki.close
        
print(snp_cnt)
# {'1': 25520618, '2': 26515477, '3': 28987834, '4': 44038378, '5': 39826839, '6': 22799217, '7': 23067848, '8': 35873282, '9': 24305699}
"""
############### END #####################################
