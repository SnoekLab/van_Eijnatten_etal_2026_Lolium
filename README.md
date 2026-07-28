# van_Eijnatten_etal_2026_Lolium

This repository contains the code used in the analysis of van Eijnatten et al., 2026 (add doi). For any questions regarding the code or data availability, please reach out to a.l.vaneijnatten@uu.nl.

PCoA_for_GWAS_and_FigureS2.R does a PCoA on each community (16S-endosphere, 16S-rhizosphere,ITS-endosphere, and ITS-rhizosphere) and saves the principal components for use in GWAS. This script also generates Figure S2.

The folder GWAS_scripts contains scripts performing the GWAS. There are separate scripts for the 16S and ITS traits and for the relative abundance, Shannon diversity, or PCoA traits.

The folder QTLs_and_candidate genes contains the scripts to score QTLs and discover condidate genes. Get_peak_SNPs_16S.R and Get_peak_SNPs_ITS.R select all peak SNPs for each significant association. Then, Make_fastasa_and_blastx.R takes the sequences of all genes in a 0.2Mbp interval around each peak SNP and puts them in a fasta for use in BLAST. The annotation of this script also contains the command line code we used to download the ITAG4.0 Arabidopsis proteome and to find A.thaliana homologs using blastx. The script Collect_AT_IDs_from_NCBI.R then uses the NCBI accession numbers output by blastx to obtain the AT identifiers and gene IDs of the the Arabidopsis homologs. Finally, Make_files_with_candidate_genes.R processes this data and outputs an excel with one row per gene, showing the 10 closest A.thaliana homologs (based on bitscore). 

The Mantel_tests folder contains the script made to make genetic correlation matrices and the script used to do mantel tests between the genetic distance (1 - correlation matrix) and the microbiome distance (bray-curtis). 

The script Calculate_Benjamini_Hochberg_pvalues.R takes the GWAS results of the traits mapping to the pleiotropic locus as input and calculates the Benjamini-Hochberg corrected p-values.

The remaining scripts generate the figures in the manuscript.

