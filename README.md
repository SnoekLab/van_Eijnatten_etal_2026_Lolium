# van_Eijnatten_etal_2026_Lolium

This repository contains the code used in the analysis of van Eijnatten et al., 2026 (add doi). For any questions regarding the code or data availability, please reach out to a.l.vaneijnatten@uu.nl.

The SNP map used in this study was based on the VCF file from Fois et al., 2021 (https://doi.org/10.3390/genes13010020). The folder Make_SNP_map contains the scripts we used to process the VCF file. Process_vcf.py loads in the vcf file and filters the SNPs based on quality and allele frequency. Make_R_object.R loads the filtered VCF in R and saves it as an R object.

PCoA_for_GWAS_and_FigureS2.R performs PCoA analysis on each community (16S-endosphere, 16S-rhizosphere,ITS-endosphere, and ITS-rhizosphere) and saves the principal components for use in GWAS. This script also generates Figure S2.

The folder GWAS_scripts contains scripts performing the GWAS. There are separate scripts for the 16S and ITS traits and for the relative abundance, Shannon diversity, and PCoA traits.

The folder QTLs_and_candidate genes contains the scripts to score QTLs and discover condidate genes. The genome assembly and gff used was made available by Bryce et al., 2015 (https://doi.org/10.1111/tpj.13037) and can be downloaded from https://ryegrassgenome.ghpc.au.dk/DOWNLOAD/. Get_peak_SNPs_16S.R and Get_peak_SNPs_ITS.R select all peak SNPs for each significant association. Then, Make_fastasa_and_blastx.R takes the sequences of all genes in a 0.2Mbp interval around each peak SNP and puts them in a fasta for use in BLAST. The command line code we used to download the ITAG4.0 Arabidopsis proteome and to find A.thaliana homologs using blastx is as follows:

conda create -n blast-env blast python=3\n
source activate blast-env\n
conda install -c conda-forge ncbi-datasets-cli\n
datasets download genome taxon "Arabidopsis thaliana" --include protein --filename arabidopsis_proteome.zip\n
unzip arabidopsis_proteome.zip -d arabidopsis_data\n
makeblastdb -in arabidopsis_data/ncbi_dataset/data/GCA_000001735.2/protein.faa -dbtype prot -out ath_protein_db\n
blastx -query ITS_seqs.fasta -db ath_protein_db -out ITS_blast_results.fmt6 -evalue 1e-5 -num_threads 15 -outfmt 6\n

Next, the script Collect_AT_IDs_from_NCBI.R then uses the NCBI accession numbers output by blastx to obtain the AT identifiers and gene IDs of the the Arabidopsis homologs. Finally, Make_files_with_candidate_genes.R processes this data and outputs an excel with one row per gene, showing the 10 closest A.thaliana homologs (based on bitscore). 

The Mantel_tests folder contains the script used to make genetic correlation matrices and the script used to do mantel tests between the genetic distance (1 - correlation matrix) and the microbiome distance (bray-curtis). 

The script Calculate_Benjamini_Hochberg_pvalues.R takes the GWAS results of the traits mapping to the pleiotropic locus as input and calculates the Benjamini-Hochberg corrected p-values.

The remaining scripts generate the figures in the manuscript.
