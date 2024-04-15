#!/bin/bash


#referneces
genome_reference_directory="Data/Dev/pdx5/Data/Referance_data/GRCh38_Ref_Fasta";
reference="GRCh38_no_alt_analysis_set.fna";
variants_reference_directory="Data/Dev/pdx5/Data/Referance_data/GRCh38_commom_variants";
known_variants="common_known_sites.vcf.gz";

#tools
code_root="/home/pdx5/Clinical_Genomics/Code/";
bwa_root="/home/packages/bwa/";
gatk_root="/home/packages/gatk-4.5.0.0/";
picard_root="/home/packages/picard.jar";
strelka_root="/home/packages/strelka-2.9.10.centos6_x86_64/bin/";
func_germline="/home/pdx5/metagenomics_packages/funcotator_dataSources.v1.8.hg38.20230908g";
func_somatic="/home/packages/gatk-4.5.0.0/funcotator_dataSources.v1.8.hg38.20230908s";
snpeff_root="/home/pdx5/metagenomics_packages/snpEff/";
sift_root="/home/pdx5/metagenomics_packages/snpEff/SnpSift.jar";
vep_root="/home/pdx5/metagenomics_packages/ensembl-vep/";
trimmomatic_root="/home/packages/Trimmomatic-0.39/trimmomatic-0.39.jar";
trimmo_adapters="/home/Data/Dev/pdx5/Data/Referance_data/trimo_all_adapters.fa";
qualimap_root="/home/packages/qualimap_v2.3/";
python2_env="py2";
java_root="/home/pdx5/metagenomics_packages/jdk-21.0.2/bin/";
bedtools_root="/home/packages/bedtools2/bin/";

#Output results root
output_root="/home/Data/Dev/pdx5/Results";
final_stats="$output_root/Final_statistics";


#GATK filters
QD="2.00";
DP="10.00";
FS="60.00";
SOR="4.0";
GQ="10";

#Strelka filters
SB="10";

#parameters
#bwa-mem2 therads
threads="32";
#picard cpu usage
cpu="16";




