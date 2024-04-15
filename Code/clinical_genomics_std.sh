#!/bin/bash
source /home/pdx5/Clinical_Genomics/Code/config.sh

# Record the start time
start_time=$(date +%s)

echo -e "\n\n\n";
echo -e "\n";
echo "#######################################";
echo "#variables and intialization#";
echo "#######################################";
echo -e "\n";


#Root Directories
: << COMMENT
#pipeline_root="";
gatk_root="/home/saurabh/Clinical_Genomics/Tools/gatk-4.5.0.0/";
picard_root="/home/saurabh/Clinical_Genomics/Tools/";
strelka_root="/home/saurabh/anaconda3/envs/py2/bin/";
func_germline="/home/saurabh/Clinical_Genomics/Tools/gatk-4.5.0.0/funcotator_dataSources.v1.8.hg38.20230908g";
func_somatic="/home/saurabh/Clinical_Genomics/Tools/gatk-4.5.0.0/funcotator_dataSources.v1.8.hg38.20230908s";
snpeff_root="/home/saurabh/Clinical_Genomics/Tools/snpEff/";
vep_root="/home/saurabh/Clinical_Genomics/Tools/ensembl-vep/";
merge_root="/home/saurabh/Clinical_Genomics/Code/";
trimmomatic_root="/home/saurabh/Clinical_Genomics/Tools/Trimmomatic-0.39/";

COMMENT

# Check if the required number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample> <read1> <read2>"
    exit 1
fi

# Assign command-line arguments to variables

sample="$1"
sample=$(basename "$sample")

read1="$2"
seq1=$(basename "$read1" .fastq)

read2="$3"
seq2=$(basename "$read2" .fastq)


# Creating a directory based on sequence names
output_var="${output_root}/${sample}_variants_output"
output_seq="${output_root}/${sample}_sequences"
gatk_output="${output_var}/gatk_output"
strelka_output="${output_var}/strelka_output"
output_anno="${output_root}/${sample}_annotation_output"
output_final="${output_root}/${sample}_final_output"

mkdir -p "$output_seq"
mkdir -p "$output_var"
mkdir -p "$gatk_output"
mkdir -p "$output_anno"
mkdir -p "$output_final"
mkdir -p $final_stats

: <<'COM'
echo -e "Commented Block"
COM

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check for $seq1 of $sample #";
echo "#######################################";
echo -e "\n";
fastqc --outdir $output_seq $read1;
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check for $seq2 of $sample #";
echo "#######################################";
echo -e "\n";
fastqc --outdir $output_seq $read2;
date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Trimming reads with Trimmomatic#";
echo "#######################################";
echo -e "\n";

java -jar "$trimmomatic_root" PE -threads $threads "$read1" "$read2" \
"$output_seq"/"${seq1}"_paired.fastq.gz "$output_seq"/"${seq1}"_unpaired.fastq.gz \
"$output_seq"/"${seq2}"_paired.fastq.gz "$output_seq"/"${seq2}"_unpaired.fastq.gz \
ILLUMINACLIP:"$trimmo_adapters":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check of "${seq1}"_paired.fastq.gz of $sample #";
echo "#######################################";
echo -e "\n";
fastqc "$output_seq"/"${seq1}"_paired.fastq.gz;
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check of "${seq2}"_paired.fastq.gz of $sample #";
echo "#######################################";
echo -e "\n";
fastqc "$output_seq"/"${seq2}"_paired.fastq.gz;
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Unzipping all the fastqc reports for "${sample}"#";
echo "#######################################";
echo -e "\n";

for zipfile in "$output_seq"/*.zip; do
    echo "unzipping"
    unzip -qq -d "$output_seq" "$zipfile"
done

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Indexig and creating dictionary for reference genome $reference#";
echo "#######################################";
echo -e "\n";

# Construct the full paths for indexed files
reference_index="${genome_reference_directory}/${reference}.fai"
reference_amb="${genome_reference_directory}/${reference}.amb"
reference_ann="${genome_reference_directory}/${reference}.ann"
reference_bwt="${genome_reference_directory}/${reference}.bwt"
reference_pac="${genome_reference_directory}/${reference}.pac"
reference_sa="${genome_reference_directory}/${reference}.sa"
known_variants_index="${variants_reference_directory}/${known_variants}.tbi"

# Check if reference genome is indexed, if not, perform indexing
if [ ! -e "$reference_index" ] || [ ! -e "$reference_amb" ] || [ ! -e "$reference_ann" ] || [ ! -e "$reference_bwt" ] || [ ! -e "$reference_pac" ] || [ ! -e "$reference_sa" ]; then
    samtools faidx "${genome_reference_directory}/${reference}"
    "${bwa_root}"bwa index "${genome_reference_directory}/${reference}"
    java -jar "${picard_root}" CreateSequenceDictionary R="${genome_reference_directory}/${reference}"
    echo "Reference genome-$reference indexed successfully."
else
    echo "Reference genome-$reference already indexed. Skipping..."
fi

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Indexing known variants file $known_variants #";
echo "#######################################";
echo -e "\n";

# Check if known variants file is indexed, if not, perform indexing
if [ ! -e "$known_variants_index" ]; then
    tabix -fp vcf "${variants_reference_directory}/${known_variants}"
    echo "Known variants file-$known_variants indexed successfully."
else
    echo "Known variants file-$known_variants already indexed. Skipping..."
fi


date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#mapping/alignment of reads with reference genome#";
echo "#######################################";
echo -e "\n";
echo "${genome_reference_directory}/${reference}"
"${bwa_root}"bwa mem -M -t $threads "${genome_reference_directory}/${reference}" "$output_seq"/"${seq1}"_paired.fastq.gz "$output_seq"/"${seq2}"_paired.fastq.gz > "$output_var"/"${sample}"_mapped.sam

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#convert sam to bam#";
echo "#######################################";
echo -e "\n";

samtools view -b "$output_var"/"${sample}"_mapped.sam -o "$output_var"/"${sample}"_mapped.bam

date +"%D %T":"##Process Completed Time##";
 
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#sorting bam file#";
echo "#######################################";
echo -e "\n";

samtools sort "$output_var"/"${sample}"_mapped.bam -o "$output_var"/"${sample}"_mapped_sorted.bam;


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#remove duplicates#";
echo "#######################################";
echo -e "\n";

java -jar "${picard_root}" MarkDuplicates \
    INPUT="$output_var/${sample}_mapped_sorted.bam" \
    OUTPUT="$output_var/${sample}_mapped_sorted_duprem.bam" \
    METRICS_FILE="$output_var/${sample}_duplicate_metrics.txt" \
    REMOVE_DUPLICATES=false 
	

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#add read groups to bam file#";
echo "#######################################";
echo -e "\n";

java -jar "${picard_root}" AddOrReplaceReadGroups \
		I="$output_var/${sample}_mapped_sorted_duprem.bam" \
		O="$output_var/${sample}_mapped_sorted_duprem_RG.bam" \
		RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM="${sample}";


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing final bam file#";
echo "#######################################";
echo -e "\n";

samtools index "$output_var/${sample}_mapped_sorted_duprem_RG.bam";

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Displaying stats of BAM file#";
echo "#######################################";
echo -e "\n";

samtools stats "$output_var/${sample}_mapped_sorted_duprem_RG.bam" > "$output_var"/"${sample}"_stats.txt;

date +"%D %T":"##Process Completed Time##";



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 1 (BQSR)#";
echo "#######################################";
echo -e "\n";

#Creating base recal table

"${gatk_root}"./gatk BaseRecalibrator \
		-I "$output_var/${sample}_mapped_sorted_duprem_RG.bam" \
		-R "${genome_reference_directory}/${reference}" \
		--known-sites "${variants_reference_directory}/${known_variants}" \
		-O "$output_var"/"${sample}"_base_recal_1.table;

date +"%D %T":"##Process Completed Time##";


#Apply BQSR
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#apply BQSR#";
echo "#######################################";
echo -e "\n";

"${gatk_root}"./gatk ApplyBQSR \
		-R "${genome_reference_directory}/${reference}" \
		-I "$output_var/${sample}_mapped_sorted_duprem_RG.bam" \
		--bqsr-recal-file "$output_var"/"${sample}"_base_recal_1.table \
		-O "$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 2 (BQSR)#";
echo "#######################################";
echo -e "\n";

#Creating base recal table 2 for covarities

"${gatk_root}"./gatk BaseRecalibrator \
		-I "$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam" \
		-R "${genome_reference_directory}/${reference}" \
		--known-sites "${variants_reference_directory}/${known_variants}" \
		-O "$output_var"/"${sample}"_base_recal_2.table;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Analyzing Covariaties#";
echo "#######################################";
echo -e "\n";

#Analyzing Covarities

"${gatk_root}"./gatk AnalyzeCovariates \
		-before "$output_var"/"${sample}"_base_recal_1.table \
		-after "$output_var"/"${sample}"_base_recal_2.table \
		-plots "$output_var"/"${sample}"_recalibration_plot.pdf ;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collecting Allignment Summary Metrics#";
echo "#######################################";
echo -e "\n";

#Collecting Alignment metrics
java -jar "${picard_root}" CollectAlignmentSummaryMetrics \
		R="${genome_reference_directory}/${reference}" \
		I="$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam" \
		O="$output_var"/"${sample}"_allignment_metrics.txt;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collect Insert size metrics#";
echo "#######################################";
echo -e "\n";

#Collecting Insert Metrics
java -jar "${picard_root}" CollectInsertSizeMetrics \
		I="$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam" \
		O="$output_var"/"${sample}"_insert_metrics.txt \
		HISTOGRAM_FILE="$output_var"/"${sample}"_insert_metrics_HISTOGRAM.pdf;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Quality check of merged bam file - ""#";
echo "#######################################";
echo -e "\n";

"${qualimap_root}"./qualimap bamqc -bam "$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam" -outdir "${output_var}"/"$sample"_bamqc

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Removing intermeditary BAM files#";
echo "#######################################";
echo -e "\n";

final_bam="$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam"

# Check if BQSR done file is present or not
if [ ! -e "$final_bam" ]; then
    echo "$final_bam final file is not present."
else
    rm "$output_var"/"${sample}"_mapped.bam
	rm "$output_var"/"${sample}"_mapped_sorted.bam
	rm "$output_var"/"${sample}"_mapped_sorted_duprem.bam
	rm "$output_var"/"${sample}"_mapped_sorted_duprem_RG.bam
    echo "$final_bam is present."
	echo "Removing intermeditary BAM files"
fi

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Variant calling using Haplotypecaller#";
echo "#######################################";
echo -e "\n";

#haplotypecaller variant calling
"${gatk_root}"./gatk HaplotypeCaller  \
		-R "${genome_reference_directory}/${reference}" \
		-I "$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam"  \
		-O "$gatk_output"/"${sample}"_raw_variants.vcf  \
		--sample-name "${sample}" ;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Hard Filtration of Variants#";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_variants.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_variants.vcf  \
		-filter-name "QD_filter" -filter "QD < 2.00"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10.00" \
		-filter-name "FS_filter" -filter "FS >  60.00" \
		-filter-name "SOR_filter" -filter "SOR > 4.0" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 10";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered variants#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants  \
		--exclude-filtered  \
		-V "$gatk_output"/"${sample}"_Filtered_variants.vcf  \
		-O "$gatk_output"/"${sample}"_Final_variants.vcf ;

date +"%D %T":"##Process Completed Time##";
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding SNPs#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants    \
		--select-type-to-include SNP  \
		-V "$gatk_output"/"${sample}"_raw_variants.vcf  \
		-O "$gatk_output"/"${sample}"_raw_SNP.vcf ;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding INDELs#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants   \
		--select-type-to-include INDEL  \
		-V "$gatk_output"/"${sample}"_raw_variants.vcf  \
		-O "$gatk_output"/"${sample}"_raw_INDEL.vcf ;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Hard Filtration of Variants(SNP)#";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_SNP.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_SNP.vcf  \
		-filter-name "QD_filter" -filter "QD < 2.00"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10.00" \
		-filter-name "FS_filter" -filter "FS >  60.00" \
		-filter-name "SOR_filter" -filter "SOR > 4.0" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 10";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered variants(SNP)#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants  \
		--exclude-filtered  \
		-V "$gatk_output"/"${sample}"_Filtered_SNP.vcf  \
		-O "$gatk_output"/"${sample}"_Final_SNP.vcf ;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Hard Filtration of Variants(INDEL)#";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_INDEL.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_INDEL.vcf  \
		-filter-name "QD_filter" -filter "QD < 2.00"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10.00" \
		-filter-name "FS_filter" -filter "FS >  60.00" \
		-filter-name "SOR_filter" -filter "SOR > 4.0" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 10";
		

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered Variants(INDEL)#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants  \
		--exclude-filtered  \
		-V "$gatk_output"/"${sample}"_Filtered_INDEL.vcf \
		-O "$gatk_output"/"${sample}"_Final_INDEL.vcf ;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Variant Calling using strelka#";
echo "#######################################";
echo -e "\n";

#variant calling using strelka

source activate $python2_env


python "${strelka_root}"configureStrelkaGermlineWorkflow.py --bam "$output_var/${sample}_mapped_sorted_duprem_RG_BQSR.bam" --referenceFasta "${genome_reference_directory}/${reference}" --runDir "${strelka_output}"

echo "running workflow"

python "${strelka_output}"/runWorkflow.py -m local


conda deactivate

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Strelka Variant filtering using bcftools#";
echo "#######################################";
echo -e "\n";

#variant filtering using bcftools
bcftools view -i 'FILTER="PASS" && FORMAT/SB < 10' "${strelka_output}"/results/variants/variants.vcf.gz -o "${strelka_output}"/results/variants/"${sample}"_strelka_filtered_variants.vcf.gz



date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Intersection variants of gatk and strelka#";
echo "#######################################";
echo -e "\n";


"${bedtools_root}"./bedtools intersect -header -a "$gatk_output"/"${sample}"_Final_SNP.vcf -b "${strelka_output}"/results/variants/"${sample}"_strelka_filtered_variants.vcf.gz > "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Stats of the instersected variants for sample $sample #";
echo "#######################################";
echo -e "\n";

bcftools stats "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf > "${output_var}"/"$sample"_variants_stats.txt

date +"%D %T":"##Process Completed Time##";



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Annoattion with Funcotator Germline database (SNP) #";
echo "#######################################";
echo -e "\n";

#Annoattion with Funcotator Germline database (SNP)
"${gatk_root}"./gatk Funcotator \
		--variant "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf \
		--reference "${genome_reference_directory}/${reference}" \
		--ref-version hg38 \
		--data-sources-path $func_germline \
		--output "$output_anno"/"${sample}"_Germline_Funcotated_variants.vcf \
		--output-file-format VCF;

date +"%D %T":"##Process Completed Time##";
#COM
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Curating Germline Funcotation values (SNP) #";
echo "#######################################";
echo -e "\n";

#Curating Germline Funcotation values (SNP)
"${gatk_root}"./gatk VariantsToTable \
		-V "$output_anno"/"${sample}"_Germline_Funcotated_variants.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
		-O "$output_anno"/"${sample}"_germline_variants.table;

cat "$output_anno"/"${sample}"_Germline_Funcotated_variants.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' | sed -e 's/##INFO=<ID=FUNCOTATION.*: //' > "$output_anno"/"${sample}"_germline_curated_variants.txt
bcftools query -H -f '[%INFO/FUNCOTATION]\n' "$output_anno"/"${sample}"_Germline_Funcotated_variants.vcf | sed 's/|/\t/g' | awk '!/^#/ {print}' >> "$output_anno"/"${sample}"_germline_curated_variants.txt
		
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Filtering the required columns from funcotator germline annotation#";
echo "#######################################";
echo -e "\n";

awk -F'\t' -v cols="Gencode_43_hugoSymbol,Gencode_43_ncbiBuild,Gencode_43_chromosome,Gencode_43_start,Gencode_43_end,Gencode_43_variantClassification,Gencode_43_secondaryVariantClassification,Gencode_43_variantType,Gencode_43_refAllele,Gencode_43_tumorSeqAllele1,Gencode_43_tumorSeqAllele2,Gencode_43_genomeChange,Gencode_43_annotationTranscript,ACMGLMMLof_LOF_Mechanism,ACMGLMMLof_Mode_of_Inheritance,ACMGLMMLof_Notes,ACMG_recommendation_Disease_Name,ClinVar_VCF_AF_ESP,ClinVar_VCF_AF_EXAC,ClinVar_VCF_AF_TGP,ClinVar_VCF_ALLELEID,ClinVar_VCF_CLNDISDB,ClinVar_VCF_CLNDISDBINCL,ClinVar_VCF_CLNDN,ClinVar_VCF_CLNDNINCL,ClinVar_VCF_CLNHGVS,ClinVar_VCF_CLNREVSTAT,ClinVar_VCF_CLNSIG,ClinVar_VCF_CLNSIGCONF,ClinVar_VCF_CLNSIGINCL,ClinVar_VCF_CLNVC,ClinVar_VCF_CLNVCSO,ClinVar_VCF_CLNVI,ClinVar_VCF_DBVARID,ClinVar_VCF_GENEINFO,ClinVar_VCF_MC,ClinVar_VCF_ORIGIN,ClinVar_VCF_RS,ClinVar_VCF_ID,ClinVar_VCF_FILTER,gnomAD_exome_AF,gnomAD_genome_AF" '
    NR==1 {
        split(cols, arr, ",");
        for (i = 1; i <= NF; i++) {
            col_index[$i] = i;
        }
    }
    {
        for (i = 1; i <= length(arr); i++) {
            printf("%s%s", $col_index[arr[i]], (i == length(arr)) ? "\n" : "\t");
        }
    }
' "$output_anno"/"${sample}"_germline_curated_variants.txt > "$output_anno"/"${sample}"_germline_curated_filtered_variants.txt



date +"%D %T":"##Process Completed Time##";


: << 'COMMENT'
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Annoattion with Funcotator Somatic database (INDEL) #";
echo "#######################################";
echo -e "\n";

#Annoattion with Funcotator Somatic database (SNP)
"${gatk_root}"./gatk Funcotator \
		--variant "$gatk_output"/"${sample}"_Final_variants.vcf \
		--reference "${genome_reference_directory}/${reference}" \
		--ref-version hg38 \
		--data-sources-path $func_somatic \
		--output "$output_anno"/"${sample}"_Somatic_Funcotated_variants.vcf \
		--output-file-format VCF;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Curating Somatic Funcotation values (INDEL) #";
echo "#######################################";
echo -e "\n";

#Curating Somatic Funcotation values (SNP)
"${gatk_root}"./gatk VariantsToTable \
		-V "$output_anno"/"${sample}"_Somatic_Funcotated_variants.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
		-O "$output_anno"/"${sample}"_somatic_snp.table;

cat "$output_anno"/"${sample}"_Somatic_Funcotated_variants.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' > "$output_anno"/"${sample}"_somatic_curated_variants.txt
bcftools query -H -f '[%INFO/FUNCOTATION]\n' "$output_anno"/"${sample}"_Somatic_Funcotated_variants.vcf | sed 's/|/\t/g' >> "$output_anno"/"${sample}"_somatic_curated_variants.txt 
		
date +"%D %T":"##Process Completed Time##";

COMMENT



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Annotation with SNPEFF #";
echo "#######################################";
echo -e "\n";

java -jar "${snpeff_root}"snpEff.jar eff GRCh38.p14 "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf > "$output_anno"/"${sample}"_snpeff_anno_variants.vcf

cat "$output_anno"/"${sample}"_snpeff_anno_variants.vcf | "${snpeff_root}"scripts/vcfEffOnePerLine.pl | java -jar "${sift_root}" extractFields - CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C LOF[*].GENE LOF[*].NUMTR LOF[*].PERC > "$output_anno"/"${sample}"_snpeff_anno_filtered_variants.vcf

date +"%D %T":"##Process Completed Time##";
: <<'TAP'


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Annotation with Ensemble VEP#";
echo "#######################################";
echo -e "\n";

"${vep_root}"./vep -i "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf -o "$output_anno"/"${sample}"_vep_anno_variants.txt --database --species "human" --everything --tab --verbose

awk '!/^##/ {print}' "$output_anno"/"${sample}"_vep_anno_variants.txt > "$output_anno"/"${sample}"_vep_anno_curated_variants.txt

awk -F'\t' -v cols="#Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,\
cDNA_position,Codons,Existing_variation,IMPACT,STRAND,FLAGS,VARIANT_CLASS,SYMBOL,BIOTYPE,\
CANONICAL,TSL,APPRIS,CCDS,SIFT,PolyPhen,EXON,INTRON,HGVSc,AF,CLIN_SIG,SOMATIC,TRANSCRIPTION_FACTORS" '
    NR==1 {
        split(cols, arr, ",");
        for (i = 1; i <= NF; i++) {
            col_index[$i] = i;
        }
    }
    {
        for (i = 1; i <= length(arr); i++) {
            printf("%s%s", $col_index[arr[i]], (i == length(arr)) ? "\n" : "\t");
        }
    }
' "$output_anno"/"${sample}"_vep_anno_curated_variants.txt > "$output_anno"/"${sample}"_vep_anno_curated_filtered_variants.txt

#COM
date +"%D %T":"##Process Completed Time##";

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' "$gatk_output"/"${sample}"_Final_SNP.vcf > "$gatk_output"/"${sample}"_Final_SNP_INFO.vcf 
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\t%FORMAT\n'  "${strelka_output}"/results/variants/"${sample}"_strelka_filtered_variants.vcf.gz > "${strelka_output}"/results/variants/"${sample}"_strelka_variants_INFO.vcf
#COM

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#merging all the annotations in tab deliminated format#";
echo "#######################################";
echo -e "\n";

python "${code_root}"merging.py "$output_anno"/"${sample}"_vep_anno_curated_filtered_variants.txt "$output_anno"/"${sample}"_snpeff_anno_filtered_variants.vcf "$output_anno"/"${sample}"_germline_curated_filtered_variants.txt "$gatk_output"/"${sample}"_Final_SNP_INFO.vcf "${strelka_output}"/results/variants/"${sample}"_strelka_variants_INFO.vcf "$output_final"/"${sample}"_final_annotations_2.tsv

date +"%D %T":"##Process Completed Time##";
TAP
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#merging all the stats in tab deliminated format#";
echo "#######################################";
echo -e "\n";

# Define the strings you want to grep and their corresponding column headers for each file type
declare -A fastqc_grep_strings=(
    ["Total Sequences"]="Total_Sequences"
    ["Sequence length"]="Sequence_length"
    ["%GC"]="GC_percentage"
)

declare -A trimmo_grep_strings=(
    ["Total Sequences"]="Total_Trimmed_Sequences"
)

declare -A bamqc_grep_strings=(
    ["number of reads"]="Number_of_reads"
    ["number of mapped reads"]="Number_of_mapped_reads"
    ["mean mapping quality"]="Mean_mapping_quality"
    ["mean coverageData"]="Mean_coverageData"
)

declare -A variants_stats_grep_strings=(
    ["number of samples:"]="Number_of_samples"
    ["number of records:"]="Number_of_records"
    ["number of SNPs:"]="Number_of_SNPs"
    ["number of indels:"]="Number_of_indels"
)

# Function to print header
print_header() {
    printf "%s\t" "Sample"
    printf "%s\t" "Read1"
    printf "%s\t" "Read2"
    for key in "Total Sequences" "Sequence length" "%GC"; do
        printf "%s\t" "${fastqc_grep_strings[$key]}"
    done
    for key in "Total Sequences"; do
        printf "%s\t" "${trimmo_grep_strings[$key]}"
    done
    for key in "number of reads" "number of mapped reads" "mean mapping quality" "mean coverageData"; do
        printf "%s\t" "${bamqc_grep_strings[$key]}"
    done
    for key in "number of samples:" "number of records:" "number of SNPs:" "number of indels:"; do
        printf "%s\t" "${variants_stats_grep_strings[$key]}"
    done
    printf "\n"
}

# Function to print data for fastqc files
print_fastqc_data() {
    local file="$1"
    for key in "Total Sequences" "Sequence length" "%GC"; do
        value=$(grep "$key" "$file" | awk -F'\t' '{print $NF}')
        printf "%s\t" "$value"
    done
}

# Function to print data for trimmo files
print_trimmo_data() {
    local file="$1"
    for key in "Total Sequences"; do
        value=$(grep "$key" "$file" | awk -F'\t' '{print $NF}')
        printf "%s\t" "$value"
    done
}

# Function to print data for bamqc files
print_bamqc_data() {
    local file="$1"
    for key in "number of reads" "number of mapped reads" "mean mapping quality" "mean coverageData"; do
        value=$(grep "$key" "$file" | awk -F= '{gsub(/^[^=]+=/,"",$0); print $0}')
        printf "%s\t" "$value"
    done
}

# Function to print data for variants stats files
print_variants_stats() {
    local file="$1"
    for key in "number of samples:" "number of records:" "number of SNPs:" "number of indels:"; do
        value=$(grep "$key" "$file" | sed 's/^[^\t]*\t[^\t]*\t[^\t]*\t//')
        printf "%s\t" "$value"
    done
}

# Print header
print_header > "$output_final/stats.tsv"

# Print sample name
printf "%s\t" "$sample" >> "$output_final/stats.tsv"

# Print sample name
printf "%s\t" "$seq1" >> "$output_final/stats.tsv"
printf "%s\t" "$seq2" >> "$output_final/stats.tsv"

# Print data for fastqc files
print_fastqc_data "$output_seq/${seq1}_fastqc/fastqc_data.txt" >> "$output_final/stats.tsv"

# Print data for trimmo files
print_trimmo_data "$output_seq/${seq1}_paired_fastqc/fastqc_data.txt" >> "$output_final/stats.tsv"

# Print data for bamqc files
print_bamqc_data "${output_var}/${sample}_bamqc/genome_results.txt" >> "$output_final/stats.tsv"

# Print data for variants stats files
print_variants_stats "${output_var}/${sample}_variants_stats.txt" >> "$output_final/stats.tsv"


date +"%D %T":"##Process Completed Time##";
: <<'COPIUM'

COPIUM

echo done 

# Record the end time
end_time=$(date +%s)

# Calculate the total running time
total_time=$((end_time - start_time))

# Convert seconds to a human-readable format (hh:mm:ss)
formatted_time=$(date -u -d @$total_time +'%H:%M:%S')

# Print the total running time
echo "Total running time: $formatted_time"

exit 0;

















