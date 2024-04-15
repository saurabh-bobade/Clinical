#!/bin/bash


# Record the start time
start_time=$(date +%s)


# Check if the required number of arguments are provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <sample> <read1> <read2> <lane> <Platform> <config_file>"
    exit 1
fi

source $6

# Assign command-line arguments to variables
sample="$1"
sample=$(basename "$sample")

read1="$2"
seq1=$(basename "$read1" .fastq)

read2="$3"
seq2=$(basename "$read2" .fastq)

lane="$4"
lane=$(basename "$lane")

platform="$5"

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
mkdir -p $final_stats

: << 'COMMENT'
echo -e "Commented Block"
COMMENT


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check for $seq1 of $sample from $lane #";
echo "#######################################";
echo -e "\n";
fastqc --outdir $output_seq $read1
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check for $seq2 of $sample from $lane #";
echo "#######################################";
echo -e "\n";
fastqc --outdir $output_seq $read2
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
echo "#FASTA quality check of "${seq1}"_paired.fastq.gz of $sample from $lane #";
echo "#######################################";
echo -e "\n";
fastqc "$output_seq"/"${seq1}"_paired.fastq.gz
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check of "${seq2}"_paired.fastq.gz of $sample from $lane #";
echo "#######################################";
echo -e "\n";
fastqc "$output_seq"/"${seq2}"_paired.fastq.gz
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
echo "#mapping/alignment of reads with reference genome for $sample from $lane #";
echo "#######################################";
echo -e "\n";


"${bwa_root}"bwa mem -M -t $threads "${genome_reference_directory}/${reference}" "$output_seq"/"${seq1}"_paired.fastq.gz "$output_seq"/"${seq2}"_paired.fastq.gz > "$output_var"/"${sample}"_"${lane}"_mapped.sam
#bwa mem -M -t 4 "$reference" "$read1" "$read2" > "$output_var"/"${sample}"_"${lane}"_mapped.sam
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#convert sam to bam for $sample from $lane #";
echo "#######################################";
echo -e "\n";

samtools view -b "$output_var"/"${sample}"_"${lane}"_mapped.sam -o "$output_var"/"${sample}"_"${lane}"_mapped.bam

date +"%D %T":"##Process Completed Time##";

 
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#sorting bam file for $sample from $lane #";
echo "#######################################";
echo -e "\n";

samtools sort "$output_var"/"${sample}"_"${lane}"_mapped.bam -o "$output_var"/"${sample}"_"${lane}"_mapped_sorted.bam


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#rMarking duplicates for $sample from $lane #";
echo "#######################################";
echo -e "\n";

java -Xmx"$cpu"g -jar "${picard_root}" MarkDuplicates \
    INPUT="$output_var"/"${sample}"_"${lane}"_mapped_sorted.bam \
    OUTPUT="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem.bam \
    METRICS_FILE="$output_var/${sample}_duplicate_metrics.txt" \
    REMOVE_DUPLICATES=false 


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#add read groups to bam file for $sample from $lane #";
echo "#######################################";
echo -e "\n";

java -Xmx"$cpu"g -jar "${picard_root}" AddOrReplaceReadGroups \
		I="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem.bam \
		O="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG.bam \
		RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM="$sample";


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing final bam file for $sample from $lane #";
echo "#######################################";
echo -e "\n";

samtools index "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG.bam

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Displaying stats of BAM file for $sample from $lane #";
echo "#######################################";
echo -e "\n";

samtools stats "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG.bam > "$output_var"/"${sample}"_"${lane}"_stats.txt

date +"%D %T":"##Process Completed Time##";



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 1 {BQSR} for $sample from $lane #";
echo "#######################################";
echo -e "\n";

#Creating base recal table

"${gatk_root}"./gatk BaseRecalibrator \
		-I "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG.bam \
		-R "${genome_reference_directory}/${reference}" \
		--known-sites "${variants_reference_directory}/${known_variants}" \
		-O "$output_var"/"${sample}"_"${lane}"_base_recal_1.table

date +"%D %T":"##Process Completed Time##";



#Apply BQSR
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#apply BQSR for $sample #";
echo "#######################################";
echo -e "\n";

"${gatk_root}"./gatk ApplyBQSR \
		-R "${genome_reference_directory}/${reference}" \
		-I "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG.bam \
		--bqsr-recal-file "$output_var"/"${sample}"_"${lane}"_base_recal_1.table \
		-O "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG_BQSR.bam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 2 {BQSR} for $sample #";
echo "#######################################";
echo -e "\n";

#Creating base recal table 2 for covarities

"${gatk_root}"./gatk BaseRecalibrator \
		-I "$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG_BQSR.bam \
		-R "${genome_reference_directory}/${reference}" \
		--known-sites "${variants_reference_directory}/${known_variants}" \
		-O "$output_var"/"${sample}"_"${lane}"_base_recal_2.table;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Analyzing Covariaties for $sample #";
echo "#######################################";
echo -e "\n";

#Analyzing Covarities

"${gatk_root}"./gatk AnalyzeCovariates \
		-before "$output_var"/"${sample}"_"${lane}"_base_recal_1.table \
		-after "$output_var"/"${sample}"_"${lane}"_base_recal_2.table \
		-plots "$output_var"/"${sample}"_"${lane}"_recalibration_plot.pdf ;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collecting Allignment Summary Metrics for $sample #";
echo "#######################################";
echo -e "\n";
#Collecting Alignment metrics
java -Xmx"$cpu"g -jar "${picard_root}" CollectAlignmentSummaryMetrics \
		R=$"${genome_reference_directory}/${reference}" \
		I="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG_BQSR.bam \
		O="$output_var"/"${sample}"_"${lane}"_allignment_metrics.txt;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collect Insert size metrics for $sample #";
echo "#######################################";
echo -e "\n";
#Collecting Insert Metrics
java -Xmx"$cpu"g -jar "${picard_root}" CollectInsertSizeMetrics \
		I="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG_BQSR.bam \
		O="$output_var"/"${sample}"_"${lane}"_insert_metrics.txt \
		HISTOGRAM_FILE="$output_var"/"${sample}"_"${lane}"_insert_metrics_HISTOGRAM.pdf;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Removing intermeditary BAM files#";
echo "#######################################";
echo -e "\n";

final_bam="$output_var"/"${sample}"_"${lane}"_mapped_sorted_duprem_RG_BQSR.bam

# Check if BQSR done file is present or not
if [ ! -e "$final_bam" ]; then
    echo "$final_bam final file is not present."
else
    rm "$output_var"/"${sample}"_${lane}_mapped.bam
	rm "$output_var"/"${sample}"_${lane}_mapped_sorted.bam
	rm "$output_var"/"${sample}"_${lane}_mapped_sorted_duprem.bam
	rm "$output_var"/"${sample}"_${lane}_mapped_sorted_duprem_RG.bam
    echo "$final_bam is present."
	echo "Removing intermeditary BAM files"
fi

date +"%D %T":"##Process Completed Time##";

# After generating the BAM file, print its path
echo "$output_var/${sample}_${lane}_mapped_sorted_duprem_RG_BQSR.bam"

echo done 

# Record the end time
end_time=$(date +%s)

# Calculate the total running time
total_time=$((end_time - start_time))

# Convert seconds to a human-readable format {hh:mm:ss}
formatted_time=$(date -u -d "@$total_time" +'%H:%M:%S')

# Print the total running time
echo "Total running time: $formatted_time"

exit 0;


















