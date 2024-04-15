#!/bin/bash

source /home/saurabh/Videos/Referance_data/code/config.sh

# Record the start time
start_time=$(date +%s)

sample="$1"
sample=$(basename "$sample")

read1="$2"
seq1=$(basename "$read1" .fastq)

read2="$3"
seq2=$(basename "$read2" .fastq)

bam="$4"

var="$5"

# Creating a directory based on sequence names
output_var="${sample}_variants_output"
output_seq="${sample}_sequences"
gatk_output="${output_var}/gatk_output"
strelka_output="${output_var}/strelka_output"
output_anno="${sample}_annotation_output"
output_final="${sample}_final_output"

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
print_header > "$output_final/stats.txt"

# Print sample name
printf "%s\t" "$sample" >> "$output_final/stats.txt"

# Print data for fastqc files
print_fastqc_data "$output_seq/${seq1}_fastqc/fastqc_data.txt" >> "$output_final/stats.txt"

# Print data for trimmo files
print_trimmo_data "$output_seq/${seq1}_paired_fastqc/fastqc_data.txt" >> "$output_final/stats.txt"

# Print data for bamqc files
print_bamqc_data "${output_var}/${sample}_bamqc/genome_results.txt" >> "$output_final/stats.txt"

# Print data for variants stats files
print_variants_stats "${output_var}/${sample}_variants_stats.txt" >> "$output_final/stats.txt"
