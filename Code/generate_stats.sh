#!/bin/bash

# Check if the input file is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <config_file>"
    exit 1
fi

source "$2"

# Record the start time
start_time=$(date +%s)


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
    printf "%s\t" "Read 1"
    printf "%s\t" "Read 2"
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

output_final="$output_root/Final_statistics"
mkdir -p $output_final
# Print header
print_header > "$output_final/stats.tsv"

# Loop through each sample in the sample sheet, skipping the header
tail -n +2 "$1" | while IFS= read -r sample_sheet_line; do
    # Extract sample name and read paths
    sample=$(echo "$sample_sheet_line" | cut -d$'\t' -f1)
    read1=$(echo "$sample_sheet_line" | cut -d$'\t' -f2)
    read2=$(echo "$sample_sheet_line" | cut -d$'\t' -f3)
    seq1=$(basename "$read1" .fastq)
    seq2=$(basename "$read2" .fastq)


    # Creating a directory based on sequence names
    output_var="${output_root}/${sample}_variants_output"
    output_seq="${output_root}/${sample}_sequences"
    gatk_output="${output_var}/gatk_output"
    strelka_output="${output_var}/strelka_output"
    output_anno="${output_root}/${sample}_annotation_output"
    
    # Specify paths to Trimmo, BAMQC, and variants stats files for the current sample
    fastqc_file="$output_seq/${seq1}_fastqc/fastqc_data.txt"
    trimmo_file="$output_seq/${seq1}_paired_fastqc/fastqc_data.txt"
    bamqc_file="$output_var/${sample}_bamqc/genome_results.txt"
    variants_stats_file="$output_var/${sample}_variants_stats.txt"


    # Print sample name
    printf "%s\t" "$sample" >> "$output_final/stats.tsv"
    printf "%s\t" "$seq1" >> "$output_final/stats.tsv"
    printf "%s\t" "$seq2" >> "$output_final/stats.tsv"

    # Print data for fastqc files
    print_fastqc_data "$fastqc_file" >> "$output_final/stats.tsv"

    # Print data for trimmo files
    print_trimmo_data "$trimmo_file" >> "$output_final/stats.tsv"

    # Print data for bamqc files
    print_bamqc_data "$bamqc_file" >> "$output_final/stats.tsv"

    # Print data for variants stats files
    print_variants_stats "$variants_stats_file" >> "$output_final/stats.tsv"
    
    # Print next line in stats.tsv
    printf "\n" >> "$output_final/stats.tsv"

done

