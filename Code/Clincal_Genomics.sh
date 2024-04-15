#!/bin/bash


#Configuration file path
config_root="/home/pdx5/Clinical_Genomics/Code/config.sh";
source $config_root

# Record the start time
start_time=$(date +%s)

# Check if clinical_genomics.sh is present in the same directory
if [ ! -f "${code_root}preprocessing_samples.sh" ]; then
    echo "Error: 'preprocessing_samples.sh' script is missing."
    exit 1
fi

# Check if variant_calling.sh is present in the same directory
if [ ! -f "${code_root}variant_calling.sh" ]; then
    echo "Error: 'variant_calling.sh' script is missing."
    exit 1
fi

# Check if samples_stats.sh is present in the same directory
if [ ! -f "${code_root}generate_stats.sh" ]; then
    echo "Error: 'generate_stats.sh' script is missing."
    exit 1
fi

# Check if samples.txt file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <samples.txt>"
    exit 1
fi

# Read the sample sheet, skipping the header
sample_sheet="$1"
{
  read  # read and discard the header
  while IFS=$'\t' read -r sample read1 read2 lane platform; do
      # Check if any required fields are missing
      if [ -z "$sample" ] || [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$lane" ] || [ -z "$platform" ]; then
          echo "Skipping line with missing fields."
          continue
      fi

      # Check if read1 and read2 files exist
      if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
          echo "Skipping '$read1' or '$read2' which doesn't exist or couldn't be read."
          continue
      fi

      # Output sample name and read file paths
      echo "Sample: $sample"
      echo "Read 1: $read1"
      echo "Read 2: $read2"
      echo "Lane: $lane"
      echo "Platform: $platform"
      echo "Reference: $reference"
      echo "Variants: $known_variants"
      echo "=========================================="

      # Process each lane for the current sample
      IFS=',' read -r -a lanes_array <<< "$lane"
      for lane_num in "${lanes_array[@]}"; do
          # Run clinical_genomics.sh for the current lane
          echo "Calling BAM files for each lane by their sample ids"
          "${code_root}"preprocessing_samples.sh "$sample" "$read1" "$read2" "$lane_num" "$platform" "$config_root"
      done

      # Output directory for sample BAM files
      output_var="${output_root}/${sample}_variants_output"

  done
} < "$sample_sheet"

# Merging BAM files and variant calling after processing all samples and lanes
processed_samples=()  # Array to store processed samples
{
  read  # read and discard the header
  while IFS=$'\t' read -r sample read1 read2 lane platform reference variants; do
      # Check if the sample has already been processed
      if [[ " ${processed_samples[@]} " =~ " ${sample} " ]]; then
          continue
      fi

      # Run variant_calling.sh with the merged BAM file
      "${code_root}"variant_calling.sh "$sample" "$lane" "$platform" "$config_root"

      # Add the processed sample to the array
      processed_samples+=("$sample")

  done
} < "$sample_sheet"

echo -e "\n\n\n";
echo -e "\n";
echo "#######################################";
echo "#MultiQC report generation on all files#";
echo "#######################################";

multiqc "$output_root"/ -e picard -e gatk -e vep -e samtools --force --outdir "$final_stats"/
echo -e "\n";

echo -e "\n\n\n";
echo -e "\n";
echo "#######################################";
echo "#gathering all the statistics in a TSV file#";
echo "#######################################";

"${code_root}"generate_stats.sh $sample_sheet "$config_root"

echo "Variant Calling and Annotation Completed for all samples"


# Record the end time
end_time=$(date +%s)

# Calculate the total running time
total_time=$((end_time - start_time))

# Convert seconds to a human-readable format {hh:mm:ss}
formatted_time=$(date -u -d "@$total_time" +'%H:%M:%S')

# Print the total running time
echo "Total running time: $formatted_time"
