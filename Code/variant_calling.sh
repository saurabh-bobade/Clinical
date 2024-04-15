#!/bin/bash


# Record the start time
start_time=$(date +%s)


# Check if the required number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sample> <lane> <Platform> <config_file>"
    exit 1
fi

source "$4"

# Assign command-line arguments to variables
sample="$1"

#bam="$2"

lane="$2"

platform="$3"

# Creating a directory based on sequence names
output_var="${output_root}/${sample}_variants_output"
output_seq="${output_root}/${sample}_sequences"
gatk_output="${output_var}/gatk_output"
strelka_output="${output_var}/strelka_output"
output_anno="${output_root}/${sample}_annotation_output"
output_final="${output_root}/${sample}_final_output"


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

# Merging BAM files for the current sample
bam_files=("${output_var}/${sample}"_*"_mapped_sorted_duprem_RG_BQSR.bam")
merged_bam="${output_var}/${sample}_merged_sorted_duprem_RG_BQSR.bam"

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Merging BAM files for $sample #";
echo "#######################################";
samtools merge -f "$merged_bam" "${bam_files[@]}"

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Indexing Merged BAM file for $sample #";
echo "#######################################";
samtools index "$merged_bam"
echo "Merged BAM file for $sample: $merged_bam"


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Quality check of merged bam file - ""#";
echo "#######################################";
echo -e "\n";

"${qualimap_root}"./qualimap bamqc -bam $merged_bam -outdir "${output_var}"/"$sample"_bamqc

date +"%D %T":"##Process Completed Time##";



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Variant calling using Haplotypecaller for $sample #";
echo "#######################################";
echo -e "\n";

#haplotypecaller variant calling
"${gatk_root}"./gatk HaplotypeCaller  \
		-R "${genome_reference_directory}/${reference}" \
		-I $merged_bam  \
		-O "$gatk_output"/"${sample}"_raw_variants.vcf  \
		--sample-name "$sample" ;


date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Hard Filtration of Variants for $sample #";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_variants.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_variants.vcf  \
		-filter-name "QD_filter" -filter "QD < $QD"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < $DP" \
		-filter-name "FS_filter" -filter "FS >  $FS" \
		-filter-name "SOR_filter" -filter "SOR > $SOR" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < $GQ";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered variants for $sample #";
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
echo "#Excluding SNPs for $sample #";
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
echo "#Excluding INDELs for $sample #";
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
echo "#Hard Filtration of Variants(SNP) for $sample #";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_SNP.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_SNP.vcf  \
		-filter-name "QD_filter" -filter "QD < $QD"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < $DP" \
		-filter-name "FS_filter" -filter "FS >  $FS" \
		-filter-name "SOR_filter" -filter "SOR > $SOR" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < $GQ";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered variants(SNP) for $sample #";
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
echo "#Hard Filtration of Variants(INDEL) for $sample #";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  \
		-R "${genome_reference_directory}/${reference}" \
		-V "$gatk_output"/"${sample}"_raw_INDEL.vcf  \
		-O "$gatk_output"/"${sample}"_Filtered_INDEL.vcf  \
		-filter-name "QD_filter" -filter "QD < $QD"  \
		-genotype-filter-name "DP_filter" -genotype-filter-expression "DP < $DP" \
		-filter-name "FS_filter" -filter "FS >  $FS" \
		-filter-name "SOR_filter" -filter "SOR > $SOR" \
		-genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < $GQ";
		

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered Variants(INDEL) for $sample #";
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
echo "#Variant Calling using strelka for $sample #";
echo "#######################################";
echo -e "\n";

#variant calling using strelka

source activate $python2_env


python "${strelka_root}"configureStrelkaGermlineWorkflow.py --bam $merged_bam --referenceFasta "${genome_reference_directory}/${reference}" --runDir "${strelka_output}"

echo "running workflow"

python "${strelka_output}"/runWorkflow.py -m local


source deactivate

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Strelka Variant filtering using bcftools for $sample #";
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
echo "#Annoattion with Funcotator Germline database (SNP) for $sample #";
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
echo "#Curating Germline Funcotation values (SNP) for $sample #";
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
echo "#Filtering the required columns from funcotator germline annotation for $sample #";
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
		--reference $reference \
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
echo "#Annotation with SNPEFF for $sample #";
echo "#######################################";
echo -e "\n";

java -Xmx"$cpu"g -jar "${snpeff_root}"snpEff.jar eff GRCh38.p14 "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf > "$output_anno"/"${sample}"_snpeff_anno_variants.vcf

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Filtering out SNPEFF annotations for $sample #";
echo "#######################################";
echo -e "\n";

cat "$output_anno"/"${sample}"_snpeff_anno_variants.vcf | "${snpeff_root}"scripts/vcfEffOnePerLine.pl | java -Xmx"$cpu"g -jar "${sift_root}" extractFields - CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C LOF[*].GENE LOF[*].NUMTR LOF[*].PERC > "$output_anno"/"${sample}"_snpeff_anno_filtered_variants.vcf

date +"%D %T":"##Process Completed Time##";

: << 'SUP'
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Annotation with Ensemble VEP for $sample #";
echo "#######################################";
echo -e "\n";

"${vep_root}"./vep -i "$output_var"/"${sample}"_intersected_variants_gatk_and_strelka.vcf -o "$output_anno"/"${sample}"_vep_anno_variants.txt --database --species "human" --everything --tab --verbose

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Filtering out Ensemble VEP annotations for $sample #";
echo "#######################################";
echo -e "\n";

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
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Extracting INFO from gatk and strelka variants $sample #";
echo "#######################################";
echo -e "\n";

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' "$gatk_output"/"${sample}"_Final_variants.vcf > "$gatk_output"/"${sample}"_Final_variants_INFO.vcf 
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\t%FORMAT\n'  "${strelka_output}"/results/variants/"${sample}"_strelka_filtered_variants.vcf.gz > "${strelka_output}"/results/variants/"${sample}"_strelka_variants_INFO.vcf

date +"%D %T":"##Process Completed Time##";
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#merging all the annotations in tab deliminated format for $sample #";
echo "#######################################";
echo -e "\n";

python "${code_root}"merging.py "$output_anno"/"${sample}"_vep_anno_curated_filtered_variants.txt "$output_anno"/"${sample}"_snpeff_anno_filtered_variants.vcf "$output_anno"/"${sample}"_germline_curated_filtered_variants.txt "$gatk_output"/"${sample}"_Final_variants_INFO.vcf "${strelka_output}"/results/variants/"${sample}"_strelka_variants_INFO.vcf "$output_anno"/"${sample}"_final_annotations.txt


SUP
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#All process Done for $sample #";
echo "#######################################";
echo -e "\n"

# Record the end time
end_time=$(date +%s)

# Calculate the total running time
total_time=$((end_time - start_time))

# Convert seconds to a human-readable format (hh:mm:ss)
formatted_time=$(date -u -d @$total_time +'%H:%M:%S')

# Print the total running time
echo "Total running time: $formatted_time"

exit 0;


