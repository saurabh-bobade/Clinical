
import pandas as pd
import sys

# Function to process a file
def process_file(file_path):
    
    # Load the annotations
    annotations = pd.read_csv(file_path, sep='\t')

    return annotations

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python merging.py vep_annotations snpeff_annotations funcotator_annotations gatk_info_file strelka_info_file output_path")
        sys.exit(1)

    # Extract file paths and output path from command-line arguments
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    file3_path = sys.argv[3]
    file4_path = sys.argv[4]
    file5_path = sys.argv[5]
    output_path = sys.argv[6]

    # Process each file
    vep = process_file(file1_path)
    snpeff = process_file(file2_path)
    funcotator = process_file(file3_path)
    gatk = process_file(file4_path)
    strelka = process_file(file5_path)

    # Identify the variant column for vep annotations
    vep_column = vep.iloc[:, 0]  # You may need to adjust this index based on your data

    # Processing vep annotations
    first_column = vep_column.str.rsplit('_', n=1).str[0]
    chrom_pos_split = first_column.str.rsplit('_', n=1).str
    vep['CHROM'] = chrom_pos_split[0]
    vep['POS'] = chrom_pos_split[1].astype(int)
    ref_alt_split = vep_column.str.rsplit('_', n=1).str[1].str.split('/', expand=True)
    vep['REF'] = ref_alt_split[0]
    vep['ALT'] = ref_alt_split[1]
    #Rearranging all the columns in vep
    column_order = ['CHROM', 'POS', 'REF', 'ALT'] + [col for col in vep.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]
    vep = vep[column_order]

    # Processing funcotator annotations
    funcotator['CHROM'] = funcotator['Gencode_43_chromosome']
    funcotator['POS'] = funcotator['Gencode_43_start']
    funcotator['REF'] = funcotator['Gencode_43_refAllele']
    funcotator['ALT'] = funcotator['Gencode_43_tumorSeqAllele2']
    
    #Rearranging all the columns in funcotator
    column_order = ['CHROM', 'POS', 'REF', 'ALT'] + [col for col in funcotator.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]
    funcotator = funcotator[column_order]
    
    # Processing gatk info
    #Changing the column names
    gatk = gatk.rename(columns={
    '# [1]CHROM': 'CHROM',
    '[2]POS': 'POS',
    '[3]REF': 'REF',
    '[4]ALT': 'ALT',
    '[5](null)': 'GATK_INFO'
    })
    
    # Processing strelka info
    #Making new column for chromosome info
    strelka['CHROM'] = strelka.index
    #Changing the column names
    strelka = strelka.rename(columns={
    '# [1]CHROM': 'POS',
    '[2]POS': 'REF',
    '[3]REF': 'ALT',
    '[4]ALT': 'Strelka_info',
    '[5](null)': 'Strelka_format1',
    '[6](null)': 'Strelka_format2'
    })
    
    #Rearranging all the columns in strelka INFO
    column_order = ['CHROM', 'POS', 'REF', 'ALT'] + [col for col in strelka.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]
    strelka = strelka[column_order]

    # Dropping duplicates from vep and snpeff annotations (consider the behavior)
    vep = vep.drop_duplicates()
    snpeff = snpeff.drop_duplicates()

    # Merging vep annotations and snpeff annotations based on common columns
    common_columns = ['CHROM', 'POS', 'REF', 'ALT']
    merged_1 = pd.merge(gatk, strelka, on=common_columns, how='inner')
    
    # Merging merged_1 and vep annotation
    merged_2 = pd.merge(merged_1, vep, on=common_columns, how='inner')
    
    # Merging merged_2 and snpeff annotation
    merged_3 = pd.merge(merged_2, snpeff, on=common_columns, how='inner')
    
    # Merging merged_3 and funcotated annotation
    merged_df = pd.merge(merged_3, funcotator, on=common_columns, how='inner')

    # Output the result in tab-delimited format to the specified output path
    merged_df.to_csv(output_path, sep='\t', index=False)
    
    
    
    
    
    
    
    
