#!/bin/bash
# Find names from file

original_file="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/20241015_LH00295_0140_B22VL37LT3/YRLP_0001/YRLP_0001_1_UN_Whole_C1_UUUUU_A13487_22VL37LT3_CGTCCATGTA_L007_R2_001.fastq.gz"


# Creating new name and output
new_name=$(echo "$original_file" | sed -E 's/.*_(A[0-9]+)_.*(R[12]).*/\1_\2.fastq.gz/')
cutadapt_output="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trim_test/$new_name"



# Trim
cutadapt -u 15 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\
    -a TGGAATTCTCGGGTGCCAAGG \
    -a "G{50}" \
    -o $cutadapt_output $original_file


# Run FastQC to generate the report and save it in the trim_test folder
fastqc $cutadapt_output -o ../data/trim_test/

