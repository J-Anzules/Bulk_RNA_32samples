#!/bin/bash

# Define the file containing the list of fastq files
file_list="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed/fastq_files_list.txt"

# Define the directory where you want to save the FastQC reports
output_dir="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed/qcreport"

# Ensure the output directory exists, if not, create it
mkdir -p "$output_dir"

# Loop through each fastq file in the list and run FastQC
while IFS= read -r fastq_file; do
    if [ -f "$fastq_file" ]; then
        echo "Processing $fastq_file with FastQC..."
        fastqc "$fastq_file" -o "$output_dir"
    else
        echo "File $fastq_file does not exist."
    fi
done < "$file_list"

echo "FastQC reports have been generated and saved to $output_dir."
