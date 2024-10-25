#!/bin/bash

# Define your input and output directories
file_locations="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed"
cut_output="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed_no_empty"

# Loop through all the fastq.gz files in the input directory
for input_file in "$file_locations"/*.fastq.gz; do
    # Get the basename (just the file name without the directory)
    base_name=$(basename "$input_file")
    
    # Define the output file path
    output_file="$cut_output/$base_name"
    
    # Run cutadapt to filter sequences and keep only those with 90 bp or longer
    cutadapt -m 90 -o "$output_file" --cores=25 "$input_file"
    
    
    echo "Processed: $input_file -> $output_file"
done
