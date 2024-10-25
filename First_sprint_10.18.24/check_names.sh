#!/bin/bash

# Define the file containing the paths
input_file="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/20241015_LH00295_0140_B22VL37LT3/fastq_files.txt"

# Declare an array to store the new names
declare -A new_names

# Loop through each line (file path) in the input file
while IFS= read -r line; do
    # Extract the A* number and R1/R2 info using sed
    new_name=$(echo "$line" | sed -E 's/.*_(A[0-9]+)_.*(R[12]).*/\1_\2.fastq.gz/')
    
    # Check if the new name is already in the array
    if [[ -n "${new_names[$new_name]}" ]]; then
        echo "Duplicate name detected: $new_name"
    else
        # If not a duplicate, add it to the array
        new_names[$new_name]=1
    fi
done < "$input_file"

echo "Check complete. If no duplicates are shown, the names are unique."
