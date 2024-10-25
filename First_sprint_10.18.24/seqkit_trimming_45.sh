#!/bin/bash
#SBATCH --job-name=seqkit_removal.sh
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/seqkit_removal.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/seqkit_removal_error.log
#SBATCH --nodes=1                 # Requesting 1 node
#SBATCH --ntasks=16               # Using 16 CPU cores
#SBATCH --mem=32G                 # Adjust memory as needed
#SBATCH --time=10:00:00           # Adjust time as needed

# Load seqkit
module load seqkit/2.3.0

# Input file list
file_list="/labs/yunroseli_grp/Jon/32_samples_10182024/code/mapping_files/list_of_trimmed.txt"
output_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed_removedempty_nless45bp/"

# Define your processing function
process_fastq() {
    input_file="$1"
    output_file="$output_dir$(basename "$input_file")"
    
    # Check if the output file already exists to avoid duplicate processing
    if [ -f "$output_file" ]; then
        echo "$output_file already exists, skipping."
        return 0
    fi
    
    # Perform the sequence trimming
    seqkit seq -m 45 "$input_file" > "$output_file"
}

# Loop through each file in the file list and process them one by one
while IFS= read -r input_file; do
    echo "Processing: $input_file"
    process_fastq "$input_file"
done < "$file_list"
