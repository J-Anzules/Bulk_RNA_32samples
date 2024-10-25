#!/bin/bash
#SBATCH --job-name=Final_trim    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Final_trim.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Final_trim_error.log  # Error log
#SBATCH -N 1                         # Use 1 node
#SBATCH --cpus-per-task=32            # Request 72 CPUs
#SBATCH --mem=64G                   # Request 1.5 TB memory
#SBATCH --time=24:00:00               # Set a time limit of 24 hours
#SBATCH -p all                        # Specify the partition (if applicable)

# Define your input and output directories
file_locations="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed/"
cut_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed3/"

# Load the cutadapt module (if needed)
module load cutadapt  # Adjust if you use modules

# Loop through all the R1 files in the input directory
for r1_file in "$file_locations"/*_R1.fastq.gz; do
    # Define the corresponding R2 file
    r2_file="${r1_file/_R1/_R2}"
    
    # Check if R2 file exists
    if [ -f "$r2_file" ]; then
        # Get the basenames (just the file names without the directory)
        r1_base=$(basename "$r1_file")
        r2_base=$(basename "$r2_file")
        
        # Define the output file paths for R1 and R2
        output_r1="$cut_output/$r1_base"
        output_r2="$cut_output/$r2_base"
        
        # Run cutadapt for paired-end reads and filter sequences with 90 bp or longer, using 72 cores
        cutadapt -m 90 --paired-output "$output_r2" -o "$output_r1" --cores=50 "$r1_file" "$r2_file"
        
        echo "Processed: $r1_file and $r2_file -> $output_r1 and $output_r2"
    else
        echo "Warning: Skipping $r1_file as the corresponding R2 file was not found."
    fi
done
