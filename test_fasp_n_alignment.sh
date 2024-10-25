#!/bin/bash
#SBATCH --job-name=trim5    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/fastp.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/fastp_error.log
#SBATCH --ntasks=1              # Run a single task
#SBATCH --cpus-per-task=27      # Request all CPUs
#SBATCH --mem=200G              # Request 200 GB of memory
#SBATCH --time=24:00:00         # Set an appropriate time limit

# Define file locations
fastqz_list="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/code/mapping_files/fastq_addresses.txt"
fastp_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/fastp"

# Load the necessary modules
module load fastp

# Read through the list of fastq files and process them in pairs
while IFS= read -r line; do
    # Check if the file is an R1 file (assuming R1 and R2 convention is used)
    if [[ "$line" == *"_R1_"* ]]; then
        r1_file="$line"
        r2_file="${line/_R1_/_R2_}"  # Infer the corresponding R2 file

        # Check if both R1 and R2 files exist
        if [[ -f "$r1_file" && -f "$r2_file" ]]; then
            echo " "
            echo " "
            echo "-------------------------"
            echo "Processing R1: $r1_file"
            echo "Processing R2: $r2_file"
            echo "-------------------------"
            echo " "
            echo " "

            # Extract the A number from the filename for output file naming
            anumber=$(echo "$r1_file" | sed -E 's/.*_(A[0-9]+)_.*R1_.*/\1/')

            # Define output filenames based on the extracted A number
            output_r1="$fastp_output/${anumber}_R1.fastq.gz"
            output_r2="$fastp_output/${anumber}_R2.fastq.gz"

            # Run fastp with default settings for paired-end reads
            fastp \
                --in1 "$r1_file" --in2 "$r2_file" \
                --out1 "$output_r1" --out2 "$output_r2" \
                --thread 27

            echo "Successfully processed: $r1_file and $r2_file -> $output_r1 and $output_r2"
        else
            echo "Warning: R1 and/or R2 files not found for $line"
        fi
    fi
done < "$fastqz_list"


