#!/bin/bash
#SBATCH --job-name=trim6    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/6trim.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/6trim_error.log
#SBATCH --nodelist=ppxhpcnode22 # Request this specific node
#SBATCH --ntasks=1              # Run a single task
#SBATCH --cpus-per-task=28      # Request 26 CPUs (available)
#SBATCH --mem=230G              # Request 213 GB of memory (available)
#SBATCH --time=24:00:00         # Set an appropriate time limit


# Define file locations
fastqz_list="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/code/mapping_files/fastq_addresses.txt"
cut_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trim6"

# Load the cutadapt module
module load pigz
module load cutadapt

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
            output_r1="$cut_output/${anumber}_R1.fastq.gz"
            output_r2="$cut_output/${anumber}_R2.fastq.gz"

            # Run cutadapt for paired-end reads
            cutadapt \
                -m 100 \
                -o "$output_r1" -p "$output_r2" \
                -j 25 \
                "$r1_file" "$r2_file"

            echo "Successfully processed: $r1_file and $r2_file -> $output_r1 and $output_r2"
        else
            echo "Warning: R1 and/or R2 files not found for $line"
        fi
    fi
done < "$fastqz_list"
