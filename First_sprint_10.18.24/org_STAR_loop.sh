#!/bin/bash
#SBATCH --job-name=STAR_loop    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_loop.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_loop_error.log  # Error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=24:00:00               # Set an appropriate time limit
#SBATCH --mem=245G                    # Adjust memory as necessary
#SBATCH --cpus-per-task=30            # 30 CPUs

# Define file locations
fastqz_list="/labs/yunroseli_grp/Jon/32_samples_10182024/code/mapping_files/fastq_addresses.txt"
output_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results"
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"

# Load STAR module
module load STAR/2.5.3a-foss-2017a

# Create a function to extract the A number and process paired-end files
process_star_alignment() {
    r1_file="$1"
    r2_file="${r1_file/_R1_/_R2_}"

    # Extract the A number from the filename
    anumber=$(echo "$r1_file" | sed -E 's/.*_(A[0-9]+)_.*R1_.*/\1/')

    # Define temporary directory for each sample
    temp_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp_${anumber}"

    # Check if both R1 and R2 files exist
    if [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
        echo "Processing: $anumber"
        star_prefix="$anumber"

        # Run STAR alignment
        STAR --genomeDir "$genome_directory" \
            --readFilesCommand zcat \
            --readFilesIn "$r1_file" "$r2_file" \
            --outTmpDir "$temp_dir" \
            --outSAMtype BAM Unsorted \
            --quantMode GeneCounts \
            --outFileNamePrefix "$output_dir/$star_prefix" \
            --runThreadN 30 \
            --limitBAMsortRAM 200000000000 \
            --genomeLoad NoSharedMemory

        # Clean up temporary directory
        rm -r "$temp_dir"
    else
        echo "Warning: Missing files for $anumber. Skipping."
    fi
}

# Loop through the list of fastq files and process paired-end reads
while IFS= read -r line; do
    # Process only R1 files
    if [[ "$line" == *"_R1_"* ]]; then
        process_star_alignment "$line"
    fi
done < "$fastqz_list"
