#!/bin/bash
#SBATCH --job-name=5STAR    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/5STAR_loop.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/5STAR_loop_error.log  # Error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=24:00:00               # Set an appropriate time limit
#SBATCH --mem=150G                    # Adjust memory as necessary
#SBATCH --cpus-per-task=20            # 30 CPUs

# Define file locations
file_locations="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trim5"
output_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_trim5"
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"
anumber_file="/labs/yunroseli_grp/Jon/32_samples_10182024/code/mapping_files/all_unique_Anumbers.txt"

# Load STAR module
module load  STAR/2.5.3a-foss-2017a

# Loop through the list of unique A numbers
while read -r anumber; do
    # Define input files (R1 and R2 for each sample)
    R1_strand="${file_locations}/${anumber}_R1.fastq.gz"
    R2_strand="${file_locations}/${anumber}_R2.fastq.gz"
    
    # Create unique temp directory for each sample
    temp_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp_${anumber}"
    
    # Check if both R1 and R2 exist for the sample
    if [ -f "$R1_strand" ] && [ -f "$R2_strand" ]; then
        echo "Processing: $anumber"
        star_prefix="$anumber"

        # Run STAR alignment
        STAR --genomeDir $genome_directory \
            --readFilesCommand zcat \
            --readFilesIn "$R1_strand" "$R2_strand" \
            --outTmpDir "$temp_dir" \
            --outSAMtype BAM Unsorted \
            --quantMode GeneCounts \
            --outFileNamePrefix "$output_dir/$star_prefix" \
            --runThreadN 10 \
            --limitBAMsortRAM 100000000000 \
            --genomeLoad NoSharedMemory

        # Clean up temporary directory
        rm -r "$temp_dir"
    else
        echo "Warning: Missing files for $anumber. Skipping."
    fi
done < "$anumber_file"
