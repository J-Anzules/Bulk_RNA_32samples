#!/bin/bash
#SBATCH --job-name=STAR_test_A13487    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_test_A13487.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_test_A13487_error.log         # Standard error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=2:00:00                     # Set an appropriate time limit (10 hours in this case)
#SBATCH --mem=64G                          # Adjust memory as necessary
#SBATCH -n 1                # Single task
#SBATCH --cpus-per-task=30   # 30 CPUs for that one task

# Define input files (R1 and R2 for the same sample)
R1_strand="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimtest3/A13487_R1.fastq.gz"
R2_strand="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimtest3/A13487_R2.fastq.gz"

# STAR related variables
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"
temp_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp"
output_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results"
star_prefix="A13487"

# Load STAR module
# module load STAR/2.7.9a
module load  STAR/2.5.3a-foss-2017a
# Run STAR alignment
STAR --genomeDir $genome_directory \
    --readFilesCommand zcat \
    --readFilesIn $R1_strand $R2_strand \
    --outTmpDir $temp_dir \
    --outSAMtype BAM Unsorted \
    --quantMode GeneCounts \
    --outFileNamePrefix $output_dir/$star_prefix \
    --runThreadN 30 \
    --limitBAMsortRAM 40000000000

rm -r /labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp
