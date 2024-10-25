#!/bin/bash
#SBATCH --job-name=STAR_test    # Job name
#SBATCH --output=/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_test.log  # Output log
#SBATCH --error=/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STAR_test_error.log  # Error log
#SBATCH -N 1                    # Number of nodes
#SBATCH --cpus-per-task=30       # Number of CPU cores
#SBATCH --time=48:00:00          # Time limit
#SBATCH --mem=45G                # Memory allocation

# Load the STAR module
# module load STAR/2.5.3a-foss-2017a
module load Salmon/1.4.0

# Define input FASTQ file and related directories
R2_strand="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed_removedempty_nless45bp/A13489_R1.fastq.gz"
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"
temp_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp"
output_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results"
star_prefix="A13491_R1"

# Run STAR alignment
STAR --genomeDir $genome_directory \
--readFilesCommand zcat \
--readFilesIn $R2_strand \
--outTmpDir $temp_dir \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts \
--outFileNamePrefix $output_dir/$star_prefix \
--runThreadN 30 \
--limitBAMsortRAM 40000000000
