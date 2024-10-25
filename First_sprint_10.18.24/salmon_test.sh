#!/bin/bash
#SBATCH --job-name=Salmon_test    # Job name
#SBATCH --output=/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Salmon_test.log  # Output log
#SBATCH --error=/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Salmon_test_error.log  # Error log
#SBATCH -N 1                      # Number of nodes
#SBATCH --cpus-per-task=30         # Number of CPU cores
#SBATCH --time=48:00:00            # Time limit
#SBATCH --mem=45G                  # Memory allocation

# Load the Salmon module (adjust the module name if necessary)
module load Salmon/1.4.0

# Define input FASTQ file and related directories
R1_strand="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed_removedempty_nless45bp/A13503_R1.fastq.gz"
R2_strand="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed_removedempty_nless45bp/A13503_R2.fastq.gz"
transcriptome_index="/path/to/salmon_transcriptome_index"  # Path to the Salmon index
output_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/Salmon_output"
sample_prefix="A13489"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Run Salmon quantification
salmon quant -i $transcriptome_index \
-l A \
-1 $R1_strand \
-2 $R2_strand \
-p 30 \
--validateMappings \
-o $output_dir/$sample_prefix

