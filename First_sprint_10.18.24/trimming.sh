#!/bin/bash
#SBATCH --job-name=trim_n_QC    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/trim_n_QC.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/trim_n_QC_error.log         # Standard error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=250:00:00                     # Set an appropriate time limit (10 hours in this case)
#SBATCH --mem=32G                          # Adjust memory as necessary
#SBATCH -n 1                # Single task
#SBATCH --cpus-per-task=30   # 30 CPUs for that one task


# Define the path to the list of fastq.gz file addresses
fastqz_list="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/fastq_addresses.txt"

module load cutadapt

# Loop through each line (fastq.gz file) in the file
while IFS= read -r original_file; do
    echo "Processing: $original_file"
    
    # Create a new name based on the original file name
    new_name=$(echo "$original_file" | sed -E 's/.*_(A[0-9]+)_.*(R[12]).*/\1_\2.fastq.gz/')
    cutadapt_output="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed3/$new_name"
    
    # Run cutadapt to trim the file
    cutadapt -u 15 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a TGGAATTCTCGGGTGCCAAGG \
        -a "G{50}" \
        -o $cutadapt_output $original_file
    
    # Check if cutadapt ran successfully
    if [ $? -ne 0 ]; then
        echo "Error running cutadapt on $original_file"
        continue
    fi
    
    
    # Check if FastQC ran successfully
    if [ $? -ne 0 ]; then
        echo "Error running FastQC on $cutadapt_output"
        continue
    fi

    echo "Processing of $original_file completed."
    
done < "$fastqz_list"