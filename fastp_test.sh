#!/bin/bash
#SBATCH --job-name=fastp_test    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/fastp_test.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/fastp_test_error.log
#SBATCH --ntasks=1              # Run a single task
#SBATCH --cpus-per-task=25      # Request all CPUs
#SBATCH --mem=100G              # Request 200 GB of memory
#SBATCH --time=24:00:00         # Set an appropriate time limit

# Define file locations
fastqz_list="/labs/yunroseli_grp/Jon/32_samples_10182024/code/mappings/test_fastq.txt"
fastp_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/fastp_test"

# STAR locations
star_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/STAR_test/"
star_temp_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/STAR_test/tmp_files"
output_dir="/labs/yunroseli_grp/Jon/32_samples_10182024/data/STAR_test/fastp_test"
# genome_directory="/labs/yunroseli_grp/Jon/Hg38/STAR_index"
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"


# In case I have do rerun a test
rm -rf "$fastp_output"
mkdir -p "$fastp_output"
mkdir -p "$star_temp_dir"
mkdir -p "$output_dir"
# Load the necessary modules
module load fastp
# module load STAR/2.7.9a
module load  STAR/2.5.3a-foss-2017a

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

        if [[ -f "$output_r1" && -f "$output_r2" ]]; then
             
             # Create unique temp directory for each sample
            tmp_dir="$star_temp_dir/tmp_${anumber}"
            echo "Processing: $anumber"
            star_prefix="${anumber}_"


            # Run STAR alignment
            STAR --genomeDir $genome_directory \
                --readFilesCommand zcat \
                --readFilesIn "$output_r1" "$output_r2" \
                --outTmpDir "$tmp_dir" \
                --outSAMtype BAM Unsorted \
                --quantMode GeneCounts \
                --outFileNamePrefix "$output_dir/$star_prefix" \
                --runThreadN 25 \
                --limitBAMsortRAM 100000000000 \
                --genomeLoad NoSharedMemory

            # Clean up temporary directory
            rm -r "$tmp_dir"
        else
            echo "fastp must have failed for either $anumber"

        fi
    fi
done < "$fastqz_list"
