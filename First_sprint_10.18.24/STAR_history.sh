#!/bin/bash
#SBATCH --job-name=STARtest.sh    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STARtest.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/STARtest_error.log         # Standard error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=250:00:00                     # Set an appropriate time limit (10 hours in this case)
#SBATCH --mem=32G                          # Adjust memory as necessary
#SBATCH -n 1                # Single task
#SBATCH --cpus-per-task=30   # 30 CPUs for that one task

# file_list="/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/trimmed/list_of_trimmed.txt"
file_list="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/code/list_of_trimmed.txt"
Original_file_list="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/20241015_LH00295_0140_B22VL37LT3/fastq_files.txt"
#TODO: rsyng this to the cluster
all_unique_identifiers="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/code/all_unique_Anumbers.txt"

# STAR related variables
genome_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index"
temp_dir="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp"

failed_identifier=""

# Store unique A numbers into an array
unique_A_numbers=$(sed -E 's|.*/(A[0-9]+)_R.*|\1|' "$file_list" | sort | uniq)

# module load  STAR/2.5.3a-foss-2017a
# Loop through the unique identifiers
for identifier in $unique_A_numbers; do
    echo "Starting processing for: $identifier"
    
    # Use grep to pull lines that contain the identifier
    matching_files=$(grep "$identifier" "$file_list")
    
    # Count the number of matching lines
    file_count=$(echo "$matching_files" | wc -l)
    
    # Check if exactly two lines are found
    if [ "$file_count" -eq 2 ]; then
        echo "Found two matching files for $identifier:"
        echo "$matching_files"

        # Split matching_files by lines and identify R1 and R2
        while IFS= read -r line; do
            if [[ "$line" == *"R1"* ]]; then
                R1_strand="$line"
            elif [[ "$line" == *"R2"* ]]; then
                R2_strand="$line"
            fi
        done <<< "$matching_files"

        echo "R1 strand: $R1_strand"
        echo "R2 strand: $R2_strand"
        star_prefix=$(echo "$identifier")

        # STAR --genomeDir $genome_directory \
        # --readFilesCommand zcat \
        # --readFilesIn $R1_strand $R2_strand \
        # --outTmpDir $temp_dir \
        # --outSAMtype BAM Unsorted \
        # --quantMode GeneCounts \
        # --outFileNamePrefix /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results/$star_prefix \
        # --runThreadN 30 \
        # --limitBAMsortRAM 40000000000




        # Here you can process the two matching files as needed
    else
        echo "Failed to find exactly two matching files for $identifier"
        # Append the identifier to the failed list
        failed_identifiers+="$identifier "
    fi

done

# Output the failed identifiers at the end
if [ -n "$failed_identifiers" ]; then
    echo "###########################################################"
    echo "###########################################################"
    echo " "
    echo " "
    echo "Identifiers with no matching files or incorrect file count:"
    echo "$failed_identifiers"
    echo " "
    echo " "
     # Save the failed identifiers to a file
    echo "$failed_identifiers" > ./failed_alignment.txt
    echo "Failed identifiers have been saved to ./failed_alignment.txt"
    echo " "
    echo " "
    echo "###########################################################"
    echo "###########################################################"



fi




    echo "###########################################################"
    echo "###########################################################"
    echo " "
    echo " "
    echo "Identifiers with no matching files or incorrect file count:"
    echo "$failed_identifiers"
    echo " "
    echo " "
    echo "###########################################################"
    echo "###########################################################"







#####################################################################
#####################################################################
#####################################################################
# # Define the original file path (input fastq)
# original_file="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/YRLP_0001/YRLP_0001_1_UN_Whole_C1_UUUUU_A13487_22VL37LT3_CGTCCATGTA_L007_R1_001.fastq.gz"

# # Assign R1 and R2 strand files
# R1_strand="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/YRLP_0001/YRLP_0001_1_UN_Whole_C1_UUUUU_A13487_22VL37LT3_CGTCCATGTA_L007_R1_001.fastq.gz"
# R2_strand="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/YRLP_0001/YRLP_0001_1_UN_Whole_C1_UUUUU_A13487_22VL37LT3_CGTCCATGTA_L007_R2_001.fastq.gz"

# # Generate new name based on the original file name
# new_name=$(echo "$original_file" | sed -E 's/.*_(A[0-9]+)_.*(R[12]).*/\1_\2.fastq.gz/')

# # Extract the A number for the STAR prefix (correcting the sed command)
# star_prefix=$(echo "$original_file" | sed -E 's/.*_(A[0-9]+).*/\1/')

# # Output the new names
# echo "New Name: $new_name"
# echo "STAR Prefix: $star_prefix"


# STAR --genomeDir $genome_directory \
# --readFilesCommand zcat \
# --readFilesIn $R1_strand $R2_strand \
# --outTmpDir $temp_dir \
# --outSAMtype BAM Unsorted \
# --quantMode GeneCounts \
# --outFileNamePrefix /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results/$star_prefix \
# --runThreadN 30 \
# --limitBAMsortRAM 40000000000


# module load  STAR/2.5.3a-foss-2017a
# STAR --genomeDir /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index \
# --readFilesCommand zcat \
# --readFilesIn $R1_strand $R2_strand \
# --outTmpDir /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/tmp \
# --outSAMtype BAM Unsorted \
# --quantMode GeneCounts \
# --outFileNamePrefix /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/data/20241015_LH00295_0140_B22VL37LT3/STAR_output/STAR_results/$star_prefix \
# --runThreadN 30 \
# --limitBAMsortRAM 40000000000