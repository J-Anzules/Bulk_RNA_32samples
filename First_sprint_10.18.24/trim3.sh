# Define your input and output directories
file_locations="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed"
cut_output="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/trimmed_no_empty_2"

# Loop through all the R1 files in the input directory
for r1_file in "$file_locations"/*_R1.fastq.gz; do
    # Define the corresponding R2 file
    r2_file="${r1_file/_R1/_R2}"
    
    # Check if R2 file exists
    if [ -f "$r2_file" ]; then
        # Get the basenames (just the file names without the directory)
        r1_base=$(basename "$r1_file")
        r2_base=$(basename "$r2_file")
        
        # Define the output file paths for R1 and R2
        output_r1="$cut_output/$r1_base"
        output_r2="$cut_output/$r2_base"
        
        # Run cutadapt for paired-end reads and filter sequences with 90 bp or longer
        cutadapt -m 90 --paired-output "$output_r2" -o "$output_r1" --cores=30 "$r1_file" "$r2_file"
        
        echo "Processed: $r1_file and $r2_file -> $output_r1 and $output_r2"
    else
        echo "Warning: Skipping $r1_file as the corresponding R2 file was not found."
    fi
done
