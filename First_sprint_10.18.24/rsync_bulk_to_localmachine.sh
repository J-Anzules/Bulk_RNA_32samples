#!/bin/bash
#SBATCH --job-name=rsync_job          # Job name
#SBATCH --output=rsync_output.log      # Output log file
#SBATCH --error=rsync_error.log        # Error log file
#SBATCH --time=01:00:00               # Time limit (hh:mm:ss)
#SBATCH --mem=1G                      # Memory allocation
#SBATCH --cpus-per-task=1             # Number of CPU cores
#SBATCH --partition=your_partition     # Specify partition (if needed)

# Load any required modules (if needed, adjust accordingly)
# module load rsync

# Variables
SRC="/coh_labs/yunroseli/fastq_DoNotTouch/20241015_LH00295_0140_B22VL37LT3"          # Source directory (remote or another local path)
DEST="/mnt/c/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data"

# Execute rsync command
rsync -avzh $SRC $DEST

echo "rsync completed successfully"

