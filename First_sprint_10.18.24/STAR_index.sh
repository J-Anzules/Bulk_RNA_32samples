#!/bin/bash
#SBATCH --job-name=Indexing.sh    # Job name
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Index.log
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/Jon/32_samples_10182024/slurm_job_output/Index_error.log         # Standard error log
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --time=250:00:00                     # Set an appropriate time limit (10 hours in this case)
#SBATCH --mem=80G                          # Adjust memory as necessary
#SBATCH -n 1                # Single task
#SBATCH --cpus-per-task=32   # 30 CPUs for that one task
#at least 10-15 GB per core is recommended for large genomes

# Set directories
genome_dir="/labs/yunroseli_grp/Jon/Hg38/genome_assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_file="/labs/yunroseli_grp/Jon/Hg38/Gene_annotation-gtf/Homo_sapiens.GRCh38.113.gtf"
regulatory_gff="/labs/yunroseli_grp/Jon/Hg38/regulatory_features/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20240230.gff"
index_output_dir="/labs/yunroseli_grp/Jon/Hg38/STAR_index/"

# Create the output directory if it doesn't exist
mkdir -p "$index_output_dir"


module load STAR/2.7.9a
# Run STAR to generate the index
STAR --runMode genomeGenerate \
    --genomeDir "$index_output_dir" \
    --genomeFastaFiles "$genome_dir" \
    --sjdbGTFfile "$gtf_file" \
    --sjdbOverhang 149 \
    --runThreadN 32 \
    --genomeSAindexNbases 14 \
    --limitGenomeGenerateRAM 80000000000  \
    --sjdbFileChrStartEnd "$regulatory_gff"


