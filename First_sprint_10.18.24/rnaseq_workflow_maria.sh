#Header for all SLURM jobs
#!/bin/bash
#SBATCH --job-name=STAR.sh    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mchosco@coh.org     # Where to send mail
#SBATCH -n 8                           # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=32G                      # Amount of memory in GB
#SBATCH --time=250:00:00               # Time limit hrs:min:sec
#SBATCH --output=STAR%j.log   # Standard output and error log

#Quality control FASTQC
fastqc -t 2 -o /path/to/output/folder *.fastq.gz

#Trim adapters TRIMGALORE
module load trimgalore/0.6.5
trim_galore *samplename.fastq.gz -q 30 -o /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/trim_results

#Align mapped reads and count reads per gene STAR
module load  STAR/2.5.3a-foss-2017a
STAR --genomeDir /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/genome/hg38_index --readFilesCommand zcat \
--readFilesIn /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/trim_results/sample#_R1_001_trimmed.fq.gz,sample#_R2_001_trimmed.fq.gz \
--outTmpDir /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/STAR_results/tmp \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts \
--outFileNamePrefix /net/nfs-irwrsrchnas01/labs/yunroseli_grp/Bladder_Cancer/Utah_Samples_Biospyder/Batch_1/STAR_results/sample#_ \
--runThreadN 30 \
--limitBAMsortRAM 8

#Count ouput format
#column 1: gene ID
#column 2: counts for unstranded RNA-seq
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

#Sort read count files and keep only first column of counts
#!/bin/bash

#List of input files
input_files=(*ReadsPerGene.out.tab)

#Loop through the list of files and run the sort command
for input_file in "${input_files[@]}"
do
  output_file=$(basename "${input_file}" .out.tab)_sorted.counts
  sort "${input_file}" | cut -f1,2 > "${output_file}"
  done


#Generate a counts matrix using using R script FYI you might only need to do this once? My last dataset was in 4 batches so I had to combine them all using this
#batch 1
setwd("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b1/")
files1  = list.files("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b1/", pattern = "*.txt")
files1
for(i in 1:length(files1)){

data1 = read.delim(files1[i], header=  F)
  if(i == 1){

    batch1 = as.data.frame(data1[,2])
  } else {

    batch1 = cbind.data.frame(batch1, data1[,2])
  }
}
colnames(batch1) = gsub("_sorted.txt", "", files1)
row.names(batch1) = data1$V1

#batch 2
setwd("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b2/")
files2  = list.files("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b2/", pattern = "*.txt")
files2
for(i in 1:length(files2)){

  data2 = read.delim(files2[i], header=  F)
  if(i == 1){

    batch2 = as.data.frame(data2[,2])
  } else {

    batch2 = cbind.data.frame(batch2, data2[,2])
  }
}
colnames(batch2) = gsub("_sorted.txt", "", files2)
row.names(batch2) = data2$V1

#batch 3
setwd("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b3/")
files3  = list.files("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b3/", pattern = "*.txt")
files3
for(i in 1:length(files3)){

  data3 = read.delim(files3[i], header=  F)
  if(i == 1){
    batch3 = as.data.frame(data3[,2])
  } else {

    batch3 = cbind.data.frame(batch3, data3[,2])
  }
}
colnames(batch3) = gsub("_sorted.txt", "", files3)
row.names(batch3) = data3$V1

#batch 4
setwd("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b4/")
files4  = list.files("/Users/mchosco/Desktop/Bladder_Cancer/count_data/sortedcounts_b4/", pattern = "*.txt")
files4
for(i in 1:length(files4)){

  data4 = read.delim(files4[i], header=  F)
  if(i == 1){

    batch4 = as.data.frame(data4[,2])
  } else {

    batch4 = cbind.data.frame(batch4, data4[,2])
  }
}
colnames(batch4) = gsub("_sorted.txt", "", files4)
row.names(batch4) = data4$V1

counts <- cbind(batch1,batch2, batch3, batch4)

write.csv(counts, file = "/Users/mchosco/Desktop/Bladder_Cancer/count_data/counts.csv")