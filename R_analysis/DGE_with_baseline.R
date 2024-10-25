library(DESeq2)
library(ggplot2)
library(pheatmap)

#--------------------------------Preparing Data--------------
CountMatrix_loc <- "~/Documents/RoseLab/Helya/counts_without_duplicates.txt"
ct_mtx <- read.delim(CountMatrix_loc, header = TRUE, row.names = 1, sep = "\t")
ct_mtx$ID.1 <- NULL

# Getting rid of some float numbers
ct_mtx <- round(ct_mtx)

# [1] "P_13D_Non_RT_0hr"  "P_13D_RT_0hr"      "P_13D_Non_RT_4hr"  "P_13D_RT_4hr"     
# [5] "P_13D_Non_RT_24hr" "P_13D_RT_24hr"     "P_6B_Non_RT_0hr"   "P_6B_RT_0hr"      
# [9] "P_6B_RT_0.5hr"     "P_6B_RT_2hr"       "P_6D_Non_RT_0hr"   "P_6D_RT_0hr"      
# [13] "P_6D_Non_RT_0.5hr" "P_6D_RT_0.5hr"     "P_6D_Non_RT_2hr"   "P_6D_RT_2hr"


# Create a metadata table with "baseline" (all 0hr samples) vs "RT_irradiated" (all irradiated samples)
sample_metadata_new <- data.frame(
  row.names = colnames(ct_mtx),
  condition = c("baseline", "baseline", "Non_RT", "RT_irradiated", 
                "Non_RT", "RT_irradiated", "baseline", "baseline", 
                "RT_irradiated", "RT_irradiated", "baseline", "baseline", 
                "RT_irradiated", "RT_irradiated", "RT_irradiated", "RT_irradiated"),
  timepoint = c("0hr", "0hr", "4hr", "4hr", "24hr", "24hr", "0hr", "0hr", 
                "0.5hr", "2hr", "0hr", "0hr", "0.5hr", "0.5hr", "2hr", "2hr")
)

####################################################################################
#           THIS IS THE AREA THAT SHOULD BE EDITED FOR VARIATION IN ANALYSIS
# Subset the count matrix for baseline and irradiated samples
# count_matrix_subset <- ct_mtx[, c("P_13D_Non_RT_0hr", "P_13D_RT_0hr",
#                                   "P_6B_Non_RT_0hr", "P_6B_RT_0hr", 
#                                   "P_6D_Non_RT_0hr", "P_6D_RT_0hr", 
#                                   "P_13D_RT_4hr", "P_6B_RT_0.5hr", 
#                                   "P_6B_RT_2hr", "P_6D_RT_0.5hr", 
#                                   "P_6D_RT_2hr")]

# ONLY KEEPING HOUR 4 HERE
count_matrix_subset <- ct_mtx[, c("P_13D_Non_RT_0hr", "P_13D_RT_0hr", 
                                  "P_6B_Non_RT_0hr", "P_6B_RT_0hr", 
                                  "P_6D_Non_RT_0hr", "P_6D_RT_0hr", 
                                  "P_13D_RT_4hr")]
####################################################################################
# Subset metadata to match the columns
metadata_subset <- sample_metadata_new[colnames(count_matrix_subset), ]

#-------------------------------Deseq df--------------------------
# Create DESeq2 dataset comparing baseline to irradiated samples
dds <- DESeqDataSetFromMatrix(countData = count_matrix_subset, 
                              colData = metadata_subset, 
                              design = ~ condition)

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Get the results comparing irradiated (RT_irradiated) samples to the baseline
res <- results(dds, contrast = c("condition", "RT_irradiated", "baseline"))

# Summary of results, showing the number of significant genes
summary(res)

# Count the number of genes with adjusted p-value < 0.05
significant_genes <- res[which(res$pvalue < 0.05), ]
nrow(significant_genes)

