library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(dplyr)

#--------------------------------Preparing Data--------------
CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv"

ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)
ct_mtx$ensembl_gene_id <- NULL

# Load your CSV file
metadata <- read.csv("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/metadata_short.csv")

# Subset the metadata based on the samples in the count matrix
metadata_deseq2 <- metadata %>%
  filter(TGen_Sample_ID %in% colnames(ct_mtx)) %>%
  select(TGen_Sample_ID, Condition, Hour)


# Reorder metadata to match the column order of ct_mtx
metadata_deseq2 <- metadata_deseq2 %>%
  filter(TGen_Sample_ID %in% colnames(ct_mtx)) %>%
  arrange(match(TGen_Sample_ID, colnames(ct_mtx)))

# Check if the metadata order matches the column names of the count matrix
all(metadata_deseq2$TGen_Sample_ID == colnames(ct_mtx))


# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = ct_mtx,                  # Your count matrix
  colData = metadata_deseq2,           # Your metadata with sample information
  design = ~ Condition                 # Design formula specifying your condition
)

# # Create the DESeq2 dataset
# dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
#                               colData = metadata, 
#                               design = ~ condition)

# head(count_matrix)
# head(metadata)
#-------------------------------DGE Analysis ----------------------
# Run the DESeq2 analysis
dds <- DESeq(dds)

# Get results for RT vs. Non-RT (default comparison: RT vs Non-RT)
res <- results(dds, contrast = c("Condition", "RT", "NRT"))

# # TODO: spot for automation
# write.csv(res, "~/Documents/RoseLab/BulkRNA_seq/data/DGE_results/first_16_samples/hour_0_vs_2.csv",
#           row.names = TRUE)
# View the top differentially expressed genes
# head(res)

# Summary of results, showing the number of significant genes
# summary(res)

# Count the number of genes with adjusted p-value < 0.05
significant_genes_padj <- res[which(res$padj < 0.05), ]
significant_genes_pval <- res[which(res$pvalue < 0.05), ]

write.csv(res, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/results_dge.csv")
write.csv(res, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/results_dge.csv")
write.csv(significant_genes_pval, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/significant_genes.csv")
head(significant_genes_padj)
head(significant_genes_pval)

nrow(significant_genes_padj)
nrow(significant_genes_pval)
# nrow(significant_genes)
# 
# summary(significant_genes$log2FoldChange)
#----------------------------------Visualization--------------------------

##-----------------------------Volcano Plot--------------------------
#Subset top up/down regulated genes
Upregulated_Genes <- subset(res, log2FoldChange > 0 )
Downregualted_Genes <- subset(res, log2FoldChange < 0)

# Choosing the number of top/down
Number_of_TopnDown = 6


# Select the top 10 upregulated and top 10 downregulated genes
top_upregulated <- Upregulated_Genes[order(Upregulated_Genes$pvalue, na.last = NA)[1:Number_of_TopnDown],]
top_downregulated <- Downregualted_Genes[order(Downregualted_Genes$pvalue, na.last = NA)[1:Number_of_TopnDown],]

# Combine the two dataframes
selection <- rbind(top_upregulated, top_downregulated)

# Create a volcano plot
volcano_plot_fig <- ggplot(as.data.frame(res) , aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = "grey", alpha = 0.6) +
  geom_point(data = subset(res, pvalue < 0.05 & abs(log2FoldChange) > 1),
             color = "red", alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  # Adding the text labels for the top genes
  geom_text(data = selection, aes(label = rownames(selection)), 
            vjust = -1, hjust = 1) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot") +
  theme_bw()
# Volcano Plot
ggsave(filename = "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/Volcano_plt_6up_n_down.png", 
       plot = volcano_plot_fig,   # The name of your ggplot object
       width = 10,       # Width of the saved image (in inches)
       height = 6,       # Height of the saved image (in inches)
       dpi = 300)        # Resolution (300 dpi is high quality)

##------------------------------Heat Map--------------------



# ~~~heatmat without a selection of a more balanced display was not imformative.
# ~~~There was a lot more downregulated genes than upregulated

# Extract normalized counts for the significant genes
normalized_counts <- counts(dds, normalized = TRUE)
sig_gene_counts <- normalized_counts[rownames(significant_genes), ]
# 
# # Create a heatmap of significant genes
# pheatmap(sig_gene_counts, cluster_rows = TRUE, cluster_cols = TRUE,
#          main = "Heatmap of Significant Genes (RT vs Non-RT at 0hr)")

# Sort significant genes by log2FoldChange to get top 10 upregulated and top 10 downregulated genes
top_10_up_hm <- significant_genes[order(-significant_genes$log2FoldChange), ][1:15, ]
top_10_down_hm <- significant_genes[order(significant_genes$log2FoldChange), ][1:15, ]

# Combine the two sets of top genes
top_20_genes_hm <- rbind(top_10_up_hm, top_10_down_hm)

# Get the gene names (rownames) of these top 20 genes
top_20_gene_names_hm <- rownames(top_20_genes_hm)

# Extract normalized counts from the DESeq2 object
normalized_counts <- counts(dds, normalized = TRUE)

# Subset the normalized counts for the top 20 upregulated and downregulated genes
top_20_counts_hm <- normalized_counts[top_20_gene_names_hm, ]

# Generate the heatmap
pheatmap(top_20_counts_hm, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Heatmap of Top 10 Upregulated and Downregulated Genes",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         filename = "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/heatmap_topUpnDown.png" #TODO: Spot for automation
         )  # Adjust color scale

##-------------------------------Saving figure -------------------------------------


# Load necessary libraries
library(DESeq2)
library(pheatmap)

# Step 1: Extract normalized counts for the significant genes
normalized_counts <- counts(dds, normalized = TRUE)

# Step 2: Sort significant genes by log2FoldChange
# Top 15 upregulated genes (highest positive log2FoldChange)
top_10_up_hm <- significant_genes[order(-significant_genes$log2FoldChange), ][1:15, ]

# Top 15 downregulated genes (lowest negative log2FoldChange)
top_10_down_hm <- significant_genes[order(significant_genes$log2FoldChange), ][1:15, ]

# Combine the two sets of top genes
top_20_genes_hm <- rbind(top_10_up_hm, top_10_down_hm)

# Get the gene names (rownames) of these top 20 genes
top_20_gene_names_hm <- rownames(top_20_genes_hm)

# Step 3: Subset the normalized counts for the top 20 genes
top_20_counts_hm <- normalized_counts[top_20_gene_names_hm, ]


# Step 4: Create a sample annotation table for RT vs NRT
# Extract the condition information from metadata_deseq2
groups <- metadata_deseq2$Condition
names(groups) <- metadata_deseq2$TGen_Sample_ID

# Create a group annotation for the heatmap
grpann <- data.frame(Group = groups)

# Map "NRT" to "Non-RT" and "RT" to "RT" for display in the heatmap
grpann$DisplayGroup <- ifelse(grpann$Group == "NRT", "Non-RT", "RT")

# Assign annotation colors for the heatmap (blue for "Non-RT", red for "RT")
annCol <- list(DisplayGroup = c("Non-RT" = "blue", "RT" = "red"))

# Ensure row names of the annotation match the column names of the count matrix
rownames(grpann) <- metadata_deseq2$TGen_Sample_ID

# Step 5: Remove control and diabetic suffixes in column names of the count matrix
colnames(top_20_counts_hm) <- gsub("_[CD]$", "", colnames(top_20_counts_hm))


# Step 5: Apply log2 transformation to the normalized counts (adding a pseudocount to avoid log(0))
log2_normalized_counts <- log2(top_20_counts_hm + 1)  # Add 1 to avoid log(0)


# Step 1: Create a named vector mapping TGen_Sample_ID to Sample_ID
id_mapping <- setNames(metadata$Sample_ID, metadata$TGen_Sample_ID)

# Step 2: Rename the columns in log2_normalized_counts based on the mapping
colnames(log2_normalized_counts) <- id_mapping[colnames(log2_normalized_counts)]

# Step 3: Verify that the column names have been updated
print(colnames(log2_normalized_counts))

# Generate the heatmap with log2-transformed values
pheatmap(log2_normalized_counts, 
         cluster_rows = TRUE,        # Cluster genes
         cluster_cols = TRUE,        # Cluster samples
         main = "Heatmap of Top 10 Upregulated and Downregulated Genes (Log2 Transformed)",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         scale = "row",              # Scale across rows (genes)
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Adjust color scale
         annotation_col = grpann["DisplayGroup"],  # Annotation for the columns
         annotation_colors = annCol,  # Color scheme for annotations
         filename = "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/heatmap_topUpnDown_log2.png", 
         width = 10, height = 8, dpi = 600  # Adjust size and resolution
)



