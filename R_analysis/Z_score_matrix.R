library(edgeR)
library(DESeq2)
library(pheatmap)

#-------------------Normalization edgeR-----------------------------------
CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv"
"C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/"
ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)
ct_mtx$ensembl_gene_id <- NULL
# norm.ct_mtx <- ct_mtx

## as DGEList
dge <- DGEList(counts=ct_mtx)

## calculate norm. factors
dge <- calcNormFactors(dge)

## get normalized counts
norm.ct_mtx <- cpm(dge)
#-------------------Normalization DESEQ2-----------------------------------

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


# Run the DESeq2 analysis
dds <- DESeq(dds)

# Step 3: Extract the normalized counts
norm.ct_mtx <- counts(dds, normalized=TRUE)

#-------------------------Changing Names----------------------------
# Load your CSV file
metadata <- read.csv("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/metadata_short.csv")

# create a mapping to the column names
id_mapping <- setNames(metadata$Sample_ID, metadata$TGen_Sample_ID)

# Rename columns of ct_mtx using this mapping
colnames(norm.ct_mtx) <- id_mapping[colnames(norm.ct_mtx)]

# Map the Condition (RT/NRT) to the Sample_ID (which are now the column names of norm.ct_mtx)
annotation_col <- data.frame(Condition = metadata$Condition)
rownames(annotation_col) <- metadata$Sample_ID


#-------------------Sequences of interest---------------------------

gene_lists <- list(
  repair = c(
    "ATP6V1C1", "ATP6V1G3", "CDH1", "CENPN", "CHMP2A", "CHMP4C", 
    "COX6A1", "COX6C", "COX7C", "COX8A", "EGF", "ELOB", 
    "EYA2", "GMNN", "GSR", "H1-0", "H1-2", "H2AC18", 
    "H2AC19", "H2AC8", "H2AJ", "H2BC10", "H2BC12", 
    "H2BC14", "H2BC21", "H2BC3", "H2BC6", "H2BC7", 
    "H2BC8", "H3C4", "H3C6", "H4C5", "H4C8","H4C8", "HMGA2", 
    "HSP90AA1", "HSP90AB1", "NCOR1", "PHC3", "PRDX5", 
    "RAB1B", "RAD21", "RAD51", "RPN1", "RPN2", 
    "RPS27A", "RPS3A", "SMC4", "SUMO3", "SYCP3", 
    "TBL1XR1", "TMED2", "TOPBP1", "TUBB4B", "UBA52", 
    "YWHAE"
  ),
  rt_resistance = c(
    "ATM", "ATR", "BTC", "CRIPTO", "DPPA4", "EOMES", 
    "EREG", "ESR1", "FABP7", "FGF1", "FGF10", 
    "FGF16", "FGF17", "FGF18", "FGF20", "FGF22", 
    "FGF23", "FGF3", "FGF4", "FGF5", "FGF6", 
    "FGF8", "FGF9", "FGFR2", "FLRT2", "FOXD3", 
    "FOXP3", "HES5", "HEY1", "HEY2", "HEYL", 
    "KEAP1", "KL", "LTA", "LTB", "MAMLD1", 
    "MT1A", "NANOG", "NRG1", "NRG2", "NRG3", 
    "NRG4", "POU5F1", "PRL", "PTCRA", "RSP03", 
    "SALL1", "SALL4", "SNAI1", "SOX2", "STAT3", 
    "TBXT", "TNFRSF11A", "TNFRSF12A", "TNFRSF13C", 
    "TNFSF11", "TNFSF12", "TNFSF14", "ZIC3"
  ),
  immune = c(
    "AP1M2", "APP", "ATP6V1G1", "BPIFB2", "CD177", 
    "CDC26", "CRISP2", "DSP", "ELOB", "GLA", 
    "HSP90AA1", "HSP90AB1", "HSPA8", "IDH1", "IL1R1", 
    "KIF5C", "LAMP2", "LCP1", "MALT1", "MIF", 
    "MX1", "NAT1", "P4HB", "PDIA3", "PPIA", 
    "RPN1", "RPN2", "RPS27A", "RPS6KA3", "S100A11", 
    "SAA1", "SEC61B", "SEC61G", "SFTPA2", "SRP14", 
    "TRIM68", "TSTD1", "TUBB4B", "UBA52", "VAMP8"
  )
  
)


# Checking if any genes are not found in the matrix
# repair[!repair %in% rownames(ct_mtx)]
# rt_resistance[!rt_resistance %in% rownames(ct_mtx)]
# immune[!immune %in% rownames(ct_mtx)]

output_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/zscores/"


#----------------------------Zscores----------------------------------------
# Update the heatmap generation function to include the sample annotations
generate_heatmap <- function(gene_list, gene_set_name, norm_ct_mtx, output_dir, annotation_col) {
  # Subset the normalized matrix to only include the genes of interest
  subset_mtx <- norm_ct_mtx[rownames(norm_ct_mtx) %in% gene_list, ]
  
  # Calculate Z-scores for each gene across all samples
  z_scores <- t(apply(subset_mtx, 1, function(x) (x - mean(x)) / sd(x)))
  
  # Remove rows with NA values
  z_scores <- z_scores[complete.cases(z_scores), ]
  
  # Generate the heatmap with sample annotation
  heatmap_file <- paste0(output_dir, gene_set_name, "_zscore_heatmap.png")
  pheatmap(z_scores,
           cluster_rows = FALSE,        # Cluster genes
           cluster_cols = FALSE,        # Cluster samples
           color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
           main = paste("Z-score Heatmap:", gene_set_name),
           annotation_col = annotation_col,  # Add annotation for RT/NRT
           annotation_colors = list(Condition = c("RT" = "red", "NRT" = "blue")),
           filename = heatmap_file,    # Save to file
           width = 10, height = 10, dpi = 300  # Adjust size and resolution
  )
  
  cat("Heatmap saved for", gene_set_name, "at:", heatmap_file, "\n")
}

# Loop over each gene set and generate the heatmaps with sample annotations
for (gene_set_name in names(gene_lists)) {
  generate_heatmap(gene_list = gene_lists[[gene_set_name]], 
                   gene_set_name = gene_set_name, 
                   norm_ct_mtx = norm.ct_mtx, 
                   output_dir = output_dir,
                   annotation_col = annotation_col)  # Pass annotation_col for RT/NRT labeling
}


#----------------------------Normlized, no Zscore-----------------------------

generate_heatmap_norm <- function(gene_list, gene_set_name, norm_ct_mtx, output_dir, annotation_col) {
  # Subset the normalized matrix to only include the genes of interest
  subset_mtx <- norm_ct_mtx[rownames(norm_ct_mtx) %in% gene_list, ]
  
  # Generate the heatmap with sample annotation
  heatmap_file <- paste0(output_dir, gene_set_name, "_norm_heatmap.png")
  pheatmap(subset_mtx,
           cluster_rows = FALSE,        # Cluster genes
           cluster_cols = FALSE,        # Cluster samples
           color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
           main = paste("Normalized Count Heatmap:", gene_set_name),
           annotation_col = annotation_col,  # Add annotation for RT/NRT
           annotation_colors = list(Condition = c("RT" = "red", "NRT" = "blue")),
           filename = heatmap_file,    # Save to file
           width = 10, height = 10, dpi = 300  # Adjust size and resolution
  )
  
  cat("Heatmap saved for", gene_set_name, "at:", heatmap_file, "\n")
}

# Loop over each gene set and generate the heatmaps with sample annotations
for (gene_set_name in names(gene_lists)) {
  generate_heatmap_norm(gene_list = gene_lists[[gene_set_name]], 
                   gene_set_name = gene_set_name, 
                   norm_ct_mtx = norm.ct_mtx, 
                   output_dir = output_dir,
                   annotation_col = annotation_col)  # Pass annotation_col for RT/NRT labeling
}

#---------------------------Log heatmap--------------------------------------
generate_heatmap_log <- function(gene_list, gene_set_name, norm_ct_mtx, output_dir, annotation_col) {
  # Subset the normalized matrix to only include the genes of interest
  subset_mtx <- norm_ct_mtx[rownames(norm_ct_mtx) %in% gene_list, ]
  
  # Apply log2 transformation with a pseudocount (to avoid log(0))
  log_subset_mtx <- log2(subset_mtx + 1)
  
  # Generate the heatmap with sample annotation
  heatmap_file <- paste0(output_dir, gene_set_name, "_log2_norm_heatmap.png")
  pheatmap(log_subset_mtx,
           cluster_rows = FALSE,        # Cluster genes
           cluster_cols = FALSE,        # Cluster samples
           color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
           main = paste("Log2 Normalized Count Heatmap:", gene_set_name),
           annotation_col = annotation_col,  # Add annotation for RT/NRT
           annotation_colors = list(Condition = c("RT" = "red", "NRT" = "blue")),
           filename = heatmap_file,    # Save to file
           width = 10, height = 10, dpi = 300  # Adjust size and resolution
  )
  
  cat("Heatmap saved for", gene_set_name, "at:", heatmap_file, "\n")
}

# Loop over each gene set and generate the heatmaps with sample annotations
for (gene_set_name in names(gene_lists)) {
  generate_heatmap_log(gene_list = gene_lists[[gene_set_name]], 
                        gene_set_name = gene_set_name, 
                        norm_ct_mtx = norm.ct_mtx, 
                        output_dir = output_dir,
                        annotation_col = annotation_col)  # Pass annotation_col for RT/NRT labeling
}

