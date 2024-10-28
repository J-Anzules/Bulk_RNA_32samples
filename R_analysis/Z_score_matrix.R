library(edgeR)
library(DESeq2)
library(pheatmap)
library(dplyr)

#-------------------Normalization edgeR-----------------------------------

CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_clean_fpStr.csv"
# CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv"

ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)
# ct_mtx$ensembl_gene_id <- NULL
# norm.ct_mtx <- ct_mtx

## as DGEList
dge <- DGEList(counts=ct_mtx)

## calculate norm. factors
dge <- calcNormFactors(dge)

## get normalized counts
norm.ct_mtx <- cpm(dge)
#-------------------Normalization DESEQ2-----------------------------------

CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_clean_fpStr.csv"
# CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv"

ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)
# ct_mtx$ensembl_gene_id <- NULL

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

##########################################################################
# TWO DIFFERENT TYPES OF NORMALIZATIONS

# Extract the normalized counts
norm.ct_mtx <- counts(dds, normalized=TRUE)

##########################################################################
# Variance-Stabilizing Transformation - Considering condition
vsd <- vst(dds, blind=FALSE)
norm.ct_mtx <- assay(vsd) # This is used for visualization

# Variance-Stabilizing Transformation - unbiased
vsd <- vst(dds, blind=TRUE)
norm.ct_mtx <- assay(vsd) # This is used for visualization


##########################################################################
# Perform the rlog transformation - Considering conditions
rld <- rlog(dds, blind = FALSE)
norm.ct_mtx <- assay(rld)

# Perform the rlog transformation - unbiased
rld <- rlog(dds, blind = TRUE)
norm.ct_mtx <- assay(rld)



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
  ),
  survival_stress_repair = c(
    "AREG", "ARG1", "ATR", "BAX", "BCL2", 
    "BCL2L1", "BCL2L15", "BRCA1", "BRCA2", 
    "CAT", "CCNG1", "CDK1", "CDK4", "CDK6", 
    "CDKN1A", "CGAS", "DDB2", "EGR1", "ERCC1", 
    "ERCC4", "ERCC5", "EZH2", "FANCG", "FDXR", 
    "FOS", "FOSB", "FOSL1", "GADD45A", "GART", "GSK3B",
    "H2AX", "HDAC1", "HIF1A", "HOXB9", "HSPA8", 
    "JUN", "KEAP1", "KLF4", "KLK3", "LIG4",
    "MDC1", "MDM2", "MKI67", "MRE11", "MTOR", "MYC", 
    "NANOG", "NFKB1", "NFKBIE", "PALB2",
    "PARP1", "PFKP", "PHPT1", "PRKDC", "PTTG1", "RAD50", 
    "RAD51", "RAD9A", "SESN1",
    "SMUG1", "SOX2", "STAT3", "STING1", "TERT", "TOP2A", 
    "TP53AIP1", "TP53INP2", "TP63",
    "TP73", "TXNRD1", "WEE1", "XRCC1", "XRCC5"
  ),
  metabolism = c(
    "ABCC5", "ACACA", "ACO1", "ACSL1", "ADHFE1",
    "ALDH6A1", "AMACR", "AMELX", "ARG2", "ATXN3",
    "CERK", "CMBL", "CNTN3", "COX6C", "COX7C",
    "CYB5A", "CYP4F2", "DBI", "DCUN1D1", "DCXR",
    "DECR1", "DEGS1", "DHCR7", "EDEM3", "EEF1A2",
    "FABP12", "FABP4", "FDFT1", "FOXO1", "FOXO3",
    "FOXO4", "GLYATL1", "GNG4", "GNMT", "GSK3A",
    "GSK3B", "H2BC14", "H4C5", "H4C8", "HACD2",
    "HMGA1", "HNF1A", "HNF1B", "HNF4A", "HNF4G",
    "HSP90AB1", "IDH1", "IGFBP1", "IGFBP3", "INPP4B",
    "INSR", "IRS1", "IRS2", "IRS4", "KCNC2",
    "MAOA", "MAPK1", "MRPS25", "N6AMT1", "NCOR1",
    "NUDT10", "ODC1", "PDX1", "PHC3", "PHGDH",
    "PIK3C2G", "PIK3R4", "PIP", "PLA2G4F", "PMP2",
    "PODXL2", "PTGR1", "PYCR1", "RAD21", "RBP4",
    "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14",
    "RPL37A", "RPL7A", "RPS3A", "RPS4X", "RPS5",
    "RPS6", "RPS7", "RPS8", "SC5D", "SCD",
    "SPTSSA", "TBL1XR1", "THRB", "TPTE", "UBP1",
    "UQCR10", "UQCRQ"
   ),
  apoptosis = c(
    "APAF1", "BAD", "BAK1", "BAX", "BBC3",
    "BCL2", "BID", "BLK", "CASP1", "CASP3",
    "CASP8", "CASP9", "CDH1", "CYCS", "DSG2",
    "DSP", "FAS", "GAPDH", "H1-0", "H1-2",
    "HTRA1", "HTRA2", "HTRA3", "HTRA4", "MAGED1",
    "MCL1", "NOXA1", "PMAIP1", "PRG3", "PXMP4",
    "RPN1", "RPN2", "RPS27A", "SMAD1", "SMAD2",
    "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7",
    "SMAD9", "STEAP3", "TNF", "TNFRSF10C", "TNFRSF10D",
    "TP53", "TP53AIP1", "TP63", "TP73", "TRADD",
    "UBA52", "YWHAE"
   ),
  house_keeping = c(
    "ACTB", "B2M", "GAPDH", "GUSB", "HMBS",
    "HPRT1", "PGK1", "PPIA", "RPL13A", "RPLP0",
    "SDHA", "TBP", "TFRC", "YWHAZ")
)



# Checking if any genes are not found in the matrix
# repair[!repair %in% rownames(ct_mtx)]
# rt_resistance[!rt_resistance %in% rownames(ct_mtx)]
# immune[!immune %in% rownames(ct_mtx)]

output_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/fasp_star_clean/heatmaps/"


#----------------------------Zscores----------------------------------------
# Update the heatmap generation function to include the sample annotations
generate_heatmap <- function(gene_list, gene_set_name, norm_ct_mtx, output_dir, annotation_col) {
  
  # Identify genes that are missing in the matrix
  missing_genes <- gene_list[!gene_list %in% rownames(norm_ct_mtx)]
  
  # Print missing genes for tracking
  if (length(missing_genes) > 0) {
    cat("Missing genes for", gene_set_name, ":", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Subset the normalized matrix to only include the genes of interest
  subset_mtx <- norm_ct_mtx[rownames(norm_ct_mtx) %in% gene_list, ]
  
  # Calculate Z-scores for each gene across all samples
  z_scores <- t(apply(subset_mtx, 1, function(x) (x - mean(x)) / sd(x)))
  
  # Remove rows with NA values
  z_scores <- z_scores[complete.cases(z_scores), ]
  
  # Generate the heatmap with sample annotation
  heatmap_file <- paste0(output_dir, gene_set_name, "_zscore_heatmap.png")
  pheatmap(z_scores,
           cluster_rows = TRUE,        # Cluster genes
           cluster_cols = TRUE,        # Cluster samples
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

