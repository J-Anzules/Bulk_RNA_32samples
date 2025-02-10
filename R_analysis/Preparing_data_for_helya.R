library(DESeq2)
library(pheatmap)
library(dplyr)

#This is Kabyas 16 bulk samples
k_ctmtx <- read.table("C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/counts_without_duplicates.txt", 
                   header = TRUE, sep = "\t", row.names = 1)
k_ctmtx$ID.1 <- NULL

head(k_ctmtx)
# Getting rid of some float numbers
k_ctmtx <- round(k_ctmtx)

# Create a metadata table that includes sample condition (RT vs Non-RT)
sample_metadata <- data.frame(
  row.names = colnames(k_ctmtx),
  condition = c("Non_RT", "RT", "Non_RT", "RT", "Non_RT", "RT", "Non_RT", "RT",
                "RT", "RT", "Non_RT", "RT", "Non_RT", "RT", "Non_RT", "RT"),
  timepoint = c("0hr", "0hr", "4hr", "4hr", "24hr", "24hr", "0hr", "0hr",
                "0.5hr", "2hr", "0hr", "0hr", "0.5hr", "0.5hr", "2hr", "2hr")
)

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = k_ctmtx,
                              colData = sample_metadata,
                              design = ~ condition)


# Run the DESeq2 analysis
dds <- DESeq(dds)

# Extract the normalized counts
norm.k_ctmtx <- counts(dds, normalized=TRUE)



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




repair_mtx <- k_ctmtx[rownames(k_ctmtx) %in% gene_lists[[1]], ]
rt_resistance_mtx <- k_ctmtx[rownames(k_ctmtx) %in% gene_lists[[2]], ]
immune_mtx <- k_ctmtx[rownames(k_ctmtx) %in% gene_lists[[3]], ]

repair_mtx_norm <- norm.k_ctmtx[rownames(norm.k_ctmtx) %in% gene_lists[[1]], ]
rt_resistance_mtx_norm <- norm.k_ctmtx[rownames(norm.k_ctmtx) %in% gene_lists[[2]], ]
immune_mtx_norm  <- norm.k_ctmtx[rownames(norm.k_ctmtx) %in% gene_lists[[3]], ]


write.csv(repair_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/repair.csv")
write.csv(rt_resistance_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/rtr.csv")
write.csv(immune_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/immune.csv")

write.csv(repair_mtx_norm, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/repair_norm.csv")
write.csv(rt_resistance_mtx_norm, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/rtr_norm.csv")
write.csv(immune_mtx_norm, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/immune_norm.csv")
write.csv(norm.k_ctmtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/Ct_mtx_full_norm.csv")




#---------------------------Preparing Data for Rose-----------------------------
# CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv"
# "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/"
# ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)
# ct_mtx$ensembl_gene_id <- NULL
# # norm.ct_mtx <- ct_mtx
# 
# ## as DGEList
# dge <- DGEList(counts=ct_mtx)
# 
# ## calculate norm. factors
# dge <- calcNormFactors(dge)
# 
# ## get normalized counts
# norm.ct_mtx <- cpm(dge)

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

# Extract the normalized counts
norm.ct_mtx <- counts(dds, normalized=TRUE)


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

#--------------------subsetting data--------------------------------------

repair_mtx <- norm.ct_mtx[rownames(norm.ct_mtx) %in% gene_lists[[1]], ]
rt_resistance_mtx <- norm.ct_mtx[rownames(norm.ct_mtx) %in% gene_lists[[2]], ]
immune_mtx <- norm.ct_mtx[rownames(norm.ct_mtx) %in% gene_lists[[3]], ]


write.csv(repair_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/zscores/ForRose/Data/repair.csv")
write.csv(rt_resistance_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/zscores/ForRose/Data/RTResistance.csv")
write.csv(immune_mtx, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/zscores/ForRose/Data/immune.csv")








