# Load necessary library
library(dplyr)

# File paths
CountMatrix_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_clean_fpStr.csv"
Metadata_loc <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/metadata_short.csv"

# Load count matrix
ct_mtx <- read.csv(CountMatrix_loc, header = TRUE, row.names = 1)

dim(ct_mtx)
# Load metadata
metadata <- read.csv(Metadata_loc)

# Define the samples to remove
samples_to_remove <- c(
  "COH22333-34_A_RT_0h",
  "COH22333-34_A_NRT_0h",
  "COH2233-23_Benign_RT_0h",
  "COH2233-23_Benign_NRT_0h",
  "COH2233-22 Benign RT 0h",
  "COH2233-22 Benign NRT 0h"
)

# Identify corresponding TGen_Sample_IDs from metadata
tgen_ids_to_remove <- metadata %>%
  filter(Sample_ID %in% samples_to_remove) %>%
  pull(TGen_Sample_ID)

# Filter count matrix (remove unwanted samples)
ct_mtx_filtered <- ct_mtx[, !(colnames(ct_mtx) %in% tgen_ids_to_remove)]

# Filter metadata (keep only relevant samples)
metadata_filtered <- metadata %>%
  filter(!(TGen_Sample_ID %in% tgen_ids_to_remove))

# Ensure metadata order matches the count matrix columns
metadata_filtered <- metadata_filtered %>%
  filter(TGen_Sample_ID %in% colnames(ct_mtx_filtered)) %>%
  arrange(match(TGen_Sample_ID, colnames(ct_mtx_filtered)))

# Verify if everything matches
if (!all(metadata_filtered$TGen_Sample_ID == colnames(ct_mtx_filtered))) {
  stop("Error: Metadata and count matrix column names do not match!")
}

# Define output directory
output_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/filtered_data_for_Helya"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the new count matrix (unnormalized)
write.csv(ct_mtx_filtered, file = file.path(output_dir, "count_matrix_filtered.csv"), row.names = TRUE)

# Save the new metadata
write.csv(metadata_filtered, file = file.path(output_dir, "metadata_filtered.csv"), row.names = FALSE)

cat("Filtered count matrix and metadata successfully saved to:", output_dir, "\n")
