
#Files location: "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_unclean/ReadsPerGene/"

#########################################
# Removing the first 4 lines from the ReadsPerGene count matrix
# Selecting the first column with the unstranded counts
# naming each sample/column based on it's Anumbers

library(biomaRt) # Converting ensembleID to geneID

#----------------------Loading and Preparing Count matrix----------------------
# Function to load the second column (unstranded counts) from each file
load_counts <- function(file) {
  # Extract the file name to be used as the column name
  sample_id <- sub("ReadsPerGene.out.tab", "", basename(file))
  
  # Read the data, skip the first 4 lines, use the first column as rownames
  data <- read.table(file, header = FALSE, skip = 4, row.names = 1)
  
  # Return only the unstranded counts (second column) as a data frame
  unstranded_counts <- data.frame(data[, 1])
  colnames(unstranded_counts) <- sample_id  # Rename the column with sample ID
  
  return(unstranded_counts)
}

# Specify the folder containing the files
folder_path <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_unclean/ReadsPerGene/"

# List all files in the folder that match the pattern "ReadsPerGene.out.tab"
files <- list.files(folder_path, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)

# Load the counts for all files
all_counts <- lapply(files, load_counts)

# Since rows are the same, bind them together by columns
combined_data <- do.call(cbind, all_counts)

# Use the rownames from one file for the gene IDs
gene_ids <- rownames(read.table(files[1], header = FALSE, skip = 4, row.names = 1))
rownames(combined_data) <- gene_ids


#----------------------Getting gene symbols----------------------

# Preparing Ensembl IDs (removing version numbers as before)
ensemble_ids <- gsub("\\.\\d+$", "", rownames(combined_data))
rownames(combined_data) <- gsub("\\.\\d+$", "", rownames(combined_data))

combined_data <- data.frame(ensembl_gene_id = rownames(combined_data), combined_data)
combined_data <- combined_data[!duplicated(combined_data$ensembl_gene_id), ]

# Use biomaRt to map Ensembl IDs to gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mappings from Ensembl IDs to gene symbols
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "ensembl_gene_id", 
  values = ensemble_ids, 
  mart = ensembl
)

# Create a mapping from Ensembl IDs to Gene Names
id_to_name <- setNames(gene_mapping$hgnc_symbol, gene_mapping$ensembl_gene_id)

# Replace row names in combined_data with gene names, handling duplicates by appending '_dpl'
combined_data_rownames <- rownames(combined_data)

# Initialize the counting variable
loop_counter <- 0

for (i in seq_along(combined_data_rownames)) {
  ensembl_id <- combined_data_rownames[i]
  
  # Get the corresponding gene symbol or keep NA if not found
  gene_symbol <- ifelse(ensembl_id %in% names(id_to_name), id_to_name[ensembl_id], NA)
  
  # Increase the counter for each loop iteration
  loop_counter <- loop_counter + 1
  
  # Every 500 iterations, print the gene symbol or Ensembl ID and the loop count
  if (loop_counter %% 500 == 0) {
    message("On iteration ", loop_counter, " with Ensembl ID: ", ensembl_id, 
            " and Gene Symbol: ", gene_symbol)
  }
  
  # If the gene symbol is NA, keep the original Ensembl ID
  if (is.na(gene_symbol)) {
    final_rownames[i] <- ensembl_id
  } else {
    # If this gene symbol has already been seen, append a suffix
    if (gene_symbol %in% gene_count) {
      gene_count[[gene_symbol]] <- gene_count[[gene_symbol]] + 1
      final_rownames[i] <- paste0(gene_symbol, "_dpl", gene_count[[gene_symbol]])
    } else {
      # If it's the first time seeing this gene symbol, just assign it
      gene_count[[gene_symbol]] <- 1
      final_rownames[i] <- gene_symbol
    }
  }
}
write.csv(final_rownames_df, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/final_rownames.csv", row.names = FALSE)



# Create a data frame to save the final_rownames (optional, if it's not already a data frame)
final_rownames_df <- data.frame(rownames = final_rownames)

# Replace the empty entries with the corresponding Ensembl IDs from the original combined_data_rownames
final_rownames[empty_indices] <- combined_data_rownames[empty_indices]


# Find the indices of the duplicates in final_rownames
duplicate_indices <- which(duplicated(final_rownames))

# Loop through each duplicate index and append "_dpl" to make them unique
for (i in duplicate_indices) {
  # Add "_dpl" to the duplicate row name to make it unique
  final_rownames[i] <- paste0(final_rownames[i], "_dpl")
}

# Set the updated row names back to the combined_data
rownames(combined_data) <- final_rownames

write.csv(combined_data, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_unclean.csv")

