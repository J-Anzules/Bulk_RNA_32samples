
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
  sample_id <- sub("_ReadsPerGene.out.tab", "", basename(file))
  
  # print("sample ID")
  # print(sample_id)
  # Read the data, skip the first 4 lines, use the first column as rownames
  data <- read.table(file, header = FALSE, skip = 4, row.names = 1)
  
  # Return only the unstranded counts (second column) as a data frame
  unstranded_counts <- data.frame(data[, 1])
  colnames(unstranded_counts) <- sample_id  # Rename the column with sample ID
  
  return(unstranded_counts)
}

# Specify the folder containing the files
folder_path <- "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_clean_fast/ReadsPerGene/"

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

# Preparing Ensembl IDs 
ensemble_ids <- gsub("\\.\\d+$", "", rownames(combined_data))
rownames(combined_data) <- gsub("\\.\\d+$", "", rownames(combined_data))


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

# Add a new column 'gene_symbols' based on row names and id_to_name mapping
combined_data$gene_symbols <- id_to_name[rownames(combined_data)]

#--------------------------Cleaning dataframe-----------------------------


# head(combined_data$gene_symbols, n = 100)



# Removing empty rows and duplicates
dim(combined_data) #58721    
combined_data <- combined_data[!(is.na(combined_data$gene_symbols) | combined_data$gene_symbols == ""), ]
dim(combined_data) #40882    
combined_data <- combined_data[!duplicated(combined_data$gene_symbols), ]
dim(combined_data) #40865    

#renaming the rows
rownames(combined_data) <- combined_data$gene_symbols

# Not needed anymore
combined_data$gene_symbols <- NULL

# Removing low quality reads
combined_data <- combined_data[rowSums(combined_data >= 30) >= 6, ]



#-------------Writing the results-------------------


write.csv(combined_data, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/count_matrix_clean_fpStr.csv")

