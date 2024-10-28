library(dplyr)
library(readr)
library(ggplot2)

#Unclean data
# path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_unclean/final.out/"

#Clean data
path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_clean_fast/final.log.out/"

file_list <- list.files(path, pattern = "Log.final.out$", full.names = TRUE)

# Function to read and process each file
process_file <- function(file) {
  # Read file as lines
  lines <- readLines(file)
  
  # Split lines into name and value parts, assuming each line is "MetricName\tMetricValue"
  data <- sapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) == 2) {
      parts[2]  # Extract the metric value
    } else {
      NA  # Handle any irregular lines
    }
  })
  
  # Create a named vector for easy conversion to a data frame row
  names(data) <- sapply(lines, function(line) strsplit(line, "\t")[[1]][1])
  data <- c(FileName = basename(file), data)  # Add the file name
  
  # Convert to data frame with one row
  as.data.frame(t(data), stringsAsFactors = FALSE)
}

# Apply the function to each file and combine results
df <- do.call(rbind, lapply(file_list, process_file))

# Display the dataframe
df <- type_convert(df)  # Optional: Convert columns to appropriate types

# Clean column names
clean_names <- function(names) {
  names <- gsub("\\s+", " ", names)       # Replace multiple spaces with a single space
  names <- gsub("\\|", "", names)         # Remove |
  names <- gsub("%", "percent", names)    # Replace % with 'percent'
  names <- gsub(":", "", names)           # Remove colons
  names <- gsub("^\\s+|\\s+$", "", names) # Trim leading and trailing whitespace
  names <- make.names(names, unique = TRUE) # Make syntactically valid names
  return(names)
}

# Apply this function to the data frame column names
colnames(df) <- clean_names(colnames(df))

# Remove "Log.final.out" from each entry in the FileName column
# df$FileName <- sub("Log.final.out$", "", df$FileName)
df$FileName <- sub("_Log.final.out$", "", df$FileName)
#----------------adding actual file name---------------------------#

# Load the metadata
metadata <- read.csv("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/metadata_short.csv")

# Rename the 'FileName' column in 'df' to 'TGen_Sample_ID' for merging
df <- merge(df, metadata[, c("TGen_Sample_ID", "Sample_ID")], 
            by.x = "FileName", by.y = "TGen_Sample_ID", all.x = TRUE)

# Rename the 'Sample_ID' column to 'SampleID' in the merged df
names(df)[names(df) == "Sample_ID"] <- "SampleID"

# df$Uniquely.mapped.reads.percent
# 

#--------------------------Histogram of percent unqiue------------------------#
# Remove the "%" symbol and convert to numeric
df$Uniquely.mapped.reads.percent <- as.numeric(sub("%", "", df$Uniquely.mapped.reads.percent))

# write.csv(df, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/final.log.out/unclean_alignment_data.csv")
write.csv(df, "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_clean_fast/clean_alignment_data.csv")
# Define the SampleID values to keep
# sample_ids_to_plot <- c(
#   "COH2233-23_Benign_NRT_0h", 
#   "COH2233-23_Benign_RT_0h", 
#   "COH22333-24_A+B_NRT_0h", 
#   "COH22333-24_A+B_RT_0h", 
#   "COH2233-22 Benign NRT 0h", 
#   "COH2233-22 Tumor RT 2h"
# )
# 
# # Subset the data frame to only include these SampleID values
# df <- df %>% filter(SampleID %in% sample_ids_to_plot)
df$

plot <- ggplot(df, aes(x = SampleID, y = Uniquely.mapped.reads.percent)) +
  geom_bar(stat = "identity") +
  theme_bw() +  # White background
  labs(
    title = "Uniquely Mapped Reads Percent per Sample",
    x = "Sample ID",
    y = "Uniquely Mapped Reads (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  )+
  ylim(0, 100)
  

# Specify the output path
# output_path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/final.log.out/unclean_unique_read_pct.png"
# output_path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/STAR_clean_fast/Figures/clean_unique_read_pct.png"
output_path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/fasp_star_clean/clean_unique_read_pct.png"
# output_path <- "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/final.log.out/unclean_subset_read_pct.png"
# Save the plot
ggsave(output_path, plot = plot, width = 10, height = 6, dpi = 300)


