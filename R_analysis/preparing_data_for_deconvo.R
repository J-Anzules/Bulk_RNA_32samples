data <- read.table("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/edgeR_ct_mtx.tsv", sep="\t", header=TRUE, row.names=1)
# str(data)
# 
# head(data)
# any(is.na(colnames(data)))

# colnames(norm.ct_mtx) <- gsub(" ", "_", colnames(data))
# colnames(norm.ct_mtx) <- gsub("\\.", "_", colnames(data))
# write.table(norm.ct_mtx, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/deseq_ct_mtx_clean.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Replace spaces with underscores in column names
colnames(data) <- gsub(" ", "_", colnames(data))
# Replace periods with underscores in column names
colnames(data) <- gsub("\\.", "_", colnames(data))
data$COH2233_22_Benign_RT_0h <- NULL
# colnames(data)[colnames(data) == "COH2233_22_Benign_RT_0h"] <- "22RT0h"#"COH2233_22BenignRT0h"

write.table(data, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/edgeR_ct_mtx_clean_no22.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# 
# colnames(data)
# subset_data <- data[, 1:22]
# colnames(subset_data)[colnames(subset_data) == "COH2233_22_Benign_RT_0h"] <- "COH2233_22_Benign_RT_0h_cleaned"
# write.table(subset_data, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/edgeR_ct_mtx_clean_subset.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# 
# subset_data <- data[, 23:30]
# dim(subset_data)
# colnames(data)[22]
# 
# # colnames(data)[colnames(data) == "COH2233_22_Benign_RT_0h"] <- "COH2233_22_Benign_RT_0h_cleaned"
# 
# write.table(subset_data, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/edgeR_ct_mtx_clean_subset_23to30.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# 
# 
# 
# #-------------------------------What's wrong with 22?------------------------
# 
# identical(colnames(data)[which(colnames(data) == "COH2233_22_Benign_RT_0h")], "COH2233_22_Benign_RT_0h")


#------------------------------Kavyas sample-------------------------

kdata <- read.table("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/counts_without_duplicates.txt", sep = "\t", header = TRUE, row.names = 1)
kdata$ID.1 <- NULL

any(is.na(kdata))
kdata[is.na(kdata)] <- 0

colnames(kdata) <- gsub(" ", "_", colnames(kdata))
# Replace periods with underscores in column names
colnames(kdata) <- gsub("\\.", "_", colnames(kdata))


colnames(kdata)
write.table(kdata, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq/data/kavya_data_noNA.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#-------------------------Generating Figure------------------------------
library(ggplot2)
library(reshape2)

data <- read.csv("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/for_deconvolution/KDATA.csv")
data$RMSE <- NULL
data$P.value <- NULL
data$Correlation <- NULL

colnames(data)
data$Mixture
# Melt the data to make it long format for ggplot2
data_long <- melt(data, id.vars = "Mixture", variable.name = "Cell_Type", value.name = "Proportion")


# Custom color palette for each cell type
# Custom color palette matching your immune cell names
custom_colors <- c(
  "B.cells.naive" = "#FFFF99",                 # Light yellow
  "B.cells.memory" = "#FFCC00",                # Gold
  "Plasma.cells" = "#FFFF66",                  # Yellow
  "T.cells.CD8" = "#FFD9B3",                   # Light peach
  "T.cells.CD4.naive" = "#D9B3FF",             # Light purple
  "T.cells.CD4.memory.resting" = "#B266FF",    # Medium purple
  "T.cells.CD4.memory.activated" = "#9900FF",  # Dark purple
  "T.cells.follicular.helper" = "#FF66FF",     # Pink
  "T.cells.regulatory..Tregs." = "#CC0099",    # Magenta
  "T.cells.gamma.delta" = "#FFB3E6",           # Light pink
  "NK.cells.resting" = "#B3A0CC",              # Light purple-grey
  "NK.cells.activated" = "#6600CC",            # Dark purple
  "Monocytes" = "#FF9999",                     # Light red
  "Macrophages.M0" = "#FFCC99",                # Light orange
  "Macrophages.M1" = "#FFB266",                # Medium orange
  "Macrophages.M2" = "#FF8000",                # Dark orange
  "Dendritic.cells.resting" = "#E6B3B3",       # Light red-brown
  "Dendritic.cells.activated" = "#FF4D4D",     # Red
  "Mast.cells.resting" = "#BFBFBF",            # Light grey
  "Mast.cells.activated" = "#E6E6E6",          # Very light grey
  "Eosinophils" = "#8C8C8C",                   # Medium grey
  "Neutrophils" = "#4D4D4D"                    # Dark grey
)



deconv_fig <- ggplot(data_long, aes(x = Mixture, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = custom_colors) +   # Apply custom color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),   # White background for plot area
    plot.background = element_rect(fill = "white", color = NA),    # White background for entire plot
    panel.grid.major = element_line(color = "gray90"),             # Light gray grid lines for visibility
    panel.grid.minor = element_line(color = "gray90")
  )


# Save the ggplot figure
ggsave(
  filename = "kavyas_16_samples.png",    # File name
  plot = deconv_fig,                          # Save the last created plot
  path = "C://Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/deconvolution/",             # Directory to save in
  width = 16,                                  # Width in inches (adjust as needed)
  height = 8,                                  # Height in inches (adjust as needed)
  units = "in",                                # Measurement unit
  dpi = 300                                    # Resolution in dots per inch
)


