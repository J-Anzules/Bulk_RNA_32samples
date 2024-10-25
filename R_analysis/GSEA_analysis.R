library(fgsea)
library(dplyr)
library(tidyverse)
##############################################################3
# Naming convention:
# res_[a]_[b]
# a b - these will be numbers representing the hours
# [a] =  hours of data that is RT
# [b] = *Should* always be 0 hour NRT
# 
# When analyzing gsea make sure to change all res_[a]_[a] for w.e. you are processing
##############################################################3

res <- read.csv("C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/data/results_dge.csv",
                    row.names = 1)

#---------------------- Data check and prep----------------------
# sum(is.na(res$padj))
# sum(is.na(res$pvalue))
# sum(is.na(res$lfcSE))
# sum(is.na(res$log2FoldChange))

# dim(res)
# dim(res)[1] - sum(is.na(res$pvalue))

# Dealing with NA
res <- na.omit(res)
# dim(res)

# sum(res == 0)

# Ensure pvalue is non-zero to avoid issues with -log10
res$pvalue[res$pvalue == 0] <- .Machine$double.xmin

#--------------------Preparing Ranked list----------------------------
res$rank <- res$log2FoldChange * -log10(res$pvalue) 
gene_ranks <- res$rank
names(gene_ranks) <- rownames(res)

# Remove any NA values that may have arisen
gene_ranks <- gene_ranks[!is.na(gene_ranks)]

# Sort the gene list in decreasing order
gene_ranks <- sort(gene_ranks, decreasing = TRUE)


#-------------------------GSEA-----------------------------------
#
# Loading GO pathways
pathways <- gmtPathways("C:/Users/jonan/Documents/1Work/RoseLab/reference_genomes/enrichment_analysis/c5.go.bp.v2024.1.Hs.symbols.gmt")

# Run fgsea for gene set enrichment analysis
fgsea_results <- fgsea(pathways = pathways, 
                       stats = gene_ranks, 
                       minSize = 15,   # Minimum size of gene sets
                       maxSize = 500,  # Maximum size of gene sets
                       nperm = 10000)  # Number of permutations

# Editing table
fgsea_results$direction <- ifelse(fgsea_results$NES > 0, "up", "down")
fgsea_results <- subset(fgsea_results, !is.na(NES))
fgsea_results$combined_score <- abs(fgsea_results$NES) * -log10(fgsea_results$pval)
fgsea_results <- fgsea_results[order(fgsea_results$combined_score, decreasing = TRUE)]
fgsea_results$pathway <- fgsea_results$pathway %>% str_replace("GO.+?_", "") %>% str_replace_all("_", " ")

# Convert the list in the 'leadingEdge' column to a character string (concatenated by commas)
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ", "))
# Generating df with only significant pathways
fgsea_results_sig <- subset(fgsea_results, pval < 0.05)

write.csv(fgsea_results_sig, "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/enrichment_results_significant.csv")
# TODO: AUTOMATION SPOT
# Write the dataframe to CSV
# write.csv(fgsea_results, "~/Documents/RoseLab/BulkRNA_seq/data/GSEA_results/first_16_samples/0_2.csv")
# TODO: I accidently saved the 0_2 sig results over the 0_0
# write.csv(fgsea_results_sig, "~/Documents/RoseLab/BulkRNA_seq/data/GSEA_results/first_16_samples/0_2_sig.csv")

# head(fgsea_results$combined_score, 20)
##------------------------Checking Results-----------------------------

# View the top enriched pathways
# fgsea_results <- fgsea_results[order(padj), ]
# fgsea_results <- fgsea_results[order(pval), ]
# 
# dim(fgsea_results)
# head(fgsea_results)

# Checking results stats
# sum(is.na(fgsea_results$padj))
# sum(fgsea_results$padj <= 0.10, na.rm = TRUE)
# sum(fgsea_results$pval <= 0.05, na.rm = TRUE)
# head(fgsea_results)


#---------------------------Finding Cell death and dna repair-----------------------
# Filter pathways containing keywords
# cell_death_pathways <- fgsea_results[grep("apoptosis|cell death|necrosis", fgsea_results$pathway, ignore.case = TRUE), ]
# dna_repair_pathways <- fgsea_results[grep("DNA repair|double-strand break repair|homologous recombination", fgsea_results$pathway, ignore.case = TRUE), ]
# 
# # View the filtered pathways
# cell_death_pathways
# dna_repair_pathways


#--------------------------Visualizing------------------------------------------
library(ggplot2)

##-------------------Plotting individual pathways ------------------------
# Plot the top enriched pathway
fgsea_results <- sort(fgsea_results$pval)
top_pathway <- fgsea_results$pathway[1]  # Adjust index as needed

plotEnrichment(pathways[[top_pathway]], gene_ranks) +
  labs(title = top_pathway)

##--------------------Plotting top up and down--------------------------


###---------------------Preparing top up and down---------------------------
top_count <- 15
topup <- fgsea_results %>% filter(ES > 0)
topup <- topup[order(topup$pval),]
topup <- topup[1:top_count,]

topdown <- fgsea_results %>% filter(ES < 0)
topdown <- topdown[order(topdown$pval),]
topdown <- topdown[1:top_count,]

top <- rbind(topup, rev(topdown))

# Clean up pathway names
# top$pathway <- top$pathway %>% str_replace("GO.+?_", "") %>% str_replace_all("_", " ")

top <- top[order(top$pval),]
top <- top %>% filter(pval<=0.05) %>% filter(!is.na(pathway))

#TODO: maybe, I should save all of the fgsea result?
# write.csv(top, "C:/Users/jonan/Documents/Tyseq/Data/ForReview/Top_updown_0_0.csv")

###---------------------Plotting---------------------------
# Summarize pathway information and create negative log p value variable for graphing
# top <- top %>% filter(pval <= 0.1)
pathg <- top %>% mutate(neglogpvalue = -log10(pval))

pathg$pathway <- tolower(pathg$pathway)

# Define a more visible color combination: coral for downregulated and skyblue for upregulated
pathfig <- ggplot(pathg, aes(x = reorder(pathway, neglogpvalue), y = neglogpvalue)) +
  # Conditionally set fill color based on the "ES" value
  geom_bar(aes(fill = ifelse(ES < 0, "Downregulated", "Upregulated")), stat = "identity") +
  coord_flip() +
  scale_x_discrete(name = "Pathways Associated with RT treatment") +
  
  # Use a manual color scale to assign custom colors to the levels
  scale_fill_manual(values = c("Downregulated" = "coral", "Upregulated" = "skyblue"),
                    name = "Regulation",  # Add a legend title
                    labels = c("Downregulated", "Upregulated")) +  # Legend labels
  ylab("-log(p value)") +
  
  # Customize theme for better appearance
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 0),
        axis.text.y = element_text(face = "bold", size = 10, angle = 0),
        legend.title = element_text(face = "bold"),
        legend.position = "right")  # Set legend position


# pathfig


# TODO: another place that needs to be edited for automation
ggsave(filename = "C:/Users/jonan/Documents/1Work/RoseLab/bulkRNAseq_32samples/Figures/From_unlceaned_data/top_up_n_down_pathways.png", 
       plot = pathfig,   # The name of your ggplot object
       width = 10,       # Width of the saved image (in inches)
       height = 6,       # Height of the saved image (in inches)
       dpi = 300)        # Resolution (300 dpi is high quality)



# TODO: I need to figure out a way to find the pathways that consistently appear between all the scenarios
# Once I have all the pathway results, write a function that identifies which pathway appears at least
# one other time in the other dataframes. If it appears in one other dataframe make a new dataframe,
# that has the following columns
# 1. pathway name
# 2. How many repetitions this pathway appears (minimum 2)
# 3. list of which test it appeared in
# 4. Maybe make a heatmap highlighting where the pathway appears and then let
#       heatmap represent either the rank that I made or p-value
#       maybe 
# 5. I could probably use plotEnrichment in a creative way
# 6. Do the same for the KEGG pathways
#################################################################3
#################################################################3
#################################################################3
#################################################################3
#----------------------------Scrap--------------------------------

