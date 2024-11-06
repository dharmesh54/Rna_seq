setwd("C:/Users/dharm/Documents/Biostate/Deseq2/")
getwd()

#Libraries  
library(dplyr)
library("readxl")
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape2)


#importing merged count file generated earlier
data <- read.csv("merged_counts.txt", header = TRUE, sep = "\t")





#filtering data that zeroes throughtout all samples
dim(data)
head(data)
data[,2:9]
k <- rowSums(data[,2:9]) >= 1
data <- data[k,]
dim(data)



#indexing gene column
row.names(data) <- data$Geneid

#preparing design data
sample <- c("Heart_ZT0_rep1", "Heart_ZT0_rep2",
            "Heart_ZT12_rep1", "Heart_ZT12_rep2",
            "Liver_ZT0_rep1", "Liver_ZT0_rep2",
            "Liver_ZT12_rep1", "Liver_ZT12_rep2")

# Tissue types
tissue <- c(rep("Heart", 4), rep("Liver", 4))

# Time points
time <- c(rep("ZT0", 2), rep("ZT12", 2), rep("ZT0", 2), rep("ZT12", 2))

# Biological replicates
replicate <- c("rep1", "rep2", "rep1", "rep2", "rep1", "rep2", "rep1", "rep2")

# Creating a data frame for the design table
design <- data.frame(sample, tissue, time, replicate)

design$tissue <- as.factor(design$tissue)
design$time <- as.factor(design$time)
design$replicate <- as.factor(design$replicate)

#eliminating the GeneId column as we have indexed already 
my_data <- select(data, -Geneid)
my_data 

design

#generating dataset (dds) from count data and design metadata
dds <- DESeqDataSetFromMatrix(countData = my_data, colData = design, design = ~tissue + time + tissue:time)
dim(dds)

colData(dds)

#running DeSeq
ddsDE <- DESeq(dds)

#storing result
#res0.05 <- results(ddsDE, alpha = 0.05)
#summary(res)

resultsNames(ddsDE)



# Contrast for Heart vs. Liver at time point 0
#res_time0 <- results(ddsDE, contrast = c("tissue", "Heart:0", "Liver:0"))

res_ZT0 = results(ddsDE, name="tissue_Liver_vs_Heart")
res_ZT12 = results(ddsDE, contrast=list(c("tissue_Liver_vs_Heart","tissueLiver.timeZT12")), test="Wald")
res_time_specific = results(ddsDE, name="time_ZT12_vs_ZT0")


summary(res_ZT0)

#writing the test file 
write.csv(res_ZT0, "Deseq_results_ZT0.csv")
write.csv(res_ZT12, "Deseq_results_ZT12.csv")

df_ZT0 = read.csv("Deseq_results_ZT0.csv")
df_ZT12 = read.csv("Deseq_results_ZT12.csv")
df_time_specific = read.csv("Deseq_results_Time_specofic.csv")




create_volcano_plot <- function(df) {
  # Rename the first column to "Geneid"
  colnames(df)[1] <- "Geneid"
  
  # Filter out rows with non-finite values in log2FoldChange or padj
  df <- df[is.finite(df$log2FoldChange) & is.finite(df$padj), ]
  df <- df[df$log2FoldChange >= -10 & df$log2FoldChange <= 10, ]
  df <- df[!is.na(df$padj), ]
  df <- df[!is.na(df$log2FoldChange), ]
  
  # Replace zero or very small padj values to avoid Inf in -log10 transformation
  df$padj[df$padj < .Machine$double.eps] <- .Machine$double.eps
  
  # Set a reasonable upper limit for y-axis
  y_limit <- min(50, max(-log10(df$padj), na.rm = TRUE))  # Capping y-axis at 50
  
  # Create a new column to categorize genes based on log2FoldChange
  df$Regulation <- ifelse(df$log2FoldChange > 1, "Upregulated",
                          ifelse(df$log2FoldChange < -1, "Downregulated", "Not Significant"))
  
  # Create the volcano plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    scale_x_continuous("log2(fold change)", 
                       breaks = seq(-10, 10, by = 1),
                       limits = c(-10, 10)) +
    scale_y_continuous("-log10(p-value)", limits = c(0, y_limit)) +
    theme_classic() +
    geom_label_repel(aes(label = Geneid), size = 3, 
                     box.padding = 0.5, point.padding = 0.3, force = 20,
                     max.overlaps = 50) +
    scale_color_manual(values = c("Upregulated" = "#FFA07A", "Downregulated" = "#00FFFF", "Not Significant" = "grey"))
  
  # Print summary of padj values
  print(summary(df$padj))
  
  # Return the plot object
  return(p)
}

#Calling the function
p <- create_volcano_plot(df_time_specific)

p

#Geneating heatmaps

# Extract normalized counts
norm_counts <- counts(ddsDE, normalized = TRUE)
write.csv(norm_counts, "normalised_count.csv")
norm_counts


# Define significance threshold for DEGs
alpha <- 0.05

# Filter significant DEGs
sig_genes_ZT0 <- rownames(res_ZT0)[which(res_ZT0$padj < alpha)]
sig_genes_ZT12 <- rownames(res_ZT12)[which(res_ZT12$padj < alpha)]

# Define annotation based on tissue types
sample_info <- as.data.frame(colData(ddsDE)[, c("tissue", "time")])

# Check that sample_info contains the correct 'time' values
unique(sample_info$time)


# Subset normalized counts for DEGs at each time point
deg_data_ZT0 <- norm_counts[sig_genes_ZT0, sample_info$time == "ZT0"]
deg_data_ZT12 <- norm_counts[sig_genes_ZT12, sample_info$time == "ZT12"]

# Define heatmap color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Heatmap for ZT0
pheatmap(deg_data_ZT0,
         color = heatmap_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sample_info[sample_info$time == "ZT0",],
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "ZT0 Time Point")

# Heatmap for ZT12
pheatmap(deg_data_ZT12,
         color = heatmap_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sample_info[sample_info$time == "ZT12",],
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "ZT12 Time Point")


#Generating heatmap for the top regulated genes  

#Reading normalised counts for generating heatmap to show the pattern of top regulated genes 
df2 = read.csv("normalised_count.csv")

# Rename the first column to "Geneid"
colnames(df2)[1] <- "Geneid"
colnames(df_ZT0)[1] <- "Geneid"
colnames(df_ZT12)[1] <- "Geneid"

# 1. Filter for significant genes from both datasets
significant_genes_ZT0 <- df_ZT0[df_ZT0$padj < 0.05, ]
significant_genes_ZT12 <- df_ZT12[df_ZT12$padj < 0.05, ]

# 2. Select top 30 DEGs based on absolute log2FoldChange for both datasets
top_genes_ZT0 <- head(significant_genes_ZT0[order(-abs(significant_genes_ZT0$log2FoldChange)), ], 20)
top_genes_ZT12 <- head(significant_genes_ZT12[order(-abs(significant_genes_ZT12$log2FoldChange)), ], 20)



# 3. Combine the gene lists
top_genes_combined <- unique(c(top_genes_ZT0$Geneid, top_genes_ZT12$Geneid))
top_genes_combined

# 4. Filter the normalized counts dataset based on the top genes
# Assuming normalized_counts is your normalized counts dataframe with 'Geneid' as the first column
filtered_counts <- df2[df2$Geneid %in% top_genes_combined, ]

# Reshape the DataFrame
melted_df <- melt(filtered_counts, id.vars = "Geneid", variable.name = "Sample", 
                  value.name = "value")

melted_df

melted_df <- melted_df %>%
  mutate(Sample = as.character(Sample))


DEGS_pattern <- melted_df %>%
  mutate(Tissue = sapply(strsplit(Sample, "_"), `[`, 1),
         Time = sapply(strsplit(Sample, "_"), `[`, 2))

DEGS_pattern <- DEGS_pattern %>%
  mutate(log2_value = log2(value + 1))  # Add 1 to avoid -Inf for zeros

line_samples <- c("Heart", "Liver")

line_positions <- match(line_samples, unique(DEGS_pattern$Tissue))

max(DEGS_pattern$log2_value, na.rm = TRUE)

# Calculate the minimum value of the 'value' column
min(DEGS_pattern$log2_value, na.rm = TRUE)


# Define custom color palette
colors <- c('#40E0D0', 'white', 'red')  # Color for the values
under_color <- '#00ced1'  # Color for values below the minimum
over_color <- '#8B0000'   # Color for values above the maximum

DEGS_pattern %>%
  ggplot(aes(x = Tissue, y = Geneid, fill = log2_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Adding vertical lines for the specified samples and terminal ends
  geom_vline(xintercept = c(0.5, line_positions - 0.5, 
                            line_positions + 0.5, 
                            length(unique(DEGS_pattern$Tissue)) + 0.5),
             linetype = "solid", color = "black", size = 0.1) +
  # Adding horizontal lines at the top and bottom
  geom_hline(yintercept = c(0.5, length(unique(DEGS_pattern)) + 0.5),
             linetype = "solid", color = "black", size = 0.1) +
  # Rotate the x-axis labels by 90 degrees
  theme(axis.text.x = element_text(angle = 90, hjust = 1, family = "Arial", size = 16), 
        axis.text.y = element_text(family = "Arial", size = 6), 
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  facet_wrap(~ Time, nrow = 1) 
  #scale_fill_gradientn(
    #colors = c(under_color, colors, over_color),
    #values = rescale(c(-6, -5, 0, 5, 6)),
    #limits = c(-5, 5),
    #oob = scales::oob_squish,
    #guide = "colorbar"
  #)



# Improved ggplot code for better color representation 
DEGS_pattern %>%
  ggplot(aes(x = Tissue, y = Geneid, fill = log2_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  
  # Adding vertical lines for specified samples and terminal ends
  geom_vline(xintercept = c(0.5, line_positions - 0.5, 
                            line_positions + 0.5, 
                            length(unique(DEGS_pattern$Tissue)) + 0.5),
             linetype = "solid", color = "black", size = 0.1) +
  # Adding horizontal lines at the top and bottom
  geom_hline(yintercept = c(0.5, length(unique(DEGS_pattern$Geneid)) + 0.5),
             linetype = "solid", color = "black", size = 0.1) +
  
  # Applying a color gradient for abundance data
  scale_fill_gradientn(
    colors = c("lightyellow", "orange", "red", "darkred"),
    values = scales::rescale(c(0, 5, 10, 15)),  # Adjust range based on data
    limits = c(0, 20),  # Adjust max limit as needed
    name = "log2(value)",  # Color bar label
    oob = scales::oob_squish  # Squish values outside of limits to min/max colors
  ) +
  
  # Theme and label adjustments for clarity
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, family = "Arial", size = 14), 
    axis.text.y = element_text(family = "Arial", size = 7),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank()  # Remove grid for cleaner look
  ) +
  
  # Adding facets for different time points
  facet_wrap(~ Time, nrow = 1)

