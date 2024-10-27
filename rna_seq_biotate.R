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


summary(res_ZT0)

#writing the test file 
write.csv(res_ZT0, "Deseq_results_ZT0.csv")
write.csv(res_ZT12, "Deseq_results_ZT12.csv")

df = read.csv("Deseq_results_ZT12.csv")

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


summary(df$padj)


#creating volcano plots

# Create a new column to categorize genes based on log2FoldChange
df$Regulation <- ifelse(df$log2FoldChange > 1, "Upregulated",
                        ifelse(df$log2FoldChange < -1, "Downregulated", "Not Significant"))

ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +  
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_x_continuous("log2(fold change)", 
                     breaks = seq(-10, 10, by = 1),
                     limits = c(-10, 10)) +  
  scale_y_continuous("-log10(p-value)", 
                     limits = c(0, max(-log10(df$padj), na.rm = TRUE))) +  
  theme_classic() +
  geom_label_repel(aes(label = Geneid), size = 3, 
                   box.padding = 0.5, point.padding = 0.3, force = 20,
                   max.overlaps = 50) +
  scale_color_manual(values = c("Upregulated" = "#FFA07A", "Downregulated" = "#00FFFF", "Not Significant" = "grey"))


# Filter significant genes (adjust this threshold as needed)
sig_genes <- df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Reshape data for heatmap
heatmap_data <- dcast(sig_genes, Geneid ~ Time_point + Tissue, value.var = "log2FoldChange")
heatmap_data[is.na(heatmap_data)] <- 0  # Replace NA with 0 for heatmap

# Melt the data for ggplot
heatmap_melted <- melt(heatmap_data, id.vars = "Geneid", variable.name = "Time_Tissue", value.name = "log2FoldChange")



# Assuming 'sig_genes' contains the significant genes and their log2FoldChange
heatmap_data <- dcast(sig_genes, Geneid ~ Tissue + Time_point, value.var = "log2FoldChange")
rownames(heatmap_data) <- heatmap_data$Geneid
heatmap_data <- heatmap_data[, -1]  # Remove Geneid column for heatmap

# Create the heatmap
pheatmap::pheatmap(heatmap_data, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)

