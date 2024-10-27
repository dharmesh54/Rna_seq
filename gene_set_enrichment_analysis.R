setwd("C:/Users/dharm/Documents/Biostate/Deseq2/")

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(biomaRt)
library(dplyr)

#reading csv files
df_ZT0 = read.csv("Deseq_results_ZT0.csv")
df_ZT12 = read.csv("Deseq_results_ZT12.csv")


colnames(df_ZT0)[1] <- "Geneid"
colnames(df_ZT12)[1] <- "Geneid"

# Apply the filter for padj <= 0.05 and log2FoldChange thresholds
filtered_ZT0 <- subset(df_ZT0, padj <= 0.05 & (log2FoldChange > 0.6 | log2FoldChange < -0.6))
filtered_ZT12 <- subset(df_ZT12, padj <= 0.05 & (log2FoldChange > 0.6 | log2FoldChange < -0.6))

# Define the function
perform_enrichment <- function(gene_symbols, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10) {
  
  # Convert gene symbols to Entrez IDs
  converted_genes <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Filter to include only matched Entrez IDs
  filtered_genes <- converted_genes[!is.na(converted_genes$ENTREZID), ]
  
  # Check if there are any matched genes
  if (nrow(filtered_genes) == 0) {
    message("No genes could be mapped to Entrez IDs.")
    return(NULL)
  }
  
  # Run the enrichment analysis using the converted Entrez IDs
  go_enrich <- enrichGO(
    gene = filtered_genes$ENTREZID,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    readable = TRUE,
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff
  )
  
  # Return the enrichment results
  return(as.data.frame(go_enrich@result))
}

gene_symbols <- filtered_ZT0$Geneid  # Replace this with your gene symbols
enrichment_result_0 <- perform_enrichment(gene_symbols)
head(enrichment_results)


enrichment_result_0$Time_points = "ZT0"
enrichment_result_12$Time_points = "ZT12"


go_enrich <- rbind2(enrichment_result_0, enrichment_result_12)
go_enrich <- na.omit(go_enrich)
head(go_enrich)


# Sort the go_enrich data frame based on p.adjust in ascending order
sorted_go_enrich <- go_enrich[order(go_enrich$p.adjust), ]

# Select the top 40 or 50 pathways
top_pathways <- head(sorted_go_enrich, 50)  # Change 50 to 40 if needed


# Create a dot plot for the top pathways
ggplot(top_pathways, aes(x = Time_points, y = reorder(Description, -log10(p.adjust)), size = -log10(p.adjust))) +
  geom_point(color = "lightcoral") +  # Use light red color for points
  scale_size(range = c(2, 10)) +  # Adjust the size range of the points
  labs(x = "Time Points", y = "Pathway", title = "Top Pathways Dot Map") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme_minimal()
