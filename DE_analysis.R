# /media/ekramah/DISK2/RNASEQ2/07_scripts/DE_analysis.R

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# --- Set working directory ---
setwd("/media/ekramah/DISK2/RNASEQ2")
cat(sprintf("Working directory set to: %s\n", getwd()))

# --- 1. Load data ---
cat("Loading gene counts data...\n")
count_data <- read.table("04_quantification/gene_counts.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
cat(sprintf("Raw data loaded. Dimensions: %d rows, %d columns.\n", nrow(count_data), ncol(count_data)))

# Remove the first 6 metadata columns from featureCounts output
count_data <- count_data[, 7:ncol(count_data)]
cat(sprintf("Data after removing metadata columns. Dimensions: %d rows, %d columns.\n", nrow(count_data), ncol(count_data)))

# --- Robust Column Name Cleaning ---
original_colnames <- colnames(count_data)
cat("Original column names from featureCounts (before cleaning):\n")
print(original_colnames)

# NEW ROBUST CLEANING STRATEGY: Extract only the SRRxxxxx pattern
cleaned_colnames <- sub(".*(SRR[0-9]{7}).*", "\\1", original_colnames)

# Assign the cleaned names back to the count_data
colnames(count_data) <- cleaned_colnames
cat("Cleaned column names.\n")

# --- DEBUGGING: Check cleaned column names for duplicates or issues ---
cat("Final cleaned sample names (column names of count_data after robust cleaning):\n")
print(colnames(count_data))
if (any(duplicated(colnames(count_data)))) {
  stop("CRITICAL ERROR: Duplicate sample names found after cleaning column names. Please check your featureCounts output and naming convention.")
}
cat("No duplicates found in cleaned column names. Proceeding...\n")

# Ensure count data is in integer mode (DESeq2 requires integers)
count_data <- round(count_data)
cat("Count data rounded to integers.\n")


# --- 2. Create sample metadata (colData) ---
sample_names_in_counts <- colnames(count_data)

# Define a mapping from sample accession to condition
conditions_map <- data.frame(
  sample = c("SRR7498002", "SRR7498004", "SRR7498006", "SRR7498008", "SRR7498010", # WildType
             "SRR7498016", "SRR7498018", "SRR7498020", "SRR7498022", "SRR7498024"), # Notch2CKO
  condition = c(rep("WildType", 5), rep("Notch2CKO", 5))
)
cat("Conditions map defined.\n")
cat("Expected samples in conditions_map:\n")
print(conditions_map$sample)

# Create colData by matching samples from count_data to conditions_map
cat("Attempting to create colData using match...\n")
match_indices <- match(sample_names_in_counts, conditions_map$sample)
cat("Result of match (indices):\n")
print(match_indices)

# Check for NAs in match_indices - this means some sample names from count_data were not found in conditions_map
if (any(is.na(match_indices))) {
  missing_samples <- sample_names_in_counts[is.na(match_indices)]
  stop(sprintf("CRITICAL ERROR: The following samples from your count data were not found in your conditions_map: %s. Please ensure exact matching.", paste(missing_samples, collapse=", ")))
}

colData <- conditions_map[match_indices, ] # Subset conditions_map based on matched order
cat("colData created before setting row names. colData$sample values:\n")
print(colData$sample)

# Check for duplicates in colData$sample BEFORE assigning to row.names
if (any(duplicated(colData$sample))) {
  duplicated_values <- colData$sample[duplicated(colData$sample)]
  stop(sprintf("CRITICAL ERROR: Duplicate values found in colData$sample before assigning row names: %s. This indicates an issue with sample name uniqueness after matching.", paste(duplicated_values, collapse=", ")))
}

rownames(colData) <- colData$sample # Set row names to sample accessions
cat("colData row names set successfully.\n")

# Final check: ensure all samples have a condition
if (any(is.na(colData$condition))) {
  stop("Error: Not all samples in colData have a valid condition after matching.")
}
cat("colData created successfully.\n")

# Convert 'condition' column to a factor, which DESeq2 requires
colData$condition <- factor(colData$condition)
# Set the reference level for comparison (e.g., 'WildType' as the baseline)
colData$condition <- relevel(colData$condition, ref = "WildType")
cat("Conditions set as factors with 'WildType' as reference.\n")


# --- 3. Construct DESeqDataSet object ---
cat("Constructing DESeqDataSet object...\n")
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)
cat("DESeqDataSet object created.\n")

# --- 4. Pre-filtering: Remove genes with very low counts ---
# Keep only genes with at least 10 reads in total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat(sprintf("Retained %d genes after filtering (out of %d total).\n", nrow(dds), nrow(count_data)))


# --- 5. Run DESeq2 analysis ---
cat("Running DESeq analysis...\n")
dds <- DESeq(dds)
cat("DESeq2 analysis completed.\n")

# --- DEBUGGING: Print available coefficients for lfcShrink ---
cat("Available coefficients for lfcShrink (resultsNames(dds)):\n")
print(resultsNames(dds))
# Expected: "Intercept", "condition_Notch2CKO_vs_WildType"


# --- 6. Get results ---
cat("Extracting results...\n")
# Get results for the comparison of 'Notch2CKO' vs 'WildType'
res <- results(dds, contrast=c("condition","Notch2CKO","WildType"))

# Apply Log2 Fold Change (LFC) shrinkage
# IMPORTANT FIX: Use 'coef' argument for 'apeglm' type shrinkage
lfc_coeff_name <- "condition_Notch2CKO_vs_WildType"
cat(sprintf("Applying apeglm shrinkage using coefficient: %s\n", lfc_coeff_name))
res <- lfcShrink(dds, coef=lfc_coeff_name, res=res, type="apeglm")

# Order the results by adjusted p-value (padj)
res_ordered <- res[order(res$padj),] # *** CORRECTED LINE: from res_ordered$padj to res$padj ***

# Convert to a data frame and add gene names
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)
res_df$is_significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, "Yes", "No")
cat(sprintf("Found %d significantly differentially expressed genes (padj < 0.05 and |log2FC| >= 1).\n", sum(res_df$is_significant == "Yes", na.rm = TRUE)))


# --- 7. Save results ---
write.csv(res_df, "05_DE_analysis/deseq2_results.csv", row.names = TRUE)
print("DESeq2 results saved to 05_DE_analysis/deseq2_results.csv")


# --- 8. Exploratory Data Analysis (EDA) and Visualization ---
cat("Generating plots...\n")

# --- A. Dispersion Plot ---
pdf("05_DE_analysis/dispersion_plot.pdf")
plotDispEsts(dds)
mtext("Dispersion Estimates", side = 3, line = 1.5)
dev.off()
print("Dispersion plot saved to 05_DE_analysis/dispersion_plot.pdf")

# --- B. MA Plot ---
pdf("05_DE_analysis/MA_plot.pdf")
plotMA(res, ylim=c(-5,5), main="MA Plot (Notch2CKO vs WildType)")
abline(h=c(-1,1), col="blue")
dev.off()
print("MA plot saved to 05_DE_analysis/MA_plot.pdf")

# --- C. Variance Stabilizing Transformation (VST) ---
vsd <- vst(dds, blind = FALSE)
cat("VST transformation complete.\n")

# --- D. Principal Component Analysis (PCA) ---
pdf("05_DE_analysis/PCA_plot.pdf")
plotPCA(vsd, intgroup="condition") +
  ggtitle("PCA Plot of Samples") +
  theme_bw() +
  coord_fixed()
dev.off()
print("PCA plot saved to 05_DE_analysis/PCA_plot.pdf")

# --- E. Heatmap of top differentially expressed genes ---
top_genes_for_heatmap_df <- res_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  top_n(50, wt = abs(log2FoldChange))

if(nrow(top_genes_for_heatmap_df) > 0) {
  top_genes_ids <- rownames(top_genes_for_heatmap_df)
  mat <- assay(vsd)[top_genes_ids, ]
  mat <- mat - rowMeans(mat)

  df_annotation <- as.data.frame(colData[, c("condition")])
  colnames(df_annotation) <- "Condition"
  rownames(df_annotation) <- colnames(mat)

  pdf("05_DE_analysis/heatmap_top_DE_genes.pdf")
  pheatmap(mat, annotation_col = df_annotation,
           show_rownames = FALSE,
           fontsize_row = 8,
           fontsize_col = 10,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = "Heatmap of Top Differentially Expressed Genes",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
  dev.off()
  print("Heatmap of top DE genes saved to 05_DE_analysis/heatmap_top_DE_genes.pdf")
} else {
  print("No significant genes found (padj < 0.05 and |log2FC| >= 1) for heatmap plotting.")
}

# --- F. Counts Plot for Top Genes ---
top_5_genes <- head(res_df[order(res_df$padj), "gene"], 5)

if (length(top_5_genes) > 0) {
  pdf("05_DE_analysis/top_gene_counts_plots.pdf")
  for (gene_id in top_5_genes) {
    plotCounts(dds, gene=gene_id, intgroup="condition",
               main=paste0("Normalized Counts for: ", gene_id),
               xlab="Condition",
               col=ifelse(colData(dds)$condition == "WildType", "blue", "red"),
               pch=19, cex=1.2)
  }
  dev.off()
  print("Counts plots for top genes saved to 05_DE_analysis/top_gene_counts_plots.pdf")
} else {
  print("No top genes found to generate individual counts plots.")
}

# --- G. Heatmap of Sample-to-Sample Distances ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=" - ")
colnames(sampleDistMatrix) <- NULL

df_annotation_samples <- as.data.frame(colData(vsd)$condition)
colnames(df_annotation_samples) <- "Condition"
rownames(df_annotation_samples) <- colnames(assay(vsd))

colors_dist <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("05_DE_analysis/sample_distance_heatmap.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors_dist,
         annotation_col=df_annotation_samples,
         main="Sample-to-Sample Distances (VST Transformed)")
dev.off()
print("Sample-to-Sample Distance Heatmap saved to 05_DE_analysis/sample_distance_heatmap.pdf")

# --- H. Histogram of P-values ---
hist_data <- na.omit(res_df$pvalue)

pdf("05_DE_analysis/pvalue_histogram.pdf")
hist(hist_data, breaks=50, col="skyblue", border="white",
     main="Histogram of P-values",
     xlab="P-value",
     ylab="Frequency")
dev.off()
print("P-value Histogram saved to 05_DE_analysis/pvalue_histogram.pdf")

print("Differential Expression Analysis in R complete. Check 05_DE_analysis for results and plots.")
