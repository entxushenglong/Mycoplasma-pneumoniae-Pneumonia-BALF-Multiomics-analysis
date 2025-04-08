# Import packages ---------------------------------------------------------------
library(ggtree)
library(dplyr)
library(ggplot2)
library(file2meco)
library(ggpval)
library(RColorBrewer)
library(corrplot)
library(export)
library(stringr)
library(microeco)
library(file2meco)
library(microeco)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(ggdendro)
library(ggnested)
library(mecodev)
library(randomForest)
library(patchwork)
library(gridExtra)
library(ggsignif)
library(apeglm)
pacman::p_load(tidyverse, microeco, magrittr)
library(ggtree)
library(dplyr)
library(ggplot2)
library("DESeq2")
library(ggpval)
library("clusterProfiler")

# MucES Differentially Expressed Genes
#import cts
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

#import coldata
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]
# Ensure both datasets are in the same order
cts <- cts[, coldata$SampleID]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ MucES)

smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- relevel(dds$MucES, ref = "nonMucES")

dds <- DESeq(dds)
resultsNames(dds)

# Shrink log fold changes using apeglm
resLFC <- lfcShrink(dds, coef = "MucES_nonMucES_vs_MucES", type = "apeglm")
resLFC
resLFC_dataframe <- as.data.frame(resLFC)

# Results for MucES vs nonMucES comparison
res_condition_MucES_vs_nonMucES <- results(dds, contrast = c("MucES", "MucES", "nonMucES"))
resdataframe <- as.data.frame(res_condition_MucES_vs_nonMucES)
resdataframe$gene_name <- rownames(resdataframe)

normalized_counts <- counts(dds, normalized = TRUE)
merged_data_frame <- cbind(normalized_counts, resdataframe)

# KEGG enrichment (upregulated)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange > 1)
MucES_up <- sig_dge
MucES_up <- MucES_up$gene_name

# KEGG enrichment (downregulated)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange < -1)
MucES_down <- sig_dge
MucES_down <- MucES_down$gene_name

# MucHS Differentially Expressed Genes
#import cts
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

#import coldata
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]
# Ensure both datasets are in the same order
cts <- cts[, coldata$SampleID]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ MucHS)

smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- relevel(dds$MucHS, ref = "nonMucHS")

dds <- DESeq(dds)
resultsNames(dds)

# Shrink log fold changes using apeglm
resLFC <- lfcShrink(dds, coef = "MucHS_nonMucHS_vs_MucHS", type = "apeglm")
resLFC
resLFC_dataframe <- as.data.frame(resLFC)

# Results for MucHS vs nonMucHS comparison
res_condition_MucHS_vs_nonMucHS <- results(dds, contrast = c("MucHS", "MucHS", "nonMucHS"))
resdataframe <- as.data.frame(res_condition_MucHS_vs_nonMucHS)
resdataframe$gene_name <- rownames(resdataframe)

normalized_counts <- counts(dds, normalized = TRUE)
merged_data_frame <- cbind(normalized_counts, resdataframe)

# KEGG enrichment (upregulated)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange > 1)
MucHS_up <- sig_dge
MucHS_up <- MucHS_up$gene_name

# KEGG enrichment (downregulated)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange < -1)
MucHS_down <- sig_dge
MucHS_down <- MucHS_down$gene_name

# Identify common upregulated and downregulated genes
up_genes <- intersect(MucES_up, MucHS_up)
down_genes <- intersect(MucES_down, MucHS_down)

# Merge and filter
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

#import coldata
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]
# Ensure both datasets are in the same order
cts <- cts[, coldata$SampleID]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Stage)

smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- relevel(dds$Stage, ref = "Mild")

dds <- DESeq(dds)
resultsNames(dds)

# Shrink log fold changes using apeglm
resLFC <- lfcShrink(dds, coef = "Stage_Severe_vs_Mild", type = "apeglm")
resLFC

# Results for Severe vs Mild Stage comparison
res_condition_Severe_vs_Mild <- results(dds, contrast = c("Stage", "Severe", "Mild"))
resdataframe <- as.data.frame(res_condition_Severe_vs_Mild)
resdataframe$gene_name <- rownames(resdataframe)

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- t(as.data.frame(normalized_counts))
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts$SampleID <- row.names(normalized_counts)

dataset <- readRDS("dataset.rds")
sample_table <- dataset$sample_table

# Merge normalized counts with sample table
normalized_counts <- merge(normalized_counts, sample_table[, c("SampleID", "MP_percentage", "MP_load", "Microbiome_load", "K", "Na", "Cal", "P", 
                                                               "totalPro", "Albumin", "Globulin", "Urea", "TotalCholesterol", "Urate", 
                                                               "ALP", "AST", "ALT", "ADAR", "CKMB", "LDH", "LDL.C", "FFA", "CHE", 
                                                               "PAB", "PLT", "Neutrocyte.", "Lymphocyte.", "Eosinophil.", "Group", 
                                                               "MucHS", "MucES")], by = "SampleID")


## Figure6C -------------------------------------------------------
library(ggpubr)
library(ggplot2)
# Define the list of genes to plot
genes <- c("CCL2", "CCL7", "CCL8", "BNIP3", "LOC105376219", "MAP3K7CL", "LY6E")

# Define color mapping
fill_colors <- c("MucES-HS" = "#E6C5FC", "MucHS" = "#FEC8AF", "None" = "#DBE3C9")
color_colors <- c("MucES-HS" = "#7209B7", "MucHS" = "#B61624", "None" = "#5FB955")

# Create an empty list to store the plots
plot_list <- list()

# Loop to generate a plot for each gene
for (gene in genes) {
  p <- ggplot(normalized_counts, aes(x = Group, y = .data[[gene]], fill = Group)) +
    geom_violin(outlier.shape = NA) +  # Remove outlier points
    geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) + 
    labs(x = "", y = gene) + 
    theme_test() + 
    theme(legend.position = "none",  # Hide legend for each individual plot
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          strip.background = element_rect(fill = NA)) + 
    geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
                map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) + 
    scale_fill_manual(values = fill_colors) + 
    scale_color_manual(values = color_colors) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
  
  # Add the plot to the list
  plot_list[[gene]] <- p
}

# Keep the legend of the first plot
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "top")

# Combine all plots using patchwork package and share one legend
DEG <- wrap_plots(plot_list, ncol = 7) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")
DEG

ggsave("Figure6C.pdf", plot = DEG, width = 20, height = 5)

## Figure6D ---------------------------------------------------------------
abund_sig <- normalized_counts[, c("BNIP3", "CCL2", "CCL7", "CCL8", "LOC105376219", "MAP3K7CL", "LY6E")]
sample_table_cor <- normalized_counts[, c("Microbiome_load", "MP_percentage", "MP_load", "K", "Na", "Cal", "P",
                                          "totalPro", "Albumin", "Globulin", "Urea", "TotalCholesterol", "Urate", 
                                          "ALP", "AST", "ALT", "CKMB", "LDH", "LDL.C", "FFA", "CHE", "PAB", "PLT", 
                                          "Neutrocyte.", "Lymphocyte.", "Eosinophil.")]

# Initialize a data frame to store correlation results
# Column names from sample_table_cor and row names from abund_sig
corr <- as.data.frame(matrix(nrow = ncol(abund_sig), ncol = ncol(sample_table_cor)))
rownames(corr) <- colnames(abund_sig)
colnames(corr) <- colnames(sample_table_cor)

corr_p <- matrix(nrow = ncol(abund_sig), ncol = ncol(sample_table_cor))
rownames(corr_p) <- colnames(abund_sig)
colnames(corr_p) <- colnames(sample_table_cor)

# Compute correlation and p-values
for(i in 1:ncol(abund_sig)) {
  for(j in 1:ncol(sample_table_cor)) {
    test_result <- cor.test(abund_sig[, i], sample_table_cor[, j], method = "spearman")
    corr[i, j] <- test_result$estimate  # Store correlation coefficient
    corr_p[i, j] <- test_result$p.value  # Store p-value
  }
}

# Find rows in corr_p with p-values less than 0.05
significant_rows <- apply(corr_p, 1, function(row) any(row < 0.05))

# Keep rows in corr_p with significant p-values
corr_p_significant <- corr_p[significant_rows, ]

# Filter the corresponding rows in corr
corr_significant <- corr[rownames(corr_p_significant), ]

corr_colors <- colorRampPalette(c("#327EBA", "white", "#DF6664"))(200)

pdf("Figure6D.pdf", width = 13, height = 8)

corrplot(as.matrix(corr),
         method = "circle",
         type = "full",
         col = corr_colors,  # Specify color palette
         p.mat = as.matrix(corr_p),
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.cex = 2,
         addrect = 4,
         cl.pos = "r", 
         tl.cex = 1, 
         tl.col = "black", 
         tl.srt = 45)
dev.off()

## Figure6 E-F ------------------------------------------------------------------
library(randomForest)
library(pROC)
library(caret)

## MucES 

# Convert MucES to factor
normalized_counts$MucES <- as.factor(normalized_counts$MucES)
# Set random seed to ensure reproducibility
set.seed(9)
# Stratified sampling based on MucES column, divide into training and test sets
partition <- createDataPartition(normalized_counts$MucES, p = 0.6, list = FALSE)

# Create training and test sets
trainSet <- normalized_counts[partition, ]
testSet <- normalized_counts[-partition, ]

# DEGs
RF_model_DEG <- randomForest(MucES ~ CCL2 + CCL7 + CCL8 + LY6E, data = trainSet)
DEG_res <- predict(RF_model_DEG, newdata = testSet, type = "prob")[, 2]

# Clinical model using columns from 23310 to 23333
RF_model_clin <- randomForest(MucES ~ ., 
                              data = trainSet[, c(23310:23333, which(names(trainSet) == "MucES"))], 
                              na.action = na.roughfix)
clin_res <- predict(RF_model_clin, newdata = testSet, type = "prob")[, 2]

# Random model using random genes
RF_model_rand <- randomForest(MucES ~ A1BG + A1CF + A2M, 
                              data = trainSet, 
                              na.action = na.roughfix)
random_res <- predict(RF_model_rand, newdata = testSet, type = "prob")[, 2]

# Generate ROC curve
roc_curve_DEG <- roc(testSet$MucES, predictor = DEG_res, levels = c("non-MucES", "MucES"))
roc_curve_clin <- roc(response = testSet$MucES, predictor = clin_res, levels = c("non-MucES", "MucES"))
roc_curve_rand <- roc(response = testSet$MucES, predictor = random_res, levels = c("non-MucES", "MucES"))

# Calculate AUC
AUC_DEG <- auc(roc_curve_DEG)
AUC_clin <- auc(roc_curve_clin)
AUC_rand <- auc(roc_curve_rand)

# Save plot to PDF
pdf("Figure6E.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(roc_curve_DEG, col = "#C5070A", main = "Prediction Model for MucES")
lines(roc_curve_clin, col = "#1E90FF")
lines(roc_curve_rand, col = "#5811BB")
# Add a legend with AUC values in the bottom right corner
legend("bottomright", 
       legend = c(paste("CCL Model AUC:", round(AUC_DEG, 2)), 
                  paste("Clinical Model AUC:", round(AUC_clin, 2)), 
                  paste("Random Gene Model AUC:", round(AUC_rand, 2))),
       col = c("#C5070A", "#1E90FF", "#5811BB"), 
       lwd = 3, 
       cex = 0.8, # Adjust legend text size
       bty = "o") # Set the legend box to a rectangle

# Close the graphics device and save the file
dev.off()

# MucHS
# Convert MucHS to factor
normalized_counts$MucHS <- as.factor(normalized_counts$MucHS)

# Set random seed to ensure reproducibility
set.seed(15)
# Stratified sampling based on MucHS column, divide into training and test sets
partition <- createDataPartition(normalized_counts$MucHS, p = 0.6, list = FALSE)
# Create training and test sets
trainSet <- normalized_counts[partition, ]
testSet <- normalized_counts[-partition, ]

# DEGs
RF_model_DEG <- randomForest(MucHS ~ CCL2 + CCL7 + CCL8 + LY6E, data = trainSet)
DEG_res <- predict(RF_model_DEG, newdata = testSet, type = "prob")[, 2]

# Clinical model using columns from 23310 to 23333
RF_model_clin <- randomForest(MucHS ~ ., 
                              data = trainSet[, c(23310:23333, which(names(trainSet) == "MucHS"))], 
                              na.action = na.roughfix)
clin_res <- predict(RF_model_clin, newdata = testSet, type = "prob")[, 2]

# Random model using random genes
RF_model_rand <- randomForest(MucHS ~ A1BG + A1CF + A2M, 
                              data = trainSet, 
                              na.action = na.roughfix)
random_res <- predict(RF_model_rand, newdata = testSet, type = "prob")[, 2]

# Generate ROC curve
roc_curve_DEG <- roc(testSet$MucHS, predictor = DEG_res, levels = c("non-MucHS", "MucHS"))
roc_curve_clin <- roc(response = testSet$MucHS, predictor = clin_res, levels = c("non-MucHS", "MucHS"))
roc_curve_rand <- roc(response = testSet$MucHS, predictor = random_res, levels = c("non-MucHS", "MucHS"))

# Calculate AUC
AUC_DEG <- auc(roc_curve_DEG)
AUC_clin <- auc(roc_curve_clin)
AUC_rand <- auc(roc_curve_rand)

# Save plot to PDF
pdf("Figure6F.pdf", width = 5, height = 5)
# Plot the first ROC curve
plot(roc_curve_DEG, col = "#C5070A", main = "Prediction Model for MucHS")
lines(roc_curve_clin, col = "#1E90FF")
lines(roc_curve_rand, col = "#5811BB")
# Add a legend with AUC values in the bottom right corner
legend("bottomright", 
       legend = c(paste("CCL Model AUC:", round(AUC_DEG, 2)), 
                  paste("Clinical Model AUC:", round(AUC_clin, 2)), 
                  paste("Random Gene Model AUC:", round(AUC_rand, 2))),
       col = c("#C5070A", "#1E90FF", "#5811BB"), 
       lwd = 3, 
       cex = 0.8, # Adjust legend text size
       bty = "o") # Set the legend box to a rectangle
# Close the graphics device and save the file
dev.off()