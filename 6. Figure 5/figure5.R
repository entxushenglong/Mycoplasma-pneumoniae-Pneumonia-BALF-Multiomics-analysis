# Figure5 -----------------------------------------------------------------
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

## Figure5 A--------------------------------------------------------------
library(apeglm)
pacman::p_load(tidyverse, microeco, magrittr)
library(ggtree)
library(dplyr)
library(ggplot2)
library("DESeq2")
library(ggpval)
library("clusterProfiler")

# Import gene count data
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

# Import column data (metadata)
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]

coldata <- coldata %>%
  mutate(Group = case_when(
    MucHS == "MucHS" & MucES == "MucES" ~ "MucES-HS",
    MucHS == "MucHS" & MucES != "MucES" ~ "MucHS",
    MucHS != "MucHS" & MucES == "MucES" ~ "MucES",
    TRUE ~ "None"))

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

# Transform data using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)


## Figure5A--------------------------------------------------------------

# PCA plot
p <- plotPCA(vsd, intgroup = c("Group"))
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Adjust colors and apply theme
p <- p + scale_color_manual(values = c("#7209B7", "#B61624", "#5FB955"))
p <- p + mytheme
ggsave("Figure5A.pdf", plot = p, height = 5, width = 6)

# Save results of MucES vs nonMucES comparison
res_condition_MucES_vs_nonMucES <- results(dds, contrast = c("MucES", "MucES", "nonMucES"))
resdataframe <- as.data.frame(res_condition_MucES_vs_nonMucES)
resdataframe$gene_name <- rownames(resdataframe)
write.csv(resdataframe, "MucES_vs_nonMucES.csv")

# Normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
merged_data_frame <- cbind(normalized_counts, resdataframe)

## Figure5B--------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)

# Extract normalized gene expression data
expression_data <- merged_data_frame[, colnames(normalized_counts)]

# Select top N most variable genes
topN <- 200  # For example, the top 200 most variable genes
variances <- apply(expression_data, 1, var)
top_genes <- head(order(variances, decreasing = TRUE), topN)
expression_data_top <- expression_data[top_genes, ]
annotation_col <- coldata[, c("Group"), drop = FALSE]

# Define group colors
group_colors <- c("MucES-HS" = "#E6C5FC", "MucHS" = "#FEC8AF", "None" = "#DBE3C9")

# Z-score normalization
expression_data_top <- t(scale(t(expression_data_top)))

# Create HeatmapAnnotation object for group information
ha <- HeatmapAnnotation(Group = coldata$Group, col = list(Group = group_colors))

# Plot heatmap
p <- Heatmap(expression_data_top,
             name = "Expression",
             top_annotation = ha,  # Add group annotation
             col = colorRamp2(c(-2, 0, 2), c("#357EB9", "white", "#DF6664")),  # Custom color gradient
             cluster_rows = T,  # Cluster rows
             cluster_columns = TRUE,  # Cluster columns
             show_row_names = FALSE,
             show_column_names = F)  # Hide column names

# Save heatmap to PDF
pdf("Figure5B_new.pdf", height = 5, width = 12)
draw(p)  # Render and save the heatmap
dev.off()

## Figure5C--------------------------------------------------------------
## KEGG Enrichment Analysis
# Volcano plot
library(ggVolcano)
merged_data_frame$gene_name <- rownames(merged_data_frame)
data <- add_regulate(merged_data_frame, log2FC_name = "log2FoldChange",
                     fdr_name = "padj", log2FC = 1, fdr = 0.05)
data <- data[!is.na(data$padj), ]

p <- ggvolcano(data, x = "log2FoldChange", y = "padj",
               label = "gene_name", output = FALSE,
               log2FC_cut = 1, pointSize = 0.1)
p
ggsave("Figure5C.pdf", p, width = 8, height = 5)


## Figure5D-------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

# KEGG Enrichment (upregulated genes)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange > 1)
write.csv(sig_dge, "MucES_vs_nonMucES_up.csv")

# Convert gene symbols to ENTREZID
genelist <- bitr(row.names(sig_dge), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)

# KEGG enrichment analysis
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr",
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 20)
p
ggsave("Figure5D_up.pdf", p, width = 8, height = 5)

## KEGG Enrichment (downregulated genes)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange < -1)
write.csv(sig_dge, "MucES_vs_nonMucES_down.csv")

# Convert gene symbols to ENTREZID
genelist <- bitr(row.names(sig_dge), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)

# KEGG enrichment analysis
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr",
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 20)
ggsave("Figure5D_down.pdf", p, width = 8, height = 5)

## Figure5E --------------------------------------------------------------
# Import count data
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

# Import column data
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]
# Ensure the order is consistent
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

resLFC <- lfcShrink(dds, coef = "MucHS_nonMucHS_vs_MucHS", type = "apeglm")
resLFC_dataframe <- as.data.frame(resLFC)

# Volcano plot for differential expression
merged_data_frame$gene_name <- rownames(merged_data_frame)
data <- add_regulate(merged_data_frame, log2FC_name = "log2FoldChange",
                     fdr_name = "padj", log2FC = 1, fdr = 0.05)
data <- data[!is.na(data$padj), ]

p <- ggvolcano(data, x = "log2FoldChange", y = "padj",
               label = "gene_name", output = FALSE,
               log2FC_cut = 1, pointSize = 0.1)
p
ggsave("Figure5E.pdf", p, width = 8, height = 5)


## Figure5F --------------------------------------------------------------
# KEGG enrichment (upregulated genes)
sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange > 1)
write.csv(sig_dge, "MucHS_vs_nonMucHS_up.csv")

# Convert gene symbols to ENTREZID
genelist <- bitr(row.names(sig_dge), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)

# KEGG enrichment analysis
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr",
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 20)
ggsave("Figure5F_up.pdf", p, width = 8, height = 5)


sig_dge <- subset(resdataframe, padj < 0.05)
sig_dge <- subset(sig_dge, log2FoldChange < -1)
write.csv(sig_dge,"MucHS_vs_nonMucHS_down.csv")
genelist <- bitr(row.names(sig_dge), fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)               


ekegg <- enrichKEGG(gene = genelist, organism = 'hsa', 
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", 
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 20)

ggsave("Figure5F_down.pdf", p, width = 8, height = 5)




## Figure5 G-H --------------------------------------------------------------

# Import gene count data
cts <- read.csv("genes.readcount.csv")
rownames(cts) <- cts$geneID
cts <- cts[,-1]

# Import column data
coldata <- read.csv("coldata.csv")
# Only keep rows with RNAseq data
coldata <- coldata[coldata$SampleID %in% colnames(cts), ]
# Ensure the order is consistent
cts <- cts[, coldata$SampleID]

# Create a new DESeqDataSet object
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

# Results for comparison between Severe and Mild stages
res_condition_Severe_vs_Mild <- results(dds, contrast = c("Stage", "Severe", "Mild"))
resdataframe <- as.data.frame(res_condition_Severe_vs_Mild)
resdataframe$gene_name <- rownames(resdataframe)

# Select genes with baseMean greater than 200
gene200 <- rownames(resdataframe[resdataframe$baseMean > 200, ])
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% gene200, ]
normalized_counts <- t(as.data.frame(normalized_counts))
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts$SampleID <- row.names(normalized_counts)

# Load sample table from dataset
dataset <- readRDS("dataset.rds")
sample_table <- dataset$sample_table

# Merge normalized counts with sample table
normalized_counts <- merge(normalized_counts, sample_table[, c("SampleID", "MP_load")], by = "SampleID")

# Initialize result data frame
result_df <- data.frame(Gene = colnames(normalized_counts)[2:10890], 
                        Spearman_r = numeric(length = 10889), 
                        p_value = numeric(length = 10889))

# Perform Spearman correlation analysis for each column
for (i in 2:10890) {
  gene_name <- colnames(normalized_counts)[i]
  correlation_test <- cor.test(normalized_counts[, i], normalized_counts[, 10891], method = "spearman")
  
  result_df$Spearman_r[i - 1] <- correlation_test$estimate
  result_df$p_value[i - 1] <- correlation_test$p.value
}

# Filter genes with Spearman_r greater than 0.4 and extract Gene names
positive_gene <- result_df$Gene[result_df$Spearman_r > 0.4]
genelist <- bitr(positive_gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)

# KEGG enrichment analysis for positive correlations
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa', 
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", 
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 10)
ggsave("Figure5G.pdf", p, height = 5, width = 8)

# Filter genes with Spearman_r less than -0.4 and extract Gene names
negative_gene <- result_df$Gene[result_df$Spearman_r < -0.4]
genelist <- bitr(negative_gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)

# KEGG enrichment analysis for negative correlations
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa', 
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", 
                    minGSSize = 1)

p <- dotplot(ekegg, showCategory = 10)
ggsave("Figure5H.pdf", p, height = 5, width = 8)