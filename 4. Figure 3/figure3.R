# import packages -------------------------------------------------------------------------
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

# import dataset file which contains microbiome data
dataset <- readRDS("dataset.rds")

# Figure3 -----------------------------------------------------------------

## Figure3 A-C --------------------------------------------------------------
sample_table <- dataset$sample_table

# Define color mapping
fill_colors <- c("MucES-HS" = "#E6C5FC", "MucHS" = "#FEC8AF", "None" = "#DBE3C9")
color_colors <- c("MucES-HS" = "#7209B7", "MucHS" = "#B61624", "None" = "#5FB955")

# Plot for Microbiome Load
p1 <- ggplot(sample_table, aes(x = Group, y = Microbiome_load, fill = Group)) +
  geom_violin() +  # Removed unsupported outlier.shape parameter
  geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) +
  labs(x = "", y = "Microbiome Load") +
  theme_test() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = NA)) + 
  geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
              map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) + 
  scale_fill_manual(values = fill_colors) + 
  scale_color_manual(values = color_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
p1

# Plot for MP Load
p2 <- ggplot(sample_table, aes(x = Group, y = MP_load, fill = Group)) + 
  geom_violin() +  # Removed unsupported outlier.shape parameter
  geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) + 
  labs(x = "", y = "MP Load") + 
  theme_test() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = NA)) + 
  geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
              map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) + 
  scale_fill_manual(values = fill_colors) + 
  scale_color_manual(values = color_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
p2

# Plot for MP Percentage
p3 <- ggplot(sample_table, aes(x = Group, y = MP_percentage, fill = Group)) + 
  geom_violin() +  # Removed unsupported outlier.shape parameter
  geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) + 
  labs(x = "", y = "MP Percentage") + 
  theme_test() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = NA)) + 
  geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
              map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) + 
  scale_fill_manual(values = fill_colors) + 
  scale_color_manual(values = color_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
p3

# Combine the plots
p <- p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(5, 5, 5)); p
ggsave("Figure3A-C.pdf", plot = p, height = 4, width = 12)

## Figure3 D-E -------------------------------------------------------------
# KW difference analysis
t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group", taxa_level = "all", filter_thres = 0.0001, p_adjust_method = "none")

p <- t1$plot_diff_abund(use_number = 1:20, add_sig = T, add_sig_label = "Significance", coord_flip = F,  
                        group_order = c("MucES-HS", "MucHS", "Normal"), color_values = c("#7209B7","#B61624","#5FB955"))
ggsave("Figure3D.pdf", plot = p, height = 5, width = 10)

# LefSe analysis
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.05, p_adjust_method = "none", lefse_subgroup = NULL)
lefse_res <- t1$res_diff
p <- t1$plot_diff_bar(use_number = 1:10, width = 0.8, group_order = c("MucES-HS", "MucHS", "Normal"), color_values = c("#7209B7","#5FB955"))
ggsave("Figure3E.pdf", plot = p, width = 6, height = 5)


## Figure3F ----------------------------------------------------------
abund_sig <- dataset$sample_table[, 63:65]
sample_table_cor <- dataset$sample_table[, c(8:45, 47:61)]
sample_table_cor <- dataset$sample_table[, c("K", "Na", "Cal", "P", "totalPro", "Albumin", "Globulin", "Urea", "TotalCholesterol", "Urate", "ALP", "AST",
                                             "ALT", "ADAR", "CKMB", "LDH", "LDL.C", "FFA", "CHE", "PAB", "PLT", "Neutrocyte.", "Lymphocyte.", "Eosinophil.")]

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

corr_colors <- colorRampPalette(c("#327EBA", "white", "#DF6664"))(200)
pdf("Figure3F.pdf", width = 15, height = 8)
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
         tl.cex = 1.5, 
         tl.col = "black", 
         tl.srt = 45)
dev.off()

## Figure3 G-J ----------------------------------------------------------------
library(psych)
# Define column names
columns <- c("P", "Albumin", "ALT", "Neutrocyte.")
# Initialize an empty list to store ggplot objects
plot_list <- list()
# Loop through each column and generate ggplot objects
for (col in columns) {
  
  corr1 <- corr.test(dataset$sample_table[[col]], dataset$sample_table$MP_load, method = "spearman", adjust = "none")
  
  p <- ggplot(dataset$sample_table, aes(x = .data[[col]], y = MP_load)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = -Inf, y = Inf, label = sprintf("spearman r = %.3f", corr1$r), hjust = 0, vjust = 3, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("p = %.3f", corr1$p), hjust = 0, vjust = 1.5, size = 5)
  
  # Add generated ggplot object to the list
  plot_list[[col]] <- p
}


pdf("Figure3G-J.pdf", height = 3, width = 13)
do.call("grid.arrange", c(plot_list, ncol = 4, nrow = 1))
dev.off()