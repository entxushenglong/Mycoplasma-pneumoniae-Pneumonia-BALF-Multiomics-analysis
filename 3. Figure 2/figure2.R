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


# Figure2 -------------------------------------------------------------------------
## Figure2A-D --------------------------------------------------------------------
# Taxa composition

# Phylum
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
p1 <- t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2")), 
                  others_color = "grey70", facet = "Group", xtext_keep = F, 
                  legend_text_italic = FALSE, ytitle_size = 17, strip_text = 14, 
                  clustering = F)

# Species
t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 10)
p2 <- t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2")), 
                  others_color = "grey70", facet = "Group", xtext_keep = T, 
                  legend_text_italic = FALSE, ytitle_size = 17, strip_text = 14, 
                  clustering = F, xtext_angle = 90)

# Taxa composition - By group mean
# Phylum
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
p3 <- t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2")), 
                  others_color = "grey70", xtext_size = 14, xtext_angle = 45)

# Species
t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 10, groupmean = "Group")
p4 <- t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2")), 
                  others_color = "grey70", xtext_size = 14, xtext_angle = 45)

phylum <- p3 + p1 + plot_layout(ncol = 2, widths = c(2, 10)); phylum
genus <- p4 + p2 + plot_layout(ncol = 2, widths = c(2, 10)); genus
combined_plot <- phylum / genus + plot_layout(heights = c(6, 6)); combined_plot
combined_plot <- combined_plot + theme(legend.text = element_text(size = 14))

ggsave("Figure2A-D.pdf", plot = combined_plot, height = 12, width = 24)

## Figure2 E-F ---------------------------------------------------------------
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
alpha_div <- t1$data_alpha
t1$cal_diff(method = "wilcox", p_adjust_method = "none")

shannon <- alpha_div[alpha_div$Measure == "Shannon", ]
simpson <- alpha_div[alpha_div$Measure == "Simpson", ]

# Define color mapping
fill_colors <- c("MucES-HS" = "#E6C5FC", "MucHS" = "#FEC8AF", "None" = "#DBE3C9")
color_colors <- c("MucES-HS" = "#7209B7", "MucHS" = "#B61624", "None" = "#5FB955")

p1 <- ggplot(shannon, aes(x = Group, y = Value, fill = Group)) +
  geom_violin() +  # Remove unsupported outlier.shape parameter
  geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) +
  labs(x = "", y = "Shannon") + theme_test() + theme(legend.position = "none", 
                                                     axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_rect(fill = NA)) +
  geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
              map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) +
  scale_fill_manual(values = fill_colors) + scale_color_manual(values = color_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))


p2 <- ggplot(simpson, aes(x = Group, y = Value, fill = Group)) +
  geom_violin() +  # Remove unsupported outlier.shape parameter
  geom_point(position = position_jitter(width = 0.2), size = 2, aes(color = Group)) +
  labs(x = "", y = "Simpson") + theme_test() + theme(legend.position = "none", 
                                                     axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_rect(fill = NA)) +
  geom_signif(comparisons = list(c("MucES-HS", "None"), c("MucES-HS", "MucHS"), c("MucHS", "None")),
              map_signif_level = FALSE, test = "wilcox.test", step_increase = 0.1) +
  scale_fill_manual(values = fill_colors) + scale_color_manual(values = color_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))


alpha <- p1 + p2 + plot_layout(ncol = 2, widths = c(3, 3)); alpha
ggsave("Figure2E-F.pdf", plot = alpha, height = 4, width = 8)

## Figure2G ----------------------------------------------------------------
t1 <- trans_venn$new(dataset_venn, ratio = "seqratio")
p <- t1$plot_venn(color_circle = c("#5FB955","#4125D0","#B61624"))
ggsave("Figure2G.pdf", plot = p, height = 6, width = 6)
rm(list = setdiff(ls(), c("dataset", "mytheme")))

## Figure2 H-I -------------------------------------------------------------
dataset$cal_betadiv(unifrac = F)

# Jaccard distance
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "jaccard")
t1$cal_ordination(method = "PCoA") # PCoA, PCA, DCA and NMDS are available
t1$cal_manova(manova_all = TRUE)
t1$res_manova

# Plot PCoA results and add annotation, using \n to create a new line, hjust = 1 for right alignment
pcoa_plot <- t1$plot_ordination(plot_color = "Group", 
                                plot_shape = "Group", 
                                color_values = c("#7209B7","#B61624", "#5FB955"),
                                plot_type = c("point", "ellipse")) +
  annotate("text", label = paste0("p = ", t1$res_manova$`Pr(>F)`[1]),  
           x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 4.5)

ggsave("Figure2H.pdf", plot = p1, height = 4, width = 5)

# Bray-Curtis distance
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# Calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# Return t1$res_group_distance
# Perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")

# View(t1$res_group_distance_diff)
write.csv(t1$res_group_distance_diff, "Bray−Curtis distance.csv")
# plot_group_order parameter can be used to adjust orders in x axis
p <- t1$plot_group_distance(boxplot_add = "mean", color_values = c("#4125D0","#B61624","#5FB955"))
ggsave("Figure2I.pdf", plot = p, height = 4, width = 4)

## Figure2J ----------------------------------------------------------------
# Import the dataset containing pre- and post-sampling data
dataset_double <- readRDS("dataset_double.rds")
# Filter rows where PatientID is not unique
dataset_double$sample_table <- dataset_double$sample_table %>% group_by(PatientID) %>% filter(n() > 1) %>% ungroup()
# Add SampleID column to row names
dataset_double$sample_table <- dataset_double$sample_table %>% as.data.frame() %>% column_to_rownames(var = "SampleID")
dataset_double$sample_table$SampleID <- row.names(dataset_double$sample_table)
dataset_double$tidy_dataset()

# Taxa composition
# Species
t1 <- trans_abund$new(dataset = dataset_double, taxrank = "Species", ntaxa = 10)
p2 <- t1$plot_bar(order_x = c("P60", "P93", "P77", "P78", "P12"), 
                  color_values = c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(8, "Dark2")),
                  others_color = "grey70", facet = "Time", xtext_keep = TRUE, 
                  legend_text_italic = FALSE, x_axis_name = "PatientID", xtext_angle = 45)

ggsave("Figure2J_update.pdf", plot = p2, height = 4, width = 6)