# Clear working space
rm(list=ls())
# Clear console
cat("\014")
# Clear all open graphics devices
while (dev.cur() > 1) dev.off()

# Imports
library(data.table)
library(dplyr)
library(ggplot2)

################################################################################
## Principal Component Analysis (PCA)
################################################################################

# Loading relevant RNAseq object
rnaseq_obj      <- readRDS("./data/RNAseq_processing/complete_data_obj/complete_RNAseq_data_object.RDS")
cpm_filt_counts <- rnaseq_obj$count_filt_cpm_df
metadata        <- rnaseq_obj$meta_df

# Extract log2CPM for selected genes, samples as rows, genes as columns
tmp.log2cpm.df           <- log2(cpm_filt_counts + 1) 
tmp.log2cpm.df           <- as.data.frame(t(tmp.log2cpm.df))

# PCA calculation
data.pca <- prcomp(tmp.log2cpm.df, center = TRUE, scale. = TRUE)

# Extract variance explained
pct_var  <- summary(data.pca)$importance[2, ] * 100

# Build plot data.frame
dfPlot   <- as.data.frame(data.pca$x) %>%
  merge(metadata, by = "row.names") %>%
  tibble::column_to_rownames(var = "Row.names")

dfPlot$condition <- as.factor(dfPlot$condition)

# Plotting
pl <- ggplot(dfPlot, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 3, aes(shape = condition)) +
  stat_ellipse(geom      = "polygon",
               level     = 0.8,
               alpha     = 0,
               linewidth = 0.7) +
  scale_color_manual(values = c("CD_PBS"  = "#d4dbff",
                                "HFD_PBS" = "#fdc812",
                                "HFD_CS"  = "#fd7812")) +
  scale_fill_manual(values  = c("CD_PBS"  = "#d4dbff",
                                "HFD_PBS" = "#fdc812",
                                "HFD_CS"  = "#fd7812")) +
  labs(x = sprintf("PC1 (%.1f%%)", pct_var[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var[2])) +
  theme_classic()

# Saving
save_dir <- "./plots/PCA/"

if (!dir.exists(save_dir)){
  dir.create(save_dir, recursive = T)
}

ggsave(paste0(save_dir, "jejunum_rnaseq_pca.pdf"), plot = pl, width = 9, height = 9, useDingbats = F)









# Clear working space
rm(list=ls())
cat("\014")
while (dev.cur() > 1) dev.off()

# Imports
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

################################################################################
## Data loading
################################################################################
rnaseq_obj      <- readRDS("./data/RNAseq_processing/complete_data_obj/complete_RNAseq_data_object.RDS")
cpm_filt_counts <- rnaseq_obj$count_filt_cpm_df        # genes x samples
metadata        <- rnaseq_obj$meta_df

# log2CPM, transpose ONCE -> samples as rows, genes as columns
tmp.log2cpm.df <- log2(cpm_filt_counts + 1)
tmp.log2cpm.df <- as.data.frame(t(tmp.log2cpm.df))     # <-- exactly one transpose

# Align metadata to the expression matrix
metadata           <- metadata[rownames(tmp.log2cpm.df), , drop = FALSE]
metadata$condition <- as.factor(metadata$condition)

################################################################################
## DIAGNOSTIC — must pass before anything else runs
################################################################################
cat("Matrix dims (should be n_samples x ~18000 genes):",
    nrow(tmp.log2cpm.df), "x", ncol(tmp.log2cpm.df), "\n")
cat("Condition levels found:", paste(levels(metadata$condition), collapse = ", "), "\n")
print(table(metadata$condition))

# Hard stops so you can't silently proceed on a broken matrix
if (ncol(tmp.log2cpm.df) < 1000)
  stop("Fewer than 1000 gene columns -> matrix is transposed the wrong way. Check the t() step.")
if (nrow(metadata) != nrow(tmp.log2cpm.df))
  stop("metadata rows do not match expression rows.")

# Map the three biological groups to whatever the actual level names are, by
# pattern — so the script works regardless of exact spelling (HFD_CS vs CS_HFD etc.)
lv <- levels(metadata$condition)
grp_CD     <- lv[grepl("CD",  lv, ignore.case = TRUE)]                       # chow
grp_HFD    <- lv[grepl("HFD", lv, ignore.case = TRUE) & !grepl("CS|SCIND", lv, ignore.case = TRUE)]  # HFD alone
grp_CS_HFD <- lv[grepl("HFD", lv, ignore.case = TRUE) &  grepl("CS|SCIND", lv, ignore.case = TRUE)]  # HFD + C. scindens

cat("\nResolved groups:\n")
cat("  CD      =", grp_CD,     "\n")
cat("  HFD     =", grp_HFD,    "\n")
cat("  CS_HFD  =", grp_CS_HFD, "\n")

if (length(grp_CD)!=1 || length(grp_HFD)!=1 || length(grp_CS_HFD)!=1)
  stop("Could not uniquely resolve the three groups from the condition levels above. Set grp_CD / grp_HFD / grp_CS_HFD manually.")

# Colours keyed to the ACTUAL level names
cond_colors <- setNames(c("#d4dbff", "#fdc812", "#fd7812"),
                        c(grp_CD, grp_HFD, grp_CS_HFD))

save_dir <- "./plots/PCA/"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

################################################################################
## PCA (all genes)
################################################################################
data.pca <- prcomp(tmp.log2cpm.df, center = TRUE, scale. = TRUE)
pct_var  <- summary(data.pca)$importance[2, ] * 100

dfPlot <- as.data.frame(data.pca$x) %>%
  merge(metadata, by = "row.names") %>%
  tibble::column_to_rownames(var = "Row.names")
dfPlot$condition <- as.factor(dfPlot$condition)
dfPlot$sample_id <- rownames(dfPlot)

pl <- ggplot(dfPlot, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 3, aes(shape = condition)) +
  stat_ellipse(geom = "polygon", level = 0.8, alpha = 0, linewidth = 0.7) +
  scale_color_manual(values = cond_colors) +
  scale_fill_manual(values  = cond_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", pct_var[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var[2])) +
  theme_classic()
ggsave(paste0(save_dir, "jejunum_rnaseq_pca.pdf"), plot = pl, width = 9, height = 9, useDingbats = FALSE)

pl_lab <- pl + geom_text_repel(aes(label = sample_id), size = 2.5,
                               max.overlaps = Inf, show.legend = FALSE)
ggsave(paste0(save_dir, "jejunum_rnaseq_pca_labelled.pdf"), plot = pl_lab, width = 9, height = 9, useDingbats = FALSE)

################################################################################
## Sample-to-sample distance (top-variable genes)
################################################################################
n_top_genes <- 500
gene_var  <- apply(tmp.log2cpm.df, 2, var)             # variance per GENE (columns)
sel_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(n_top_genes, length(gene_var))]
dist_mat_input <- tmp.log2cpm.df[, sel_genes, drop = FALSE]
cat(sprintf("\nDistance analysis on top %d variable genes.\n", length(sel_genes)))

sample_dist     <- dist(dist_mat_input, method = "euclidean")
sample_dist_mat <- as.matrix(sample_dist)

grp <- metadata$condition
D   <- sample_dist_mat
mean_between <- function(a, b) mean(D[grp == a, grp == b])
mean_within  <- function(a) { s <- D[grp == a, grp == a]; mean(s[lower.tri(s)]) }

cat("\n--- Mean sample-to-sample Euclidean distances (top-variable genes) ---\n")
cat(sprintf("CS_HFD vs CD   : %.2f\n", mean_between(grp_CS_HFD, grp_CD)))
cat(sprintf("CS_HFD vs HFD  : %.2f\n", mean_between(grp_CS_HFD, grp_HFD)))
cat(sprintf("CD     vs HFD  : %.2f\n", mean_between(grp_CD,     grp_HFD)))
cat(sprintf("within CD      : %.2f\n", mean_within(grp_CD)))
cat(sprintf("within CS_HFD  : %.2f\n", mean_within(grp_CS_HFD)))
cat(sprintf("within HFD     : %.2f\n", mean_within(grp_HFD)))

d_cs_cd  <- mean_between(grp_CS_HFD, grp_CD)
d_cs_hfd <- mean_between(grp_CS_HFD, grp_HFD)
ratio    <- d_cs_cd / d_cs_hfd
cat(sprintf("\nCS_HFD-CD / CS_HFD-HFD ratio = %.2f  ", ratio))
if (is.na(ratio)) {
  cat("=> still NA: check the resolved group names above.\n")
} else if (ratio < 0.85) {
  cat("=> CS_HFD closer to CD.\n")
} else if (ratio > 1.15) {
  cat("=> CS_HFD closer to HFD.\n")
} else {
  cat("=> CS_HFD intermediate between CD and HFD.\n")
}

################################################################################
## Heatmaps + dendrogram
################################################################################
annot_df  <- data.frame(condition = metadata$condition, row.names = rownames(metadata))
annot_col <- list(condition = cond_colors[levels(metadata$condition)])

pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         clustering_method        = "complete",
         annotation_row = annot_df, annotation_col = annot_df,
         annotation_colors = annot_col,
         color = colorRampPalette(c("#08306b","#4292c6","#deebf7","#ffffff"))(100),
         main = sprintf("Sample-to-sample Euclidean distance (top %d variable genes)", n_top_genes),
         filename = paste0(save_dir, "jejunum_rnaseq_sample_distance_heatmap.pdf"),
         width = 9, height = 9)

sample_cor      <- cor(t(dist_mat_input), method = "pearson")
sample_cor_dist <- as.dist(1 - sample_cor)
pheatmap(sample_cor,
         clustering_distance_rows = sample_cor_dist,
         clustering_distance_cols = sample_cor_dist,
         clustering_method        = "complete",
         annotation_row = annot_df, annotation_col = annot_df,
         annotation_colors = annot_col,
         color = colorRampPalette(c("#ffffff","#deebf7","#4292c6","#08306b"))(100),
         main = sprintf("Sample-to-sample Pearson correlation (top %d variable genes)", n_top_genes),
         filename = paste0(save_dir, "jejunum_rnaseq_sample_correlation_heatmap.pdf"),
         width = 9, height = 9)

pdf(paste0(save_dir, "jejunum_rnaseq_sample_distance_dendrogram.pdf"), width = 9, height = 6, useDingbats = FALSE)
plot(hclust(sample_dist, method = "complete"),
     main = sprintf("Hierarchical clustering (Euclidean, top %d variable genes)", n_top_genes),
     xlab = "", sub = "")
dev.off()







