# Clear working space
rm(list=ls())
# Clear console
cat("\014")
# Clear all open graphics devices
while (dev.cur() > 1) dev.off()

# Imports
library(limma)
library(edgeR)
library(ggplot2)
library(data.table)

################################################################################
## Differential Expression Analysis (DEA)
################################################################################

# Load RNAseq object
rnaseq_obj <- readRDS("./data/RNAseq_processing/complete_data_obj/complete_RNAseq_data_object.RDS")
counts.df  <- rnaseq_obj$count_filt_df
meta.df    <- rnaseq_obj$meta_df
stopifnot(identical(rownames(meta.df), colnames(counts.df)))

# Run DEA
y <- edgeR::DGEList(counts = counts.df, samples = meta.df)
y <- edgeR::calcNormFactors(y)

cond_order                <- c("HFD_CS", "HFD_PBS", "CD_PBS")
stopifnot(all(cond_order %in% unique(meta.df$condition)))
cond_comb_df              <- as.data.frame(t(combn(cond_order, 2)))
colnames(cond_comb_df)    <- c("cond1", "cond2")
cond_comb_df$contrast_id  <- paste0(cond_comb_df$cond1, "__vs__", cond_comb_df$cond2)
cond_comb_df$contrast_def <- paste0("(", cond_comb_df$cond1, ")-(", cond_comb_df$cond2, ")")
all_contrasts             <- cond_comb_df$contrast_def
names(all_contrasts)      <- cond_comb_df$contrast_id
all_contrasts

designFormula  <- as.formula("~0+condition")
mm             <- model.matrix(designFormula, data = y$samples)
colnames(mm)   <- gsub("^condition", "", colnames(mm))

yv        <- voom(y, mm, plot = TRUE)
contrasts <- lapply(all_contrasts, function(k) {
  eval(parse(text = paste0("out <- makeContrasts(", unname(k), ", levels = mm)")))
  colnames(out) <- names(k)
  out
})

fitModel   <- limma::lmFit(yv, design = mm)
top_tables <- lapply(names(contrasts), function(k) {
  fit_obj.k               <- limma::eBayes(limma::contrasts.fit(fitModel, contrasts = contrasts[[k]]))
  table.k                 <- limma::topTable(fit_obj.k, sort.by = "none", n = Inf)
  table.k$ensembl_gene_id <- rownames(table.k)
  table.k                 <- merge(table.k, rnaseq_obj$gene_conv_df, by.x = "ensembl_gene_id", by.y = "gene_id")
  table.k$contrast_id     <- k
  table.k                 <- table.k[!is.na(table.k$logFC), ]
  table.k                 <- table.k[order(table.k$adj.P.Val, -1 * abs(table.k$logFC)), ]
  table.k
})
names(top_tables) <- names(contrasts)

################################################################################
## Saving
################################################################################

saveDir <- "./data/DEA/"

if (!dir.exists(saveDir)){
  dir.create(saveDir, recursive = TRUE)
}

all_top_tables.df <- Reduce(function(...) merge(..., by = c("ensembl_gene_id", "gene_name", "gene_biotype"), all = TRUE), lapply(names(top_tables), function(k) {
  out <- top_tables[[k]]
  colnames(out)[colnames(out) %in% c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")] <-
    paste0(colnames(out)[colnames(out) %in% c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")], "___", k)
  out[, !(colnames(out) %in% c("contrast_id"))]
}))

cols_order        <- c("ensembl_gene_id", "gene_name", "gene_biotype",
                       sort(colnames(all_top_tables.df)[!(colnames(all_top_tables.df) %in% c("ensembl_gene_id", "gene_name", "gene_biotype"))]))
all_top_tables.df <- dplyr::select(all_top_tables.df, dplyr::all_of(cols_order))

saveRDS(top_tables, paste0(saveDir, "limma_top_tables_list.RDS"))

for (k in names(top_tables)){
  write.table(top_tables[[k]], paste0(saveDir, "limma_top_table___", k, ".tsv"), row.names = FALSE, sep = "\t")
}
saveRDS(all_top_tables.df, paste0(saveDir, "limma_top_table__all_merged.RDS"))
write.table(all_top_tables.df, paste0(saveDir, "limma_top_table__all_merged.tsv"), row.names = FALSE, sep = "\t")


