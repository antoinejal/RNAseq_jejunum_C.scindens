# Clear working space
rm(list=ls())
# Clear console
cat("\014")
# Clear all open graphics devices
while (dev.cur() > 1) dev.off()

# Imports
library(data.table)

################################################################################
## Over-representation Analysis
################################################################################

# Load DEA results from RDS
dea_list        <- readRDS("./data/DEA/limma_top_tables_list.RDS")
names(dea_list) <- gsub("\\.", "___", names(dea_list))

# Build all pairwise contrast combinations
comb.df           <- as.data.frame(t(combn(names(dea_list), 2)))
colnames(comb.df) <- c("contrast1", "contrast2")

gene_groups <- lapply(1:nrow(comb.df), function(x){

  cat(x, "/", nrow(comb.df), "...\n")
  
  # Load DEA tables for both contrasts
  df.dea.1                             <- dea_list[[comb.df$contrast1[x]]]
  df.dea.1                             <- dplyr::select(df.dea.1, dplyr::all_of(c("ensembl_gene_id", "gene_name", "logFC", "adj.P.Val")))
  colnames(df.dea.1)                   <- c("ensembl_gene_id", "gene_name", "logFC_contrast1", "adj.P.Val_contrast1")
  df.dea.1$direction_contrast1         <- paste0(ifelse(df.dea.1$logFC_contrast1 > 0, "Up", "Down"), "1")
  df.dea.1$significance_contrast1      <- paste0(ifelse(df.dea.1$adj.P.Val_contrast1 < 0.05, "Sig", "Not_Sig"), "1")
  
  df.dea.2                             <- dea_list[[comb.df$contrast2[x]]]
  df.dea.2                             <- dplyr::select(df.dea.2, dplyr::all_of(c("ensembl_gene_id", "gene_name", "logFC", "adj.P.Val")))
  colnames(df.dea.2)                   <- c("ensembl_gene_id", "gene_name", "logFC_contrast2", "adj.P.Val_contrast2")
  df.dea.2$direction_contrast2         <- paste0(ifelse(df.dea.2$logFC_contrast2 > 0, "Up", "Down"), "2")
  df.dea.2$significance_contrast2      <- paste0(ifelse(df.dea.2$adj.P.Val_contrast2 < 0.05, "Sig", "Not_Sig"), "2")
  
  df.dea.all <- merge(df.dea.1, df.dea.2, by = c("ensembl_gene_id", "gene_name"))
  df.dea.all <- df.dea.all[!(df.dea.all$significance_contrast1 == "Not_Sig1" & df.dea.all$significance_contrast2 == "Not_Sig2"), ]
  
  df.dea.all$color  <- paste0(df.dea.all$direction_contrast1, "__", df.dea.all$significance_contrast1, "____", df.dea.all$direction_contrast2, "__", df.dea.all$significance_contrast2)
  gene_groups       <- split(df.dea.all$ensembl_gene_id, df.dea.all$color)
  names(gene_groups)
  
  # Also group genes where significant in only one of the 2 contrasts
  grouping_categories <- list(c("Down1__Sig1____Down2__Not_Sig2", "Down1__Not_Sig1____Down2__Sig2"),
                              c("Up1__Sig1____Up2__Not_Sig2", "Up1__Not_Sig1____Up2__Sig2"),
                              c("Down1__Sig1____Up2__Not_Sig2", "Down1__Not_Sig1____Up2__Sig2"),
                              c("Up1__Sig1____Down2__Not_Sig2", "Up1__Not_Sig1____Down2__Sig2"))
  for(i in 1:length(grouping_categories)){
    if(all(grouping_categories[[i]] %in% names(gene_groups))){
      tmp <- unlist(gene_groups[names(gene_groups) %in% grouping_categories[[i]]])
      tmp <- list(tmp)
      names(tmp)  <- paste(grouping_categories[[i]], collapse = "__and__")
      gene_groups <- c(gene_groups, tmp)
    }
  }
  gene_groups[lengths(gene_groups) >= 5]
})
names(gene_groups) <- paste0(comb.df$contrast1, "__vs__", comb.df$contrast2)


# Saving
saveDir <- "./data/ORA/"

if (!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

gs_list <- readRDS("./data/input_data/msigdbr_gs_collections/gs_collection_list.RDS")

# Build universe from all DEA files
universe.genes <- lapply(dea_list, function(x) unique(x$ensembl_gene_id))
universe.genes <- sort(unique(unlist(universe.genes)))

overwrite <- TRUE
nbCores   <- 30

tt        <- lapply(names(gene_groups), function(x){

  message("Running ORA for ", x, "...")
  
  out_file <- paste0(saveDir, "/ORA_results_contrasts_logFC_comparison_gene_groups___", x, ".RDS")
  if (file.exists(out_file) & !overwrite){
    cat("Results already saved. Skipping...\n")
    return()
  }
  
  out.list <- lapply(names(gene_groups[[x]]), function(z){
    
    out.list <- parallel::mclapply(mc.cores = nbCores, X = names(gs_list), FUN = function(y){ 
      
      tryCatch({
        ora <- clusterProfiler::enricher(gene          = gene_groups[[x]][[z]],
                                         pvalueCutoff  = 1,
                                         pAdjustMethod = "BH",
                                         universe      = universe.genes,
                                         qvalueCutoff  = 1,
                                         maxGSSize     = 800,
                                         TERM2GENE     = dplyr::select(gs_list[[y]], dplyr::all_of(c("gs_name", "ensembl_gene"))))
        
        ora_df <- ora@result
        if (nrow(ora_df) == 0){
          return()
        }
        termsMeta  <- unique(dplyr::select(gs_list[[y]], dplyr::all_of(c("gs_name"))))
        ora_df     <- merge(ora_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = T, all.y = F)
        ora_df     <- ora_df[order(ora_df$pvalue), ]
        
        ora_df$GeneRatio_num <- unlist(lapply(ora_df$GeneRatio, function(w) eval(parse(text = w))))
        ora_df$BgRatio_num   <- unlist(lapply(ora_df$BgRatio, function(w) eval(parse(text = w))))
        
        ora_df$contrast_comparison_id <- x
        ora_df$gene_group_id          <- z
        ora_df$gs_collection          <- y
        rownames(ora_df)             <- ora_df$ID
        ora_df
      }, error = function(e){
        NULL
      })
      
    })
    if(all(unlist(lapply(out.list, is.null)))){
      return()
    }
    do.call(rbind, out.list)
  })
  if(all(unlist(lapply(out.list, is.null)))){
    return()
  }

  out.df  <- unique(rbindlist(out.list, use.names = T, fill = T))
  out.df  <- as.data.frame(out.df)
  
  saveRDS(out.df, out_file)
  write.table(out.df, gsub("RDS$", "tsv", out_file), sep = "\t", row.names = F)
})


