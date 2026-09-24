# Clear working space
rm(list=ls())
# Clear console
cat("\014")
# Clear all open graphics devices
while (dev.cur() > 1) dev.off()

# Imports
library(msigdbr)
library(data.table)

################################################################################
## Prepare gene set collections for functional enrichment
################################################################################

saveDir.gs <- "./data/input_data/msigdbr_gs_collections"

if(!dir.exists(saveDir.gs)){
  dir.create(saveDir.gs, recursive = T)
}

if(!file.exists(paste0(saveDir.gs, "/gs_collection_list.RDS"))){
  msigdbr.geneSets <- list(KEGG         = rbind(msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY"), msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_MEDICUS")),
                           REACTOME     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
                           WIKIPATHWAYS = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS"),
                           CGP          = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP"),
                           CP           = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP"),
                           BIOCARTA     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA"),
                           PID          = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:PID"),
                           TF           = rbind(msigdbr::msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD"), msigdbr::msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_LEGACY")),
                           BP           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
                           CC           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"),
                           MF           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF"),
                           HPO          = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "HPO"),
                           IMMUNESIGDB  = msigdbr::msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB"),
                           CELLTYPE     = msigdbr::msigdbr(species = "Mus musculus", category = "C8"),
                           HALLMARK     = msigdbr::msigdbr(species = "Mus musculus", category = "H"))
  
  # Add new datasets collected from enrichr
  gs_paths <- c("CellMarker_2024"                              = "./data/input_data/enrichr_gs_collections/CellMarker_2024.txt",
                "DisGeNET"                                     = "./data/input_data/enrichr_gs_collections/DisGeNET.txt",
                "Elsevier"                                     = "./data/input_data/enrichr_gs_collections/Elsevier_Pathway_Collection.txt",
                "GWAS_Catalog_2025"                            = "./data/input_data/enrichr_gs_collections/GWAS_Catalog_2025.txt",
                "HubMAP_ASCT_2022"                             = "./data/input_data/enrichr_gs_collections/HuBMAP_ASCTplusB_augmented_2022.txt",
                "Jensen_Compartments"                          = "./data/input_data/enrichr_gs_collections/Jensen_COMPARTMENTS.txt",
                "Jensen_Diseases"                              = "./data/input_data/enrichr_gs_collections/Jensen_DISEASES.txt",
                "KOMP2_2022"                                   = "./data/input_data/enrichr_gs_collections/KOMP2_Mouse_Phenotypes_2022.txt",
                "MGI_2024"                                     = "./data/input_data/enrichr_gs_collections/MGI_Mammalian_Phenotype_Level_4_2024.txt",
                "PanglaoDB_2021"                               = "./data/input_data/enrichr_gs_collections/PanglaoDB_Augmented_2021.txt",
                "Tabula_Sapiens"                               = "./data/input_data/enrichr_gs_collections/Tabula_Sapiens.txt")
  
  if(!file.exists("./data/input_data/enrichr_gs_collections/human_gene_table.tsv")){
    mart          <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attr          <- biomaRt::listAttributes(mart)
    human_genes   <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol", "external_synonym", "uniprot_gn_symbol"),
                                    mart       = mart,
                                    verbose    = T)
    human_genes[human_genes == ""] <- NA
    write.table(human_genes, "./data/input_data/enrichr_gs_collections/human_gene_table.tsv", sep = "\t", row.names = F)
  } else {
    human_genes <- fread("./data/input_data/enrichr_gs_collections/human_gene_table.tsv", data.table = F, stringsAsFactors = F)
  }
  human_genes <- human_genes[, colnames(human_genes) != "entrezgene_id"]
  human_genes <- type.convert(reshape2::melt(human_genes, id.vars = c("ensembl_gene_id")), as.is = T)
  human_genes <- unique(human_genes[!is.na(human_genes$value) & human_genes$value != "", ])
  
  if (!file.exists("./data/input_data/enrichr_gs_collections/mouse_gene_table.tsv")){
    mart          <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    attr          <- biomaRt::listAttributes(mart)
    mouse_genes   <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "mgi_symbol", "external_synonym", "uniprot_gn_symbol"),
                                    mart       = mart,
                                    verbose    = T)
    mouse_genes[mouse_genes == ""] <- NA
    write.table(mouse_genes, "./data/input_data/enrichr_gs_collections/mouse_gene_table.tsv", sep = "\t", row.names = F)
  } else {
    mouse_genes <- fread("./data/input_data/enrichr_gs_collections/mouse_gene_table.tsv", data.table = F, stringsAsFactors = F)
  }
  mouse_genes <- mouse_genes[, colnames(mouse_genes) != "entrezgene_id"]
  mouse_genes <- type.convert(reshape2::melt(mouse_genes, id.vars = c("ensembl_gene_id")), as.is = T)
  mouse_genes <- unique(mouse_genes[!is.na(mouse_genes$value) & mouse_genes$value != "", ])
  
  genesConv             <- fread("~/rcp_storage/common/Users/vonalven/BXD_Heart/Data/input_data/alliancegenome_geneTable_conversion/ORTHOLOGY-ALLIANCE_COMBINED_14_02_2024__formatted.tsv", header = T, data.table = F, stringsAsFactors = F)
  human2mouse           <- genesConv[genesConv$Gene1SpeciesName == "Homo sapiens" & genesConv$Gene2SpeciesName == "Mus musculus", ]
  human2mouse           <- unique(dplyr::select(human2mouse, c("Gene1EnsemblID", "Gene2EnsemblID")))
  colnames(human2mouse) <- c("hsapiens_gene_ensembl_id", "mmusculus_gene_ensembl_id")
  for(i in colnames(human2mouse)){
    human2mouse <- tidyr::separate_rows(human2mouse, c(i), sep = "\\,")
  }
  
  # By default, biomartr stores gene sets with human gene names, even if the gene sets are mouse-specific
  gs_list <- lapply(gs_paths, function(x){
    
    print(x)
    df.x               <- clusterProfiler::read.gmt(x)
    df.x$human_gene_id <- plyr::mapvalues(tolower(df.x$gene), from = tolower(human_genes$value), to = human_genes$ensembl_gene_id, warn_missing = F)
    df.x               <- df.x[grepl("^ENS", df.x$human_gene_id), ]
    df.x               <- merge(df.x, human2mouse, by.x = "human_gene_id", by.y = "hsapiens_gene_ensembl_id")
    df.x               <- unique(dplyr::select(df.x, dplyr::all_of(c("term", "mmusculus_gene_ensembl_id"))))
    df.x               <- df.x[!is.na(df.x$mmusculus_gene_ensembl_id) & df.x$mmusculus_gene_ensembl_id != "", ]
    colnames(df.x)     <- c("gs_name", "ensembl_gene")
    df.x               <- df.x[order(df.x$gs_name, df.x$ensembl_gene), ]
    df.x
  })
  
  
  custom_gs              <- openxlsx::read.xlsx("./data/input_data/enrichr_gs_collections/Custom dataset_haber_nature_2017_Zagoren_Development_2025_gene_signature_small_intestine2.xlsx")
  colnames(custom_gs)    <- paste0("custom__", tolower(gsub("\\(|\\)", "", gsub("\\.", "_", colnames(custom_gs)))))
  custom_gs              <- type.convert(reshape2::melt(as.matrix(custom_gs)), as.is = T)
  custom_gs              <- custom_gs[, -1]
  custom_gs$ensembl_gene <- plyr::mapvalues(tolower(custom_gs$value), from = tolower(mouse_genes$value), to = mouse_genes$ensembl_gene_id, warn_missing = F)
  custom_gs              <- custom_gs[, c(1, 3)]
  colnames(custom_gs)    <- c("gs_name", "ensembl_gene")
  custom_gs              <- custom_gs[!is.na(custom_gs$ensembl_gene), ]
  custom_gs              <- custom_gs[grepl("^ENS", custom_gs$ensembl_gene), ]
  custom_gs              <- list("Custom_Cell_Signatures" = custom_gs)
  
  msigdbr.geneSets <- c(msigdbr.geneSets, gs_list, custom_gs)
  saveRDS(msigdbr.geneSets, paste0(saveDir.gs, "/gs_collection_list.RDS"))
}

################################################################################
## Gene Set Enrichment Analysis (GSEA)
################################################################################

saveDir <- "./data/GSEA/"

if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = TRUE)
}

rnaseq_obj           <- readRDS("./data/RNAseq_processing/complete_data_obj/complete_RNAseq_data_object.RDS")
gene_table           <- rnaseq_obj$gene_conv_df
rownames(gene_table) <- gene_table$gene_id
gs_list              <- readRDS("./data/input_data/msigdbr_gs_collections/gs_collection_list.RDS")

# Load DEA results from RDS
dea_list        <- readRDS("./data/DEA/limma_top_tables_list.RDS")
names(dea_list) <- gsub("\\.", "___", names(dea_list))

nbCores <- 30

tt <- lapply(rev(names(dea_list)), function(x){
  
  # Extract contrast name for tracking
  contrast.x <- x
  
  message("Running GSEA analysis for contrast: ", contrast.x, "...")
  
  # Load the DEA table for this contrast
  dea_tbl <- as.data.frame(dea_list[[x]])
  
  out.list <- parallel::mclapply(mc.cores = nbCores, X = names(gs_list), FUN = function(y){
    
    tmp             <- dea_tbl
    tmp             <- tmp[!is.na(tmp$logFC), ]
    tmp             <- tmp[order(tmp$logFC, decreasing = TRUE), ]
    geneList        <- tmp$logFC
    names(geneList) <- tmp$ensembl_gene_id
    
    gsea <- clusterProfiler::GSEA(geneList      = geneList,
                                  nPermSimple   = 100000,
                                  pvalueCutoff  = 1,
                                  maxGSSize     = 800,
                                  pAdjustMethod = "BH",
                                  TERM2GENE     = dplyr::select(gs_list[[y]], dplyr::all_of(c("gs_name", "ensembl_gene"))))
    gsea_df <- gsea@result

    if(nrow(gsea_df) == 0){
      return()
    }

    termsMeta  <- unique(dplyr::select(gs_list[[y]], dplyr::all_of(c("gs_name"))))
    gsea_df    <- merge(gsea_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = TRUE, all.y = FALSE)
    gsea_df    <- gsea_df[order(gsea_df$pvalue), ]

    gsea_df$core_enrichment_gene_name <- unlist(lapply(gsea_df$core_enrichment, function(z){
      elements.z <- unlist(strsplit(z, "\\/"))
      paste(unique(gene_table[elements.z[elements.z %in% rownames(gene_table)], ]$gene_name), collapse = "/")
    }))
    
    rownames(gsea_df)         <- gsea_df$ID
    gsea_df$gs_collection     <- y
    gsea_df$contrast_id       <- contrast.x
    gsea_df
    
  })

  if (all(unlist(lapply(out.list, is.null)))){
    return()
  }

  out.df           <- unique(data.table::rbindlist(out.list, use.names = TRUE, fill = TRUE))
  out.df$coreSize  <- unlist(lapply(out.df$core_enrichment, function(x) length(unlist(strsplit(x, "/")))))
  out.df$geneRatio <- out.df$coreSize / out.df$setSize
  out.df           <- as.data.frame(out.df)
  
  saveRDS(out.df, paste0(saveDir, "GSEA_results_contrast_", contrast.x, ".RDS"))
  write.table(out.df, paste0(saveDir, "GSEA_results_contrast_", contrast.x, ".tsv"), sep = "\t", row.names = FALSE)
})
