# Clear working space
rm(list=ls())
# Clear console
cat("\014")
# Clear all open graphics devices
while (dev.cur() > 1) dev.off()

# Imports
library(stringr)
library(parallel)
library(data.table)
library(edgeR)

################################################################################
## QC of RNA-seq samples
################################################################################

# Define and create output directory
RNAseq_QC_Path <- "./data/RNAseq_processing/RNAseq_QC/samples"
if(!dir.exists(RNAseq_QC_Path)){
  dir.create(RNAseq_QC_Path, recursive = T)
}

nbCores     <- 30
FastQC_path <- "./tools/FastQC/fastqc"

# Samples
samples <- c("146", "150", "154", "155", "158", "161", "162", "164", "165", "168", "169", "175", "176", "179", "181")

# List fasta files, then assign names to have a named vector, where names are sample names
fqfiles         <- list.files("~/rcp_storage_archive/RAW/RNAseq/F24A430002208_MUSfseqR_Antoine_Cscindens_HFD", full.names = T, pattern = "fq\\.gz", recursive = T)
names(fqfiles)  <- gsub(".*\\/|\\.fq\\.gz", "", fqfiles)

# Keep only samples in the samples vector
fqfiles         <- fqfiles[gsub("_[12]$", "", names(fqfiles)) %in% samples]

# Loop over fasta files
nn <- parallel::mclapply(mc.cores = nbCores, X = names(fqfiles), FUN = function(ff){

  # Create saving sub-directory for current sample (current file)
  fqcdir <- paste0(RNAseq_QC_Path, "/", ff)
  if(!dir.exists(fqcdir)){
    dir.create(fqcdir, recursive = T)
  }
  
  # Check if output QC file already exists. If it is the case, skip the iteration and return NA (the code in a 
  # function stops after the return, that is the output of the function - in this case empty and the function 
  # stops without really doing anything)
  if(file.exists(paste0(fqcdir, "/", ff, "_fastqc.html"))){
    return()
  }
  
  # build fastQC command.
  fqcommand <- paste0(FastQC_path, " " , fqfiles[[ff]]," --threads 1 --outdir ", fqcdir)
  
  # run the command
  system(fqcommand)
})

################################################################################
## Aggregate FastQC and generate MultiQC report
################################################################################

RNAseq_QC_Path            <- "./data/RNAseq_processing/RNAseq_QC/samples"
RNAseq_QC_Path_aggregated <- "./data/RNAseq_processing/RNAseq_QC/summary"

if(!dir.exists(RNAseq_QC_Path_aggregated)){
  dir.create(RNAseq_QC_Path_aggregated, recursive = T)
}

if (!dir.exists("./tools/multiqc")) {
  system("pip install --upgrade --target ./tools/ multiqc")
}

command <- paste("PYTHONPATH=./tools python -m multiqc", normalizePath(RNAseq_QC_Path), "--outdir", normalizePath(RNAseq_QC_Path_aggregated), sep = " ")
system(command)

################################################################################
## Download reference genome files
################################################################################

saveDir <- "./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102"

if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

gtf_ftp_url              <- "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
primary_assembly_ftp_url <- "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"

if(!file.exists(paste0(saveDir, "/", gsub(".*\\/", "", gtf_ftp_url)))){
  system(paste0("wget -P ", saveDir, " ", gtf_ftp_url))
}

if(!file.exists(paste0(saveDir, "/", gsub(".*\\/", "", primary_assembly_ftp_url)))){
  system(paste0("wget -P ", saveDir, " ", primary_assembly_ftp_url))
}

################################################################################
## Generate genome (build STAR indexed genome from reference genome) 
## Use GRCm38.102
## The reads have length 150, therefore sjdbOverhang is set to 150 - 1
################################################################################

saveDir <- "./data/RNAseq_processing/reference_genome_STAR_indexed"

if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

if(!file.exists("./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.102.gtf")){
  R.utils::gunzip("./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.102.gtf.gz", remove = F)
}

if(!file.exists("./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa")){
  R.utils::gunzip("./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz", remove = F)
}

gtfref          <- "./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.102.gtf"
fastaref        <- "./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa"
starAlignerPath <- "./tools/STAR_2.7.11b/Linux_x86_64_static/STAR"
stargenome.ref  <- saveDir
nThreads        <- 50
readLength      <- 150
command         <- paste0(starAlignerPath,
                          " --runThreadN ", nThreads, 
                          " --runMode genomeGenerate",
                          " --genomeDir ", stargenome.ref, 
                          " --genomeFastaFiles ", fastaref, 
                          " --sjdbGTFfile ", gtfref, 
                          " --sjdbOverhang ", readLength - 1)
if(!file.exists(paste0(stargenome.ref, "/SA"))){
  system(command)
}

################################################################################
## STAR mapping / alignment
################################################################################

# STAR mapping requires to temporarily copy locally files otherwise the mapping fails
# We create a temporary folder on lispsrv1 then once the mapping is done we remove it
tmpLocalFolder_fasta  <- "./tmp_fastaFiles"
tmpLocalFolder_align  <- "./tmp_starAlign"
starAlignedPath       <- "./data/RNAseq_processing/star_alignment/GRCm38_release-102"

for (i in c(tmpLocalFolder_fasta, tmpLocalFolder_align, starAlignedPath)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}

starAlignerPath   <- "./tools/STAR_2.7.11b/Linux_x86_64_static/STAR"
indexedGenomePath <- "./data/RNAseq_processing/reference_genome_STAR_indexed"

samples_list <- sort(list.files("~/rcp_storage_archive/RAW/RNAseq/F24A430002208_MUSfseqR_Antoine_Cscindens_HFD", full.names = T, recursive = T, pattern = "fq\\.gz$"))
samples_list <- split(samples_list, gsub(".*\\/|_.*", "", samples_list))

# Samples
samples <- c("146", "150", "154", "155", "158", "161", "162", "164", "165", "168", "169", "175", "176", "179", "181")

# Keep only samples in the samples vector
samples_list <- samples_list[gsub("_[12]$", "", names(samples_list)) %in% samples]

nThreads <- 30
nn       <- mclapply(mc.cores = 2, X = names(samples_list), FUN = function(ss){

  tmpLocalFolder_align.ss <- paste0(tmpLocalFolder_align, "/", ss)
  
  if(!dir.exists(tmpLocalFolder_align.ss)){
    dir.create(tmpLocalFolder_align.ss, recursive = T)
  }
  
  fastaFiles <- stringr::str_sort(samples_list[[ss]], numeric = T) # sort string numerically to always have _1 first and _2 second
  
  # We always should have exactly 2 files (forward and reverse read) for one sample since we have paired end data
  if(length(fastaFiles) != 2){
    cat("Less than 2 fasta files detected for sample ", ss, "!! --> Aborting star mapping and skipping iteration...\n")
    return(NA)
  }
  
  # Skip if already mapped
  outputFile <- paste0(starAlignedPath, "/refalign_", ss, "_ReadsPerGene.out.tab")
  # if(file.exists(outputFile)){
  #   cat(paste0("Sample ", ss, " already mapped succesfully\n\n"))
  #   return(NA)
  # }
  
  # Get path of fasta files after copying them locally
  fastaFiles.local <- gsub(".*\\/", paste0(tmpLocalFolder_fasta, "/"), fastaFiles)
  
  cat("Making local copy of fastq files...\n")
  if(!all(file.exists(fastaFiles.local))){
    system(paste0("cp ", paste(fastaFiles, collapse = " "), " ", tmpLocalFolder_fasta))
  }
  
  
  cat("Running STAR alignment...\n")
  alignCommand <- paste0(starAlignerPath, 
                         " --runThreadN ", nThreads, 
                         " --genomeDir ", indexedGenomePath,
                         " --quantMode TranscriptomeSAM GeneCounts",
                         " --readFilesCommand zcat",
                         " --readFilesIn ", paste(fastaFiles.local, collapse = " "),
                         " --outFileNamePrefix ", tmpLocalFolder_align.ss, "/refAlign_", ss, "_")
  
  alignCommand <- paste0("bash -c \"ulimit -n 524288; exec ", alignCommand, "\"")
  system(alignCommand)
  
  # After mapping and saving files locally, move them to the lispnas1 folder (first list files, then build moving command then run command)
  cat("Alignment complete, transferring output to lispnas1 server...\n")
  starOutputs  <- list.files(tmpLocalFolder_align.ss, pattern = paste0("refAlign_", ss, "_"), full.names = T)
  
  # Check if all output files that should be there are present
  files.check <- paste0(tmpLocalFolder_align.ss, "/refAlign_", ss, "_", c("Log.final.out",
                                                                          "Log.out",
                                                                          "Log.progress.out",
                                                                          "ReadsPerGene.out.tab",
                                                                          "SJ.out.tab", 
                                                                          "Aligned.toTranscriptome.out.bam"))
  if(!all(files.check %in% starOutputs)){
    system(paste0("rm -r ", paste(starOutputs, collapse = " ")))

    # Remove temporary files copied locally to run star mapping analysis
    cat("Cleaning temporary files\n")
    system(paste0("rm ", paste(fastaFiles.local, collapse = " ")))

    return(NA)
  }
  
  moveCommands <- paste0("mv ", starOutputs[starOutputs %in% files.check], " ", starAlignedPath)
  sapply(moveCommands, system)
  
  # Remove temporary files copied locally to run star mapping analysis
  cat("Cleaning temporary files\n")
  system(paste0("rm -r ", tmpLocalFolder_align.ss))
  
  cat("Sample ", ss, " mapped successfully\n\n")
  
})

# Delete temporary local folders
unlink(tmpLocalFolder_fasta, recursive = T)
unlink(tmpLocalFolder_align, recursive = T)


# Check file size. This is just to see if there are not empty output files (just a few bites)
fileSize <- file.info(list.files(starAlignedPath, recursive = T, full.names = T))

################################################################################
## STAR alignment QC
################################################################################

starAlignedPath     <- "./data/RNAseq_processing/star_alignment/GRCm38_release-102"
star_qc_folderPath  <- "./data/RNAseq_processing/starAlignment_QC/GRCm38_release-102"
tmpFolderPath       <- paste0(star_qc_folderPath, "/tmp_logs")

for(i in c(star_qc_folderPath, tmpFolderPath)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}

# Get star mapping log files and copy them in a temporary folder
logFiles <- list.files(starAlignedPath, pattern = "Log.final.out$", recursive = T, full.names = T)
system(paste0("cp ", paste(logFiles, collapse = " "), " ", tmpFolderPath))

# Aggregate log files in a QC report with multiQC
if (!dir.exists("./tools/multiqc")) {
  system("pip install --upgrade --target ./tools/ multiqc")
}

command <- paste("PYTHONPATH=./tools python -m multiqc", normalizePath(tmpFolderPath), "--outdir", normalizePath(star_qc_folderPath), sep = " ")
system(command)

# Delete temporary log files folder
unlink(tmpFolderPath, recursive = T)

################################################################################
## Get STAR RNAseq gene-level counts
################################################################################

# Get reads per gene (rpg) file paths
rpg_paths        <- list.files("./data/RNAseq_processing/star_alignment/GRCm38_release-102", pattern = "ReadsPerGene", full.names = T, recursive = T)
names(rpg_paths) <- gsub(".*refAlign_|_ReadsPerGene.*", "", rpg_paths)


# The ReadsPerGene.out.tab output files of STAR (from option –quantMode GeneCounts) 
# contain 4 columns that correspond to different counts per gene calculated according 
# to the protocol’s strandedness
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA ((htseq-count 
#           option -s yes)). This corresponds to library stranded "forward", i.e.
#           the second cDNA strand was the one sequenced, i.e. it is the same as 
#           the original mRNA sample
# column 4: counts for the 2nd read strand aligned with RNA (the most common 
#           protocol nowadays, htseq-count option -s reverse). This corresponds to 
#           library stranded "reverse", i.e. the first cDNA strand is the one 
#           sequenced, so its sequence is the complement of the original mRNA
# If we have stranded data and choose one of the columns 3 or 4, the other 
# column (4 or 3) will give the count of antisense reads.
# https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html
# https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html 
# https://sydney-informatics-hub.github.io/training-RNAseq/06-CounttableToR/index.html
#
# Here we have have a stranded library (indicated in the BGI reports). If we do the sum
# of counts for each column in a ReadsPerGene file we get more counts in the last column
# of reads that map on both strands, i.e. reverse design. We therefore select column 4.
# We test here if it is the same for all samples

strandness.test <- lapply(names(rpg_paths), function(x){

  tmp.df <- fread(rpg_paths[[x]])
  tmp.df <- apply(tmp.df[-(1:4), -1], 2, sum)
  tmp.df <- as.data.frame(t(tmp.df))
  colnames(tmp.df) <- c("tot_unstranded", "tot_forward", "tot_reverse")
  tmp.df <- cbind(sample = x, tmp.df)
  tmp.df
})

strandness.test.df                   <- do.call(rbind, strandness.test)
strandness.test.df$is_reverse_design <- strandness.test.df$tot_reverse > strandness.test.df$tot_forward
strandness.test.df$rev_forw_ratio    <- strandness.test.df$tot_reverse / strandness.test.df$tot_forward

table(strandness.test.df$is_reverse_design)
hist(strandness.test.df$rev_forw_ratio, breaks = 20)

# Now build count table
allReads <- lapply(names(rpg_paths), FUN = function(x){
  
  # Select 1st and 4th columns
  out           <- fread(rpg_paths[x], data.table = F, stringsAsFactors = F)
  out           <- out[, c(1, 4)]
  colnames(out) <- c("gene", paste0("sample__", x))
  out
})

allGenes    <- unique(unlist(lapply(allReads, function(x){as.character(x$gene)})))
countMatrix <- do.call(cbind, lapply(allReads, function(x){
  out <- x[match(allGenes, x$gene), 2, drop = F]
  
  if(!identical(x[match(allGenes, x$gene), 1], allGenes)){
    stop("ERROR")
  }
  
  out
}))

rownames(countMatrix) <- allGenes
countMatrix           <- countMatrix[-(1:4), ]

dirSave <- "./data/RNAseq_processing/countMatrix/geneLevel"

if(!dir.exists(dirSave)){
  dir.create(dirSave, recursive = T)
}

# Saving count matrix
saveRDS(countMatrix, paste0(dirSave, "/STAR_RNAseq_count_matrix.RDS"))
write.table(countMatrix, paste0(dirSave, "/STAR_RNAseq_count_matrix.tsv"), sep = "\t", row.names = T, quote = F)

################################################################################
## Prepare gene conversion table from gtf files to convert between 
## Ensembl gene id and gene name 
################################################################################

saveDir <- "./data/RNAseq_processing/geneConversionTables"

if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

gtfref <- list.files("./data/input_data/reference_genome/Mus_musculus/GRCm38_release-102", recursive = T, full.names = T, pattern = "gtf.gz")
if(!file.exists(gsub(".gz$", "", gtfref))){
  R.utils::gunzip(gtfref, remove = F) # deflate
}

gtf    <- rtracklayer::import(str_replace(gtfref, ".gz", ""))
gtf_df <- type.convert(as.data.frame(gtf), as.is = T)

# Check column types and correct
gtf_df$score <- as.numeric(gtf_df$score)

geneConversionTable <- unique(dplyr::select(gtf_df, c("gene_id", "gene_name", "gene_biotype")))

write.table(geneConversionTable, paste0(saveDir, "/geneConversionTable_GRCm38_release-102.csv"), row.names = F, sep = ",")
saveRDS(geneConversionTable, paste0(saveDir, "/geneConversionTable_GRCm38_release-102.RDS"))

write.table(gtf_df, paste0(saveDir, "/geneConversionTable_GRCm38_release-102_complete_gtf.csv"), row.names = F, sep = ",")
saveRDS(gtf_df, paste0(saveDir, "/geneConversionTable_GRCm38_release-102_complete_gtf.RDS"))

################################################################################
## Format metadata
################################################################################

meta.df           <- openxlsx::read.xlsx("./data/input_data/metadata/raw_meta.xlsx")
colnames(meta.df) <- meta.df[1, ]
meta.df           <- dplyr::select(meta.df, c("Mouse ID", "Cage", "Group"))
meta.df           <- type.convert(meta.df[-1, ], as.is = T)

colnames(meta.df) <- c("mouse_id", "cage_id", "condition")

cols_num <- colnames(meta.df)[!(colnames(meta.df) %in% c("mouse_id", "cage_id", "condition"))]

for(i in cols_num){
  meta.df[[i]] <- as.numeric(meta.df[[i]])
}

meta.df$condition <- gsub(" ", "_", meta.df$condition)
meta.df           <- cbind(sample_id = gsub("KSC-008", "sample__", meta.df$mouse_id), meta.df)
meta.df           <- meta.df[!is.na(meta.df$sample_id), ]
rownames(meta.df) <- meta.df$sample_id

# Samples
samples <- c("146", "150", "154", "155", "158", "161", "162", "164", "165", "168", "169", "175", "176", "179", "181")

# Keep only samples in the samples vector
meta.df <- dplyr::filter(meta.df, gsub("sample__", "", sample_id) %in% samples)

saveRDS(meta.df, "./data/input_data/metadata/formatted_meta.RDS")

################################################################################
## Prepare object grouping counts, metadata, gene conversion table, ...
## Match dimnames
## Useful for all downstream analyses
################################################################################

saveDir <- "./data/RNAseq_processing/complete_data_obj"

if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

gene_conv.df <- readRDS("./data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102.RDS")
counts.df    <- readRDS("./data/RNAseq_processing/countMatrix/geneLevel/STAR_RNAseq_count_matrix.RDS")
meta.df      <- readRDS("./data/input_data/metadata/formatted_meta.RDS")

stopifnot(all(colnames(counts.df) %in% rownames(meta.df)))
common_samples <- intersect(colnames(counts.df), rownames(meta.df))
counts.df      <- counts.df[, common_samples]
meta.df        <- meta.df[common_samples, ]

# Initialize DGE object
y              <- edgeR::DGEList(counts = counts.df, samples = meta.df)

# Filter out lowly expressed genes
keep.exprs     <- edgeR::filterByExpr(y, group = y$samples$condition)
y              <- y[keep.exprs, , keep.lib.sizes=FALSE]

counts.df.filt <- as.data.frame(y$counts)
table(keep.exprs)

# Sanity checks
stopifnot(identical(colnames(counts.df), rownames(meta.df)))
stopifnot(identical(colnames(counts.df.filt), rownames(meta.df)))

data.list <- list(meta_df           = meta.df,
                  count_df          = counts.df,
                  count_cpm_df      = as.data.frame(edgeR::cpm(counts.df)),
                  count_filt_df     = counts.df.filt,
                  count_filt_cpm_df = as.data.frame(edgeR::cpm(counts.df.filt)),
                  gene_conv_df      = gene_conv.df)

saveRDS(data.list, paste0(saveDir, "/complete_RNAseq_data_object.RDS"))

