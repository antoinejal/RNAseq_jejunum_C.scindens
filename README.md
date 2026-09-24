# RNAseq_jejunum_C.scindens

## Summary

Intestinal lipid absorption plays a central role in systemic energy homeostasis. Using a precision microbial approach that modulates intestinal bile acid (BA) composition while preserving overall BA levels, we uncovered an unexpected role for the BA-7a-dehydroxylating bacterium Clostridium scindens in mitigating lipid absorption and diet-induced obesity. This protective effect was evident across distinct mouse models, from simplified gnotobiotic Oligo-MM12 mice to conventional mice with a complex microbiota. Mechanistically, Clostridium scindens colonization preserved the spatial zonation of lipid homeostasis genes along the crypt-villus axis, and normalized goblet cell numbers, both of which were disrupted by high-fat diet feeding. These coordinated changes maintained the mucus barrier and epithelial architecture essential for normal intestinal function. Overall, these findings reveal that colonization of a single BA-transforming commensal bacterium can counteract the maladaptive changes triggered by chronic fat exposure and protect against diet-induced obesity, providing a mechanistic framework for future therapeutic interventions in metabolic disease.

## File Description

Below, you can find descriptions of the files found under **scripts/**. The files were ran in the oredered in which they are described. 

| File | Description |
|---|---|
| `RNA_seq_processing.R` | Sample processing: quality control using FastQC, reference genome uploading and indexing followed by STAR alignment. Creation of a RNAseq object for downstream analyses. |
| `PCA.R` | Principal Component Analysis. |
| `DEA.R` | Differentially expressed gene analysis. |
| `GSEA.R` | Gene set enrichment analysis. |
| `ORA.R` | Over representation analysis. |
