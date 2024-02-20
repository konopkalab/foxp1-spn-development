# INTRONIC READS RATIO

### Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

### Compute Intronic Reads Ratio
```{R}
## Libraries
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)
library(harmony)
library(scCustomize)
library(DropletQC)
library(ggpubr)
options(future.globals.maxSize= 50000 * 1024^2)


## GTF 
REF_ANNOTATION <- "/work/Neuroinformatics_Core/akulk1/RESOURCES/DATABASES/GENCODEvM17/gencode.vM17.protein_coding.gtf"


## BAM PATHS
## WT_CTL
WT_CTL_P18_C1_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK/WT_1/outs/possorted_genome_bam.bam"
WT_CTL_P18_C2_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK/WT_2/outs/possorted_genome_bam.bam"
WT_CTL_P18_C3_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK/WT_3/outs/possorted_genome_bam.bam"

## CKO_CTL
CKO_CTL_P18_3BF1_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/CTL_NK3BF1/outs/possorted_genome_bam.bam"
CKO_CTL_P18_3BM3_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/CTL_NK3BM3/outs/possorted_genome_bam.bam"
CKO_CTL_P18_3BM4_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/CTL_NK3BM4/outs/possorted_genome_bam.bam"

## CKO_RES
CKO_RES_P18_3BM8_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/RES_NK3BM8/outs/possorted_genome_bam.bam"
CKO_RES_P18_3CF5_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/RES_NK3CF5/outs/possorted_genome_bam.bam"
CKO_RES_P18_3CF8_BAM <- "/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_JAN2022/RES_NK3CF8/outs/possorted_genome_bam.bam"


## ANNOTATED SEURAT OBJECT
load("SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED.RData")
# seuObj

as.data.frame(table(seuObj$GenoVirAgeSample))
#               Var1  Freq
# 7    WT_CTL_P18_C1  9507
# 8    WT_CTL_P18_C2 10996
# 9    WT_CTL_P18_C3  9206

# 1 CKO_CTL_P18_3BF1 10834
# 2 CKO_CTL_P18_3BM3 11390
# 3 CKO_CTL_P18_3BM4  7243

# 4 CKO_RES_P18_3BM8 20586
# 5 CKO_RES_P18_3CF5 10169
# 6 CKO_RES_P18_3CF8 10073

WT_CTL_P18_C1_CB <- gsub("P18_CTL_C1_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "WT_CTL_P18_C1",]))
WT_CTL_P18_C2_CB <- gsub("P18_CTL_C2_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "WT_CTL_P18_C2",]))
WT_CTL_P18_C3_CB <- gsub("P18_CTL_C3_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "WT_CTL_P18_C3",]))
CKO_CTL_P18_3BF1_CB <- gsub("P18_CTL_3BF1_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_CTL_P18_3BF1",]))
CKO_CTL_P18_3BM3_CB <- gsub("P18_CTL_3BM3_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_CTL_P18_3BM3",]))
CKO_CTL_P18_3BM4_CB <- gsub("P18_CTL_3BM4_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_CTL_P18_3BM4",]))
CKO_RES_P18_3BM8_CB <- gsub("P18_RES_3BM8_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_RES_P18_3BM8",]))
CKO_RES_P18_3CF5_CB <- gsub("P18_RES_3CF5_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_RES_P18_3CF5",]))
CKO_RES_P18_3CF8_CB <- gsub("P18_RES_3CF8_", "", row.names(seuObj@meta.data[seuObj$GenoVirAgeSample == "CKO_RES_P18_3CF8",]))


## CALCULATE INTRONIC READS RATIO
getIntronRatio <- function(bampath, cellbarcodes, prefx)
        {
        intron_ratio <- nuclear_fraction_annotation(
                annotation_path = REF_ANNOTATION, 
                bam = bampath, 
                bam_index = paste0(bampath, ".bai"),
                barcodes = cellbarcodes,
                tiles = 1,
                cores = 36,
                verbose = TRUE)
        write.table(intron_ratio, paste("NK_INTRONIC_READS_RATIO_P18_", prefx, ".txt", sep=""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
        return(intron_ratio)
        }

WT_CTL_P18_C1_IR <- getIntronRatio(WT_CTL_P18_C1_BAM, WT_CTL_P18_C1_CB, "WT_CTL_P18_C1")
WT_CTL_P18_C2_IR <- getIntronRatio(WT_CTL_P18_C2_BAM, WT_CTL_P18_C2_CB, "WT_CTL_P18_C2")
WT_CTL_P18_C3_IR <- getIntronRatio(WT_CTL_P18_C3_BAM, WT_CTL_P18_C3_CB, "WT_CTL_P18_C3")
CKO_CTL_P18_3BF1_IR <- getIntronRatio(CKO_CTL_P18_3BF1_BAM, CKO_CTL_P18_3BF1_CB, "CKO_CTL_P18_3BF1")
CKO_CTL_P18_3BM3_IR <- getIntronRatio(CKO_CTL_P18_3BM3_BAM, CKO_CTL_P18_3BM3_CB, "CKO_CTL_P18_3BM3")
CKO_CTL_P18_3BM4_IR <- getIntronRatio(CKO_CTL_P18_3BM4_BAM, CKO_CTL_P18_3BM4_CB, "CKO_CTL_P18_3BM4")
CKO_RES_P18_3BM8_IR <- getIntronRatio(CKO_RES_P18_3BM8_BAM, CKO_RES_P18_3BM8_CB, "CKO_RES_P18_3BM8")
CKO_RES_P18_3CF5_IR <- getIntronRatio(CKO_RES_P18_3CF5_BAM, CKO_RES_P18_3CF5_CB, "CKO_RES_P18_3CF5")
CKO_RES_P18_3CF8_IR <- getIntronRatio(CKO_RES_P18_3CF8_BAM, CKO_RES_P18_3CF8_CB, "CKO_RES_P18_3CF8")

save(WT_CTL_P18_C1_IR, WT_CTL_P18_C2_IR, WT_CTL_P18_C3_IR, CKO_CTL_P18_3BF1_IR, CKO_CTL_P18_3BM3_IR, CKO_CTL_P18_3BM4_IR, CKO_RES_P18_3BM8_IR, CKO_RES_P18_3CF5_IR, CKO_RES_P18_3CF8_IR, file = "NK_INTRONIC_READS_RATIO_P18.RData")

```

