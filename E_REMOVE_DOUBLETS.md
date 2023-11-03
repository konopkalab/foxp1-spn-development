# Removing Doublets using DoubletFinder


## WT-CTL

### Load modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

### Run Seurat Pipeline (without filtering data)
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## load dataset from 10X run
p18.ctl.c1 <- Read10X_h5("CELLBENDER_SEP2022/WT_1_CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.ctl.c2 <- Read10X_h5("CELLBENDER_SEP2022/WT_2_CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.ctl.c3 <- Read10X_h5("CELLBENDER_SEP2022/WT_3_CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)

## adding sample & condition info to column names
colnames(p18.ctl.c1) <- paste("P18_CTL_C1", colnames(p18.ctl.c1), sep = "_")
colnames(p18.ctl.c2) <- paste("P18_CTL_C2", colnames(p18.ctl.c2), sep = "_")
colnames(p18.ctl.c3) <- paste("P18_CTL_C3", colnames(p18.ctl.c3), sep = "_")

## fetch genes or rows corresponding to gencode protein coding gene symbols
p18.ctl.c1.temp <- p18.ctl.c1[row.names(p18.ctl.c1) %in% ref$GeneSymbol,]
p18.ctl.c2.temp <- p18.ctl.c2[row.names(p18.ctl.c2) %in% ref$GeneSymbol,]
p18.ctl.c3.temp <- p18.ctl.c3[row.names(p18.ctl.c3) %in% ref$GeneSymbol,]

## make a data from from matrix
NK_P18_CTL_C1 <- as.data.frame(as.matrix(p18.ctl.c1.temp))
NK_P18_CTL_C2 <- as.data.frame(as.matrix(p18.ctl.c2.temp))
NK_P18_CTL_C3 <- as.data.frame(as.matrix(p18.ctl.c3.temp))

## add gene symbol as a column
NK_P18_CTL_C1$Genes <- row.names(NK_P18_CTL_C1)
NK_P18_CTL_C2$Genes <- row.names(NK_P18_CTL_C2)
NK_P18_CTL_C3$Genes <- row.names(NK_P18_CTL_C3)

## combine individual tables into a giant data frame
dataCombined <- list("NK_P18_CTL_C1" = NK_P18_CTL_C1, "NK_P18_CTL_C2" = NK_P18_CTL_C2, "NK_P18_CTL_C3" = NK_P18_CTL_C3)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

nk.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

nk.data[is.na(nk.data)] <- 0
p18.ctl.data <- nk.data

## Save Counts as RData
save(p18.ctl.data, file = "NK_P18_WT_CTL_CB_COUNTS.RData")

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p18.ctl.data, project = "NK_P18_WT_CTL_CB")

## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 4, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "SEURAT_NK_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

## Data Normalization
seuObj <- NormalizeData(seuObj)

## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)
ggsave(filename = "SEURAT_NK_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

## Data Scaling
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)
ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## Identify significant PCs
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## Data Clustering
seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

## Save RData
save(seuObj, file = "SEURAT_NK_P18_WT_CTL_CB.RData")
```

### Run DoubletFinder to identify potential doublets
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("SEURAT_NK_P18_WT_CTL_CB.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## pK Identification (no ground-truth)
sweep.ctl.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.ctl.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) 
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = colnames(seuObjDF@meta.data)[grepl("^pANN", colnames(seuObjDF@meta.data))], sct = FALSE)

dfcols <- c(colnames(seuObjDF@meta.data)[grepl("^DF.classification", colnames(seuObjDF@meta.data))])

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta[,dfcols[[1]]] == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp[,dfcols[[2]]] == "Singlet",]
seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

## Save RData
save(seuObjDF, selpcs, selpk, nExp_poi.adj, file = "SEURAT_NK_P18_WT_CTL_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "SEURAT_NK_P18_WT_CTL_CB_DF_FILT_META.RData")

```


<p>&nbsp;</p>


## CKO-CTL

### Load modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.0.2-gccmkl
```

### Run Seurat Pipeline (without filtering data)
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## Load CellBender count tables
p18.ctl.3cf1 <- Read10X_h5("CELLBENDER_JAN2022/NK3CF1_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.ctl.3bm4 <- Read10X_h5("CELLBENDER_JAN2022/NK3BM4_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.ctl.3bf1 <- Read10X_h5("CELLBENDER_JAN2022/NK3BF1_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.ctl.3bm3 <- Read10X_h5("CELLBENDER_JAN2022/NK3BM3_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)

## adding sample & condition info to column names
colnames(p18.ctl.3cf1) <- paste("P18_CTL_3CF1", colnames(p18.ctl.3cf1), sep = "_")
colnames(p18.ctl.3bm4) <- paste("P18_CTL_3BM4", colnames(p18.ctl.3bm4), sep = "_")
colnames(p18.ctl.3bf1) <- paste("P18_CTL_3BF1", colnames(p18.ctl.3bf1), sep = "_")
colnames(p18.ctl.3bm3) <- paste("P18_CTL_3BM3", colnames(p18.ctl.3bm3), sep = "_")

## fetch genes or rows corresponding to gencode protein coding gene symbols
p18.ctl.3cf1.temp <- p18.ctl.3cf1[row.names(p18.ctl.3cf1) %in% ref$GeneSymbol,]
p18.ctl.3bm4.temp <- p18.ctl.3bm4[row.names(p18.ctl.3bm4) %in% ref$GeneSymbol,]
p18.ctl.3bf1.temp <- p18.ctl.3bf1[row.names(p18.ctl.3bf1) %in% ref$GeneSymbol,]
p18.ctl.3bm3.temp <- p18.ctl.3bm3[row.names(p18.ctl.3bm3) %in% ref$GeneSymbol,]

## make a data from from matrix
NK_P18_CTL_3CF1 <- as.data.frame(as.matrix(p18.ctl.3cf1.temp))
NK_P18_CTL_3BM4 <- as.data.frame(as.matrix(p18.ctl.3bm4.temp))
NK_P18_CTL_3BF1 <- as.data.frame(as.matrix(p18.ctl.3bf1.temp))
NK_P18_CTL_3BM3 <- as.data.frame(as.matrix(p18.ctl.3bm3.temp))

## add gene symbol as a column
NK_P18_CTL_3CF1$Genes <- row.names(NK_P18_CTL_3CF1)
NK_P18_CTL_3BM4$Genes <- row.names(NK_P18_CTL_3BM4)
NK_P18_CTL_3BF1$Genes <- row.names(NK_P18_CTL_3BF1)
NK_P18_CTL_3BM3$Genes <- row.names(NK_P18_CTL_3BM3)

## combine individual tables into a giant data frame
dataCombined <- list("NK_P18_CTL_3CF1" = NK_P18_CTL_3CF1, "NK_P18_CTL_3BM4" = NK_P18_CTL_3BM4, "NK_P18_CTL_3BF1" = NK_P18_CTL_3BF1, "NK_P18_CTL_3BM3" = NK_P18_CTL_3BM3)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

nk.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

nk.data[is.na(nk.data)] <- 0
p18.ctl.data <- nk.data

## Save Counts as RData
save(p18.ctl.data, file = "NK_P18_CKO_CTL_CB_COUNTS.RData")

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p18.ctl.data, project = "NK_P18_CKO_CTL_CB")

## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 4, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "SEURAT_NK_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

## Data Normalization
seuObj <- NormalizeData(seuObj)

## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)
ggsave(filename = "SEURAT_NK_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

## Data Scaling
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)
ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## Identify significant PCs
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## Data Clustering
seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

## Save RData
save(seuObj, file = "SEURAT_NK_P18_CKO_CTL_CB.RData")
```

### Run DoubletFinder to identify potential doublets
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("SEURAT_NK_P18_CKO_CTL_CB.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## pK Identification (no ground-truth)
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) 
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_4550", sct = FALSE)

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_4550 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_4264 == "Singlet",]
seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

## Save RData
save(seuObjDF, file = "SEURAT_NK_P18_CKO_CTL_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "SEURAT_NK_P18_CKO_CTL_CB_DF_FILT_META.RData")

```


<p>&nbsp;</p>



## CKO-RES

### Load modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.0.2-gccmkl
```

### Run Seurat Pipeline (without filtering data)
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
set.seed(10)

## reference gene symbol list
ref <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## load dataset from 10X run
p18.res.3cf5 <- Read10X_h5("CELLBENDER_JAN2022/NK3CF5_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.res.3cf8 <- Read10X_h5("CELLBENDER_JAN2022/NK3CF8_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.res.3bf5 <- Read10X_h5("CELLBENDER_JAN2022/NK3BF5_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)
p18.res.3bm8 <- Read10X_h5("CELLBENDER_JAN2022/NK3BM8_raw_feature_bc_matrix/CellBender_Out_filtered.h5", use.names = TRUE, unique.features = TRUE)

## adding sample & condition info to column names
colnames(p18.res.3cf5) <- paste("P18_RES_3CF5", colnames(p18.res.3cf5), sep = "_")
colnames(p18.res.3cf8) <- paste("P18_RES_3CF8", colnames(p18.res.3cf8), sep = "_")
colnames(p18.res.3bf5) <- paste("P18_RES_3BF5", colnames(p18.res.3bf5), sep = "_")
colnames(p18.res.3bm8) <- paste("P18_RES_3BM8", colnames(p18.res.3bm8), sep = "_")

## fetch genes or rows corresponding to gencode protein coding gene symbols
p18.res.3cf5.temp <- p18.res.3cf5[row.names(p18.res.3cf5) %in% ref$GeneSymbol,]
p18.res.3cf8.temp <- p18.res.3cf8[row.names(p18.res.3cf8) %in% ref$GeneSymbol,]
p18.res.3bf5.temp <- p18.res.3bf5[row.names(p18.res.3bf5) %in% ref$GeneSymbol,]
p18.res.3bm8.temp <- p18.res.3bm8[row.names(p18.res.3bm8) %in% ref$GeneSymbol,]

## make a data from from matrix
NK_P18_RES_3CF5 <- as.data.frame(as.matrix(p18.res.3cf5.temp))
NK_P18_RES_3CF8 <- as.data.frame(as.matrix(p18.res.3cf8.temp))
NK_P18_RES_3BF5 <- as.data.frame(as.matrix(p18.res.3bf5.temp))
NK_P18_RES_3BM8 <- as.data.frame(as.matrix(p18.res.3bm8.temp))

## add gene symbol as a column
NK_P18_RES_3CF5$Genes <- row.names(NK_P18_RES_3CF5)
NK_P18_RES_3CF8$Genes <- row.names(NK_P18_RES_3CF8)
NK_P18_RES_3BF5$Genes <- row.names(NK_P18_RES_3BF5)
NK_P18_RES_3BM8$Genes <- row.names(NK_P18_RES_3BM8)

## combine individual tables into a giant data frame
dataCombined <- list("NK_P18_RES_3CF5" = NK_P18_RES_3CF5, "NK_P18_RES_3CF8" = NK_P18_RES_3CF8, "NK_P18_RES_3BF5" = NK_P18_RES_3BF5, "NK_P18_RES_3BM8" = NK_P18_RES_3BM8)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

nk.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

nk.data[is.na(nk.data)] <- 0
p18.res.data <- nk.data

## Save Counts as RData
save(p18.res.data, file = "NK_P18_CKO_RES_CB_COUNTS.RData")

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = p18.res.data, project = "NK_P18_CKO_RES_CB")

## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 4, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Age", "Genotype", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaAge <- metaData$Age
names(metaAge) <- row.names(metaData)

metaGeno <- metaData$Genotype
names(metaGeno) <- row.names(metaData)

seuObj$Genotype <- metaGeno
seuObj$Sample <- metaSample
seuObj$Age <- metaAge
seuObj$GenoAgeSample <- paste(seuObj$Genotype, seuObj$Age, seuObj$Sample, sep = "_")

## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "GenoAgeSample"

## Visualize Data QC
p9 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "SEURAT_NK_QC_1.pdf", plot = p9, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

## Data Normalization
seuObj <- NormalizeData(seuObj)

## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top90 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top90, repel = TRUE)
ggsave(filename = "SEURAT_NK_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

## Data Scaling
seuObj <- ScaleData(seuObj, vars.to.regress = c("nCount_RNA"))

## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)
ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## Identify significant PCs
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## Data Clustering
seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

## Save RData
save(seuObj, file = "SEURAT_NK_P18_CKO_RES_CB.RData")

```

### Run DoubletFinder to identify potential doublets
```{R}
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

load("SEURAT_NK_P18_CKO_RES_CB.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1

## pK Identification (no ground-truth)
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(seuObj@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_6131", sct = FALSE)

seuMeta <- as.data.frame(seuObjDF@meta.data)
seuMetaFiltTemp <- seuMeta[seuMeta$DF.classifications_0.25_0.005_6131 == "Singlet",]
seuMetaFilt <- seuMetaFiltTemp[seuMetaFiltTemp$DF.classifications_0.25_0.005_5810 == "Singlet",]
seuMetaFilt$DoubletFinder <- rep("Retained", nrow(seuMetaFilt))

df.meta <- seuMetaFilt$DoubletFinder
names(df.meta) <- row.names(seuMetaFilt)

seuObjDF$DoubletFinder <- df.meta
seuObjDF$DoubletFinder[is.na(seuObjDF$DoubletFinder)] <- "Discarded"

## Save RData
save(seuObjDF, file = "SEURAT_NK_P18_CKO_RES_CB_DF.RData")
save(seuMeta, seuMetaFilt, file = "SEURAT_NK_P18_CKO_RES_CB_DF_FILT_META.RData")

```

<p>&nbsp;</p>

