# SPN SUB-CLUSTERING

### Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

### SPN Sub-clustering
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
options(future.globals.maxSize= 50000 * 1024^2)

## load seurat object
load("./../SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED.RData")


Idents(seuObj) <- "CellType"
print(table(seuObj@meta.data$CellType))
#             Astrocytes               Cortical            Endothelial 
#                  14639                    925                   2717 
#           Interneurons              Microglia Neurogenic_Progenitors 
#                    720                   4497                   4578 
#       Oligodendrocytes            Progenitors                   SPNs 
#                  19095                   1634                  41865 
#           Undetermined 
#                   9334

# dSPN <- c(0, 1, 3) # 11048 + 9407 + 6325
# iSPN <- c(4, 9, 25, 32) # 4998 + 3490 + 1001 + 652
# eSPN <- c(10, 18) # 2691 + 1639

# ident2keep <- sort(c(dSPN, iSPN, eSPN))

subObj <- subset(seuObj, idents = "SPNs", slot = "counts")
subObj$GenoVir <- paste(subObj$Genotype, subObj$Virus, sep = "_")

## preprocess
subObj <- NormalizeData(subObj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE, npcs = 100)

## jackstraw analysis
## Compute and Score Jackstraw
subObj <- JackStraw(subObj, num.replicate = 100, dims = 100)
subObj <- ScoreJackStraw(subObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(subObj, ndims = 100)
p7 <- JackStrawPlot(subObj, dims = 1:100)

ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## IDENTIFY PCS
pcScores <- subObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1


## harmony
subObj <- RunHarmony(subObj, group.by.vars = c("Batch", "Sex"), reduction = "pca", assay.use = "RNA", max.iter.harmony = 100)

## clustering
subObj <- FindNeighbors(subObj, dims = 1:selpcs, reduction = "harmony")
subObj <- FindClusters(subObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), reduction = "harmony")
subObj <- RunUMAP(subObj, dims = 1:selpcs, reduction = "harmony")

## SAVE RDATA
save(subObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY.RData")



## CLUSTREE
seutree <- clustree(subObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SEURAT_NK_INTEGRATE_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## UMAP FOR ALL RESOLUTIONS
ResolutionList <- grep("RNA_snn_res", colnames(subObj@meta.data), value = TRUE)

for (Resolution in ResolutionList)
    {
    print(paste("====> RESOLUTION ", Resolution, sep = ""))

    pdf(paste0("SEURAT_NK_INTEGRATE_UMAP_RES_", Resolution, ".pdf"), width=7, height=6)
    g <- DimPlot(object = subObj, label = TRUE, reduction = "umap", group.by = Resolution)
    print(g)
    dev.off()

    pdf(paste0("SEURAT_NK_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width=7, height=3)
    v <- VlnPlot(object = subObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
    print(v)
    dev.off()
    }



## UMAP FOR MARKER GENES
mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cspg4", "Mag", "Cx3cr1", "Flt1", "Slc17a7", "Chat", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b", "Oprm1", "Isl1", "Pdyn", "Lypd1", "Nnat", "Ebf1", "Epha4", "Mef2c")

fpl1 <- FeaturePlot_scCustom(seurat_object = subObj, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 15, height = 24, units = "in", dpi = 150)

fpl2 <- FeaturePlot_scCustom(seurat_object = subObj, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 15, height = 24, units = "in", dpi = 150)


plotCluUMAP1 <- DimPlot_scCustom(seurat_object = subObj, group.by = "GenoVir", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOTYPE.pdf", plot = plotCluUMAP1, width = 10, height = 8, units = "in", dpi = 300)

plotCluUMAP1b <- DimPlot_scCustom(seurat_object = subObj, group.by = "GenoVir", split.by = "GenoVir", num_columns = 3, pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOTYPE_FACET.pdf", plot = plotCluUMAP1b, width = 24, height = 8, units = "in", dpi = 300)

plotCluUMAP1c <- DimPlot_scCustom(seurat_object = subObj, group.by = "GenoVirAgeSample", split.by = "GenoVirAgeSample", num_columns = 3, pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOVIRAGESAMPLE_FACET.pdf", plot = plotCluUMAP1c, width = 25, height = 20, units = "in", dpi = 300)

plotCluUMAP1d <- DimPlot_scCustom(seurat_object = subObj, group.by = "GenoVirAgeSample", pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOVIRAGESAMPLE.pdf", plot = plotCluUMAP1d, width = 8, height = 6, units = "in", dpi = 300)

plotCluUMAP1e <- DimPlot(subObj, group.by = c("Sex", "Batch", "Genotype", "Virus", "GenoVirAgeSample"), pt.size = 0.01)
ggsave("NK_SEURAT_UMAP_COVARIATES.pdf", plot = plotCluUMAP1e, width = 20, height = 16, units = "in", dpi = 300)



## RESOLUTION 1.2

my_levels <- seq(0, max(as.numeric(subObj$RNA_snn_res.1.2)) - 1)
subObj$RNA_snn_res.1.2 <- factor(x = subObj$RNA_snn_res.1.2, levels = my_levels)

vpl2 <- Stacked_VlnPlot(seurat_object = subObj, features = mygenes, group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_ViolinPlot_1.2_A.pdf", plot = vpl2 , width = 12, height = 18, units = "in", dpi = 150, useDingbats=FALSE)

vpl3 <- Stacked_VlnPlot(seurat_object = subObj, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_ViolinPlot_1.2_B.pdf", plot = vpl3 , width = 12, height = 3, units = "in", dpi = 150, useDingbats=FALSE)

dpl2 <- DotPlot_scCustom(seurat_object = subObj, features = rev(mygenes), colors_use = viridis_plasma_dark_high, group.by = "RNA_snn_res.1.2", x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_DotPlot_1.2.pdf", plot = dpl2, width = 12, height = 16, units = "in", dpi = 150, useDingbats=FALSE)


## RES 1.2
## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
cellsPerCluster <- as.data.frame.matrix(table(subObj@meta.data$RNA_snn_res.1.2, subObj@meta.data$GenoVir))
cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
write.table(cellsPerCluster, paste("NK_SEURAT_GENOTYPE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

cellsPerCluster2 <- melt(cellsPerCluster)
colnames(cellsPerCluster2) <- c("CLUSTER", "GENOTYPE", "CELLS")
p7 <- ggplot(cellsPerCluster2) +
        geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOTYPE), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = paste("NK_SEURAT_GENOTYPE_PER_CLUSTER_1.2.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)




## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
cellsPerCluster <- as.data.frame.matrix(table(subObj@meta.data$RNA_snn_res.1.2, subObj@meta.data$GenoVirAgeSample))
cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
cellsPerCluster <- cellsPerCluster[order(cellsPerCluster$Cluster),]
write.table(cellsPerCluster, paste("NK_SEURAT_GENOVIRAGESAMPLE_PER_CLUSTER_1.2.txt", sep = "_"), row.names = F, col.names = T, quote = F, sep = "\t")

cellsPerCluster2 <- melt(cellsPerCluster)
colnames(cellsPerCluster2) <- c("CLUSTER", "GENOVIRAGESAMPLE", "CELLS")
p7 <- ggplot(cellsPerCluster2) +
        geom_bar(aes(x = CLUSTER, y = CELLS, fill = GENOVIRAGESAMPLE), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = paste("NK_SEURAT_GENOVIRAGESAMPLE_PER_CLUSTER_1.2.pdf", sep = "_"), plot = p7, width = 9, height = 3, units = "in", dpi = 150)

```


### Cluster Markers
```{R}
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
library(scCustomize)
library(viridis)

load("SEURAT_NK_P18_CKO_RES_HARMONY.RData")
# subObj

seuObj.markers <- FindAllMarkers(subObj, only.pos = TRUE) #, min.pct = 0, logfc.threshold = 0)
write.table(subObj.markers, "NK_SEURAT_INTEGRATION_CLUST_1.2_CLUSTER_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

```


### WPRE MCHERRY
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
options(future.globals.maxSize= 50000 * 1024^2)


##------------------------------------
## load Seurat Harmony object
load("SEURAT_NK_P18_CKO_RES_HARMONY.RData")
# subObj
seuObj <- subObj


## WPRE | SEPT 2022
c1wpre <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK_WPREMCHERRY/WT_1_Amp/outs/filtered_feature_bc_matrix")
colnames(c1wpre) <- paste("P18_CTL_C1", colnames(c1wpre), sep = "_")
c1wpre2 <- as.data.frame(c1wpre["WPRE",])
colnames(c1wpre2) <- "WPRE"

c2wpre <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK_WPREMCHERRY/WT_2_Amp/outs/filtered_feature_bc_matrix")
colnames(c2wpre) <- paste("P18_CTL_C2", colnames(c2wpre), sep = "_")
c2wpre2 <- as.data.frame(c2wpre["WPRE",])
colnames(c2wpre2) <- "WPRE"

c3wpre <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/02_COUNT_NK_WPREMCHERRY/WT_3_Amp/outs/filtered_feature_bc_matrix")
colnames(c3wpre) <- paste("P18_CTL_C3", colnames(c3wpre), sep = "_")
c3wpre2 <- as.data.frame(c3wpre["WPRE",])
colnames(c3wpre2) <- "WPRE"

wpreDataSept2022 <- rbind(c1wpre2, c2wpre2, c3wpre2)



## WPRE | JAN 2022
nk3bf1c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3BF1c/outs/filtered_feature_bc_matrix")
colnames(nk3bf1c) <- paste("P18_CTL_3BF1", colnames(nk3bf1c), sep = "_")
nk3bf1c2 <- as.data.frame(nk3bf1c["WPRE",])
colnames(nk3bf1c2) <- "WPRE"

nk3bm3c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3BM3c/outs/filtered_feature_bc_matrix")
colnames(nk3bm3c) <- paste("P18_CTL_3BM3", colnames(nk3bm3c), sep = "_")
nk3bm3c2 <- as.data.frame(nk3bm3c["WPRE",])
colnames(nk3bm3c2) <- "WPRE"

nk3bm4c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3BM4c/outs/filtered_feature_bc_matrix")
colnames(nk3bm4c) <- paste("P18_CTL_3BM4", colnames(nk3bm4c), sep = "_")
nk3bm4c2 <- as.data.frame(nk3bm4c["WPRE",])
colnames(nk3bm4c2) <- "WPRE"

wpreDataJan2022 <- rbind(nk3bf1c2, nk3bm3c2, nk3bm4c2)



## MCHERRY | JAN 2022
nk3bm8c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3BM8c/outs/filtered_feature_bc_matrix")
colnames(nk3bm8c) <- paste("P18_RES_3BM8", colnames(nk3bm8c), sep = "_")
nk3bm8c2 <- as.data.frame(nk3bm8c["MCHERRY",])
colnames(nk3bm8c2) <- "MCHERRY"

nk3cf5c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3CF5c/outs/filtered_feature_bc_matrix")
colnames(nk3cf5c) <- paste("P18_RES_3CF5", colnames(nk3cf5c), sep = "_")
nk3cf5c2 <- as.data.frame(nk3cf5c["MCHERRY",])
colnames(nk3cf5c2) <- "MCHERRY"

nk3cf8c <- Read10X("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/03_COUNT/NK_MAY2022_WPRE_MCHERRY/NK3CF8c/outs/filtered_feature_bc_matrix")
colnames(nk3cf8c) <- paste("P18_RES_3CF8", colnames(nk3cf8c), sep = "_")
nk3cf8c2 <- as.data.frame(nk3cf8c["MCHERRY",])
colnames(nk3cf8c2) <- "MCHERRY"

mcherryDataJan2022 <- rbind(nk3bm8c2, nk3cf5c2, nk3cf8c2)


wpreData <- rbind(wpreDataJan2022, wpreDataSept2022)
mcherryData <- mcherryDataJan2022


## Update Seurat Object
metaTemp <- as.data.frame(seuObj@meta.data)

wpreCommon <- merge(metaTemp, wpreData, by = "row.names", all.x = TRUE)
mcherryCommon <- merge(metaTemp, mcherryData, by = "row.names", all.x = TRUE)


metaWPRE <- wpreCommon$WPRE
names(metaWPRE) <- wpreCommon$Row.names

metaMCHERRY <- mcherryCommon$MCHERRY
names(metaMCHERRY) <- mcherryCommon$Row.names

seuObj$WPRE <- metaWPRE
seuObj$WPRE[is.na(seuObj$WPRE)] <- 0
seuObj$MCHERRY <- metaMCHERRY
seuObj$MCHERRY[is.na(seuObj$MCHERRY)] <- 0


# pumap1a <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "WPRE", raster = TRUE, label = FALSE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$WPRE))), palette = "varibow")) #+ NoLegend()
# ggsave("NK_SEURAT_PLOT_UMAP_WPRE.PDF", plot = pumap1a, width = 6.5, height = 6, units = "in", dpi = 150)

# pumap2a <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "MCHERRY", raster = TRUE, label = FALSE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$MCHERRY))), palette = "varibow")) #+ NoLegend()
# ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY.PDF", plot = pumap2a, width = 8.5, height = 6, units = "in", dpi = 150)


## WPRE | MCHERRY >= 1
seuObj$WPRE_GTEQ_1 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE)
seuObj$MCHERRY_GTEQ_1 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_1", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_1.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_1", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_1.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)


## WPRE | MCHERRY >= 2
seuObj$WPRE_GTEQ_2 <- gsub("^1$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_2 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_2)
seuObj$MCHERRY_GTEQ_2 <- gsub("^1$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_2 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_2)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_2", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_2.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_2", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_2.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)


## WPRE | MCHERRY >= 3
seuObj$WPRE_GTEQ_3 <- gsub("^1$|^2$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_3 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_3)
seuObj$MCHERRY_GTEQ_3 <- gsub("^1$|^2$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_3 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_3)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_3", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_3.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_3", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_3.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)


## WPRE | MCHERRY >= 4
seuObj$WPRE_GTEQ_4 <- gsub("^1$|^2$|^3$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_4 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_4)
seuObj$MCHERRY_GTEQ_4 <- gsub("^1$|^2$|^3$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_4 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_4)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_4", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_4.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_4", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_4.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)


## WPRE | MCHERRY >= 5
seuObj$WPRE_GTEQ_5 <- gsub("^1$|^2$|^3$|^4$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_5 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_5)
seuObj$MCHERRY_GTEQ_5 <- gsub("^1$|^2$|^3$|^4$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_5 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_5)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_5", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_5.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_5", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_5.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)



## WPRE | MCHERRY >= 6
seuObj$WPRE_GTEQ_6 <- gsub("^1$|^2$|^3$|^4$|^5$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_6 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_6)
seuObj$MCHERRY_GTEQ_6 <- gsub("^1$|^2$|^3$|^4$|^5$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_6 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_6)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_6", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_6.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_6", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_6.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)



## WPRE | MCHERRY >= 7
seuObj$WPRE_GTEQ_7 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_7 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_7)
seuObj$MCHERRY_GTEQ_7 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_7 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_7)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_7", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_7.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_7", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_7.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)



## WPRE | MCHERRY >= 8
seuObj$WPRE_GTEQ_8 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_8 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_8)
seuObj$MCHERRY_GTEQ_8 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_8 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_8)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_8", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_8.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_8", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_8.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)



## WPRE | MCHERRY >= 9
seuObj$WPRE_GTEQ_9 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$|^8$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_9 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_9)
seuObj$MCHERRY_GTEQ_9 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$|^8$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_9 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_9)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_9", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_9.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_9", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_9.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)



## WPRE | MCHERRY >= 10
seuObj$WPRE_GTEQ_10 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$|^8$|^9$", "0", seuObj$WPRE)
seuObj$WPRE_GTEQ_10 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$WPRE_GTEQ_10)
seuObj$MCHERRY_GTEQ_10 <- gsub("^1$|^2$|^3$|^4$|^5$|^6$|^7$|^8$|^9$", "0", seuObj$MCHERRY)
seuObj$MCHERRY_GTEQ_10 <- gsub("[1-9]|[1-9][0-9]|[1-9][0-9][0-9]", "1", seuObj$MCHERRY_GTEQ_10)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "WPRE_GTEQ_10", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "blue", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_WPRE_GTEQ_10.PDF", plot = pumap1b, width = 6.5, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.05, group.by = "MCHERRY_GTEQ_10", raster = TRUE, label = FALSE, label.size = 2, colors_use = c("grey", "red", "white")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_MCHERRY_GTEQ_10.PDF", plot = pumap2b, width = 6.5, height = 6, units = "in", dpi = 150)


save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY.RData")

```


### Update Annotation
```{R}
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
library(scCustomize)
library(viridis)


load("SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY.RData")
# seuObj

seuObj$Annotation <- gsub("^0$|^1$|^2$|^3$|^4$|^5$|^11$|^17$", "dSPN", seuObj$RNA_snn_res.1.2)
seuObj$Annotation <- gsub("^6$|^7$|^10$|^12$", "iSPN", seuObj$Annotation)
seuObj$Annotation <- gsub("^8$|^13$", "eSPN", seuObj$Annotation)
seuObj$Annotation <- gsub("^9$|^14$|^15$|^16$|^18$|^19$", "Undetermined", seuObj$Annotation)

seuObj$Annotation_Cluster <- paste(seuObj$Annotation, seuObj$RNA_snn_res.1.2, sep = "_")

save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED.RData")


wpreCellType <- as.data.frame.matrix(table(seuObj$Annotation, seuObj$WPRE))
mcherryCellType <- as.data.frame.matrix(table(seuObj$Annotation, seuObj$MCHERRY))

wpreCellTypeCluster <- as.data.frame.matrix(table(seuObj$Annotation_Cluster, seuObj$WPRE))
mcherryCellTypeCluster <- as.data.frame.matrix(table(seuObj$Annotation_Cluster, seuObj$MCHERRY))

wpreGenoVirAgeSample <- as.data.frame.matrix(table(seuObj$GenoVirAgeSample, seuObj$WPRE))
mcherryGenoVirAgeSample <- as.data.frame.matrix(table(seuObj$GenoVirAgeSample, seuObj$MCHERRY))


write.table(wpreCellType, "NK_CELLTYPE_WPRE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(wpreCellType, "NK_CELLTYPE_MCHERRY.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(wpreCellTypeCluster, "NK_CELLTYPE_CLUSTER_WPRE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mcherryCellTypeCluster, "NK_CELLTYPE_CLUSTER_MCHERRY.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(wpreGenoVirAgeSample, "NK_SAMPLE_WPRE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mcherryGenoVirAgeSample, "NK_SAMPLE_MCHERRY.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

```


### Intronic Reads Ratio
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
library(ggpubr)
options(future.globals.maxSize= 50000 * 1024^2)


## Load Seurat object with intronic reads ratio
load("./../SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED_INTRONIC_RATIO.RData")
# seuObj

metaAll <- as.data.frame(seuObj@meta.data)
rm(seuObj)

## Load annotated Seurat object
load("SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED.RData")
# seuObj

metaTemp <- as.data.frame(seuObj@meta.data)

metaData <- merge(metaAll, metaTemp, by = "row.names", all.y = TRUE)
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL

metaIRR <- metaData$INTRONIC_READS_RATIO
names(metaIRR) <- row.names(metaData)

seuObj$INTRONIC_READS_RATIO <- metaIRR

save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED_INTRONIC_RATIO.RData")


## Plots
Idents(seuObj) <- "Annotation"

seuObjSel <- subset(seuObj, idents = "Undetermined", invert = TRUE)

pumap1a <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.01, group.by = "Annotation", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("MAY2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_CELLTYPE_CLEANED.PDF", plot = pumap1a, width = 4, height = 4, units = "in", dpi = 150)

pumap1b <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.01, group.by = "Annotation_Cluster", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("MAY2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_CELLTYPE_CLUSTER_CLEANED.PDF", plot = pumap1b, width = 4, height = 4, units = "in", dpi = 150)

metaNew <- as.data.frame(seuObjSel@meta.data)

p1 <- ggboxplot(metaNew, x = 'Annotation_Cluster', y = 'INTRONIC_READS_RATIO', color = 'Annotation') + rotate_x_text(90)
ggsave("MAY2023_SEURAT_NK_P18_SPN_NuclearFraction_Per_CellType.pdf", p1, width = 6, height = 3, units = "in", dpi = 300)


p1 <- Stacked_VlnPlot(seurat_object = seuObjSel, group.by = "Annotation_Cluster", features = c("nCount_RNA", "nFeature_RNA", "INTRONIC_READS_RATIO"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "MAY2023_SEURAT_NK_P18_SPN_nUMI_nGenes_IRR_Per_CellType.pdf", plot = p1, width = 10, height = 6, units = "in", dpi = 300)

p2 <- Stacked_VlnPlot(seurat_object = seuObjSel, group.by = "Annotation_Cluster", split.by = "GenoVir", features = c("nCount_RNA", "nFeature_RNA", "INTRONIC_READS_RATIO"), colors_use = c("blue2", "green4", "red3"), x_lab_rotate = TRUE, plot_legend = TRUE)
ggsave(filename = "MAY2023_SEURAT_NK_P18_SPN_nUMI_nGenes_IRR_Per_CellType_GenoVir.pdf", plot = p2, width = 15, height = 6, units = "in", dpi = 300)

p3 <- Stacked_VlnPlot(seurat_object = seuObjSel, features = "Foxp1", group.by = "Annotation_Cluster", split.by = "GenoVir", colors_use = c("blue2", "green4", "red3"), x_lab_rotate = TRUE, plot_legend = TRUE)
ggsave(filename = "MAY2023_SEURAT_NK_P18_SPN_nUMI_nGenes_IRR_Per_CellType_GenoVir_Foxp1.pdf", plot = p3, width = 15, height = 3, units = "in", dpi = 300)


```


### ADDITIONAL PLOTS
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
library(ggpubr)
library(paletteer)
options(future.globals.maxSize= 50000 * 1024^2)


##------------------------------------
load("./../SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED_INTRONIC_RATIO.RData")
# seuObj

##------------------------------------
## Plots
Idents(seuObj) <- "Annotation"

seuObjSel <- subset(seuObj, idents = "Undetermined", invert = TRUE)

seuObjSel$AnnotationA <- gsub("dSPN", "D1-SPN", seuObjSel$Annotation)
seuObjSel$AnnotationA <- gsub("iSPN", "D2-SPN", seuObjSel$AnnotationA)
seuObjSel$AnnotationA <- gsub("eSPN", "eSPN", seuObjSel$AnnotationA)

seuObjSel$AnnotationB <- gsub("^dSPN_0$", "D1-00", seuObjSel$Annotation_Cluster)
seuObjSel$AnnotationB <- gsub("^dSPN_1$", "D1-01", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_2$", "D1-02", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_3$", "D1-03", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_4$", "D1-04", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_5$", "D1-05", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_11$", "D1-11", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^dSPN_17$", "D1-17", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^eSPN_8$", "eSPN-08", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^eSPN_13$", "eSPN-13", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^iSPN_6$", "D2-06", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^iSPN_7$", "D2-07", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^iSPN_10$", "D2-10", seuObjSel$AnnotationB)
seuObjSel$AnnotationB <- gsub("^iSPN_12$", "D2-12", seuObjSel$AnnotationB)


mycol1 <- c("red", "blue", "green")
pumap1a <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "AnnotationA", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol1) + NoLegend()
ggsave("JUN2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_CELLTYPE_CLEANED.PDF", plot = pumap1a, width = 4, height = 4, units = "in", dpi = 300)


mycol2 <- c(paletteer_c("grDevices::Reds 3", 12)[1:8], paletteer_c("grDevices::Blues 3", 7)[1:4], paletteer_c("grDevices::Greens 3", 4)[1:2])
pumap1b <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "AnnotationB", raster = TRUE, raster.dpi = c(1024, 1024), label = TRUE, label.size = 2, colors_use = mycol2) + NoLegend()
ggsave("JUN2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_CELLTYPE_CLUSTER_CLEANED.PDF", plot = pumap1b, width = 4, height = 4, units = "in", dpi = 300)

pumap1c <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "AnnotationB", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, label.size = 2, colors_use = mycol2) + NoLegend()
ggsave("JUN2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_CELLTYPE_CLUSTER_CLEANED2.PDF", plot = pumap1c, width = 4, height = 4, units = "in", dpi = 300)


seuObjSel$GenoVirA <- gsub("WT_CTL", "A_CTL", seuObjSel$GenoVir)
seuObjSel$GenoVirA <- gsub("CKO_CTL", "B_CKO", seuObjSel$GenoVirA)
seuObjSel$GenoVirA <- gsub("CKO_RES", "C_RES", seuObjSel$GenoVirA)

mycol3 <- c("#646464", "#FA7D64", "#0A7D7D")
pumap1d <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "GenoVirA", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, label.size = 2, colors_use = mycol3) + NoLegend()
ggsave("JUN2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_GENOVIR_CLEANED.PDF", plot = pumap1d, width = 4, height = 4, units = "in", dpi = 300)

pumap1e <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "GenoVirA", split.by = "GenoVirA", split_seurat = TRUE, raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, label.size = 2, colors_use = mycol3) + NoLegend()
ggsave("JUN2023_SEURAT_NK_P18_CKO_RES_HARMONY_SPN_UMAP_GENOVIR_CLEANED2.PDF", plot = pumap1e, width = 10, height = 4, units = "in", dpi = 300)

```


