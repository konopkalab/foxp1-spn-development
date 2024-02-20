# CLUSTERING WITH SEURAT AND HARMONY

## Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

## Run Harmony
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
## Load seurat object for 
## Individual genotype
#### WT-CTL
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/D_SEP2022/05_SEURAT_AFTER_CB_DF/A_P18_CONTROL/NK_P18_CTL_CB_DF_CLUST.RData")
wt.ctl.seurat <- seuObj
rm(seuObj)

## CKO-CTL
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/08_SEURAT_AFTER_CB_DF/A_P18_CONTROL/NK_P18_CTL_CB_DF_CLUST.RData")
cko.ctl.seurat.temp <- seuObj
rm(seuObj)

## CKO-RES
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/08_SEURAT_AFTER_CB_DF/B_P18_RESCUE/NK_P18_RES_CB_DF_CLUST.RData")
cko.res.seurat.temp <- seuObj
rm(seuObj)


##------------------------------------
## remove outlier samples
cko.ctl.seurat <- subset(cko.ctl.seurat.temp, subset = GenoAgeSample %in% c("CTL_P18_3CF1"), invert = TRUE)
table(cko.ctl.seurat$GenoAgeSample)
# CTL_P18_3BF1 CTL_P18_3BM3 CTL_P18_3BM4 
#        10834        11390         7243

cko.res.seurat <- subset(cko.res.seurat.temp, subset = GenoAgeSample %in% c("RES_P18_3BF5"), invert = TRUE)
table(cko.res.seurat$GenoAgeSample)
# RES_P18_3BM8 RES_P18_3CF5 RES_P18_3CF8 
#        20586        10169        10073

##------------------------------------
## re-organize seurat objects
wt.ctl.seurat$Virus <- wt.ctl.seurat$Genotype
wt.ctl.seurat$Genotype <- rep("WT", nrow(wt.ctl.seurat@meta.data))
wt.ctl.seurat$GenoVirAgeSample <- paste(wt.ctl.seurat$Genotype, wt.ctl.seurat$Virus, wt.ctl.seurat$Age, wt.ctl.seurat$Sample, sep = "_")

table(wt.ctl.seurat@meta.data$GenoVirAgeSample)
# WT_CTL_P18_C1 WT_CTL_P18_C2 WT_CTL_P18_C3 
#          9507         10996          9206


cko.ctl.seurat$Virus <- cko.ctl.seurat$Genotype
cko.ctl.seurat$Genotype <- rep("CKO", nrow(cko.ctl.seurat@meta.data))
cko.ctl.seurat$GenoVirAgeSample <- paste(cko.ctl.seurat$Genotype, cko.ctl.seurat$Virus, cko.ctl.seurat$Age, cko.ctl.seurat$Sample, sep = "_")

table(cko.ctl.seurat@meta.data$GenoVirAgeSample)
# CKO_CTL_P18_3BF1 CKO_CTL_P18_3BM3 CKO_CTL_P18_3BM4 
#            10834            11390             7243

cko.res.seurat$Virus <- cko.res.seurat$Genotype
cko.res.seurat$Genotype <- rep("CKO", nrow(cko.res.seurat@meta.data))
cko.res.seurat$GenoVirAgeSample <- paste(cko.res.seurat$Genotype, cko.res.seurat$Virus, cko.res.seurat$Age, cko.res.seurat$Sample, sep = "_")

table(cko.res.seurat@meta.data$GenoVirAgeSample)
# CKO_RES_P18_3BM8 CKO_RES_P18_3CF5 CKO_RES_P18_3CF8 
#            20586            10169            10073


##------------------------------------
## add additional meta information
wt.ctl.seurat$Sex <- gsub("^WT_CTL_P18_C1$", "F", wt.ctl.seurat$GenoVirAgeSample)
wt.ctl.seurat$Sex <- gsub("^WT_CTL_P18_C2$|^WT_CTL_P18_C3$", "M", wt.ctl.seurat$Sex)

cko.ctl.seurat$Sex <- gsub("^CKO_CTL_P18_3BF1$", "F", cko.ctl.seurat$GenoVirAgeSample)
cko.ctl.seurat$Sex <- gsub("^CKO_CTL_P18_3BM3$|^CKO_CTL_P18_3BM4$", "M", cko.ctl.seurat$Sex)

cko.res.seurat$Sex <- gsub("^CKO_RES_P18_3CF5$", "F", cko.res.seurat$GenoVirAgeSample)
cko.res.seurat$Sex <- gsub("^CKO_RES_P18_3BM8$|^CKO_RES_P18_3CF8$", "M", cko.res.seurat$Sex)

wt.ctl.seurat$Batch <- rep(2, nrow(wt.ctl.seurat@meta.data))
cko.ctl.seurat$Batch <- rep(1, nrow(cko.ctl.seurat@meta.data))
cko.res.seurat$Batch <- rep(1, nrow(cko.res.seurat@meta.data))


##------------------------------------
## merge seurat object for genotypes
seuObj <- merge(x = wt.ctl.seurat, y = c(cko.ctl.seurat, cko.res.seurat))

table(seuObj$GenoVirAgeSample)
# CKO_CTL_P18_3BF1 CKO_CTL_P18_3BM3 CKO_CTL_P18_3BM4 CKO_RES_P18_3BM8 
#            10834            11390             7243            20586 
# CKO_RES_P18_3CF5 CKO_RES_P18_3CF8    WT_CTL_P18_C1    WT_CTL_P18_C2 
#            10169            10073             9507            10996 
#    WT_CTL_P18_C3 
#             9206

##------------------------------------
## preprocess
seuObj <- NormalizeData(seuObj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE, npcs = 100)


##-------------------------------------------------------
## jackstraw analysis
## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)

ggsave(filename = "SEURAT_NK_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_NK_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)


##------------------------------------
## harmony
seuObj <- RunHarmony(seuObj, group.by.vars = c("Batch", "Sex"), reduction = "pca", assay.use = "RNA", max.iter.harmony = 100)


##------------------------------------
## clustering
seuObj <- FindNeighbors(seuObj, dims = 1:selpcs, reduction = "harmony")
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), reduction = "harmony")
seuObj <- RunUMAP(seuObj, dims = 1:selpcs, reduction = "harmony")


##-------------------------------------------------------
## SAVE RDATA
save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY.RData")



##------------------------------------
## Plots

## CLUSTREE
seutree <- clustree(seuObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SEURAT_NK_INTEGRATE_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## UMAP FOR ALL RESOLUTIONS
ResolutionList <- grep("RNA_snn_res", colnames(seuObj@meta.data), value = TRUE)

for (Resolution in ResolutionList)
    {
    print(paste("====> RESOLUTION ", Resolution, sep = ""))

    pdf(paste0("SEURAT_NK_INTEGRATE_UMAP_RES_", Resolution, ".pdf"), width=7, height=6)
    g <- DimPlot(object = seuObj, label = TRUE, reduction = "umap", group.by = Resolution)
    print(g)
    dev.off()

    pdf(paste0("SEURAT_NK_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width=7, height=3)
    v <- VlnPlot(object = seuObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
    print(v)
    dev.off()
    }


## UMAP FOR MARKER GENES
mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cspg4", "Mag", "Cx3cr1", "Flt1", "Slc17a7", "Chat", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b", "Oprm1", "Isl1", "Pdyn", "Lypd1", "Nnat", "Ebf1", "Epha4", "Mef2c")

fpl1 <- FeaturePlot_scCustom(seurat_object = seuObj, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 15, height = 24, units = "in", dpi = 150)

fpl2 <- FeaturePlot_scCustom(seurat_object = seuObj, features = mygenes, reduction = "umap", pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 15, height = 24, units = "in", dpi = 150)


plotCluUMAP1 <- DimPlot_scCustom(seurat_object = seuObj, group.by = "Genotype", pt.size = 0.1, reduction = "umap", label = FALSE, raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOTYPE.pdf", plot = plotCluUMAP1, width = 10, height = 8, units = "in", dpi = 300)

plotCluUMAP1b <- DimPlot_scCustom(seurat_object = seuObj, group.by = "Genotype", split.by = "Genotype", num_columns = 4, pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOTYPE_FACET.pdf", plot = plotCluUMAP1b, width = 32, height = 8, units = "in", dpi = 300)

plotCluUMAP1c <- DimPlot_scCustom(seurat_object = seuObj, group.by = "GenoVirAgeSample", split.by = "GenoVirAgeSample", num_columns = 3, pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOVIRAGESAMPLE_FACET.pdf", plot = plotCluUMAP1c, width = 25, height = 20, units = "in", dpi = 300)

plotCluUMAP1d <- DimPlot_scCustom(seurat_object = seuObj, group.by = "GenoVirAgeSample", pt.size = 0.1, reduction = "umap", raster = TRUE)
ggsave(filename = "NK_SEURAT_UMAP_GENOVIRAGESAMPLE.pdf", plot = plotCluUMAP1d, width = 8, height = 6, units = "in", dpi = 300)

plotCluUMAP1e <- DimPlot(seuObj, group.by = c("Sex", "Batch", "Genotype", "Virus", "GenoVirAgeSample"), pt.size = 0.01)
ggsave("NK_SEURAT_UMAP_COVARIATES.pdf", plot = plotCluUMAP1e, width = 20, height = 16, units = "in", dpi = 300)



# RESOLUTION 1.2

my_levels <- seq(0, max(as.numeric(seuObj$RNA_snn_res.1.2)) - 1)
seuObj$RNA_snn_res.1.2 <- factor(x = seuObj$RNA_snn_res.1.2, levels = my_levels)

vpl2 <- Stacked_VlnPlot(seurat_object = seuObj, features = mygenes, group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_ViolinPlot_1.2_A.pdf", plot = vpl2 , width = 12, height = 18, units = "in", dpi = 150, useDingbats=FALSE)

vpl3 <- Stacked_VlnPlot(seurat_object = seuObj, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.1.2", pt.size = 0, x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_ViolinPlot_1.2_B.pdf", plot = vpl3 , width = 12, height = 3, units = "in", dpi = 150, useDingbats=FALSE)

dpl2 <- DotPlot_scCustom(seurat_object = seuObj, features = rev(mygenes), colors_use = viridis_plasma_dark_high, group.by = "RNA_snn_res.1.2", x_lab_rotate = TRUE)
ggsave(filename = "NK_SEURAT_PLOT_DotPlot_1.2.pdf", plot = dpl2, width = 12, height = 16, units = "in", dpi = 150, useDingbats=FALSE)


## RES 1.2
## BARPLOT SHOWING CELLS PER CLUSTER PER GENOTYPE | RES 1.2
cellsPerCluster <- as.data.frame.matrix(table(seuObj@meta.data$RNA_snn_res.1.2, seuObj@meta.data$Genotype))
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
cellsPerCluster <- as.data.frame.matrix(table(seuObj@meta.data$RNA_snn_res.1.2, seuObj@meta.data$GenoVirAgeSample))
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
library(harmony)
library(scCustomize)
options(future.globals.maxSize= 50000 * 1024^2)

## Load Seurat Object
load("SEURAT_NK_P18_CKO_RES_HARMONY.RData")

Idents(seuObj) <- "RNA_snn_res.1.2"

seuObj.markers <- FindAllMarkers(seuObj, only.pos = TRUE) #, min.pct = 0, logfc.threshold = 0)
write.table(seuObj.markers, "SEURAT_NK_P18_CKO_RES_HARMONY_RES_1.2_CLUSTER_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
save(seuObj.markers, file = "SEURAT_NK_P18_CKO_RES_HARMONY_RES_1.2_CLUSTER_MARKERS.RData")

```
