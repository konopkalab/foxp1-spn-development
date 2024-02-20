# CLUSTERING WITH SEURAT AND HARMONY

## Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

## Harmony
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

### Annotate Clusters
```{R}

##------------------------------------
## Annotate Clusters
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(WGCNA)


cluMarkers <- read.table("SEURAT_NK_P18_CKO_RES_HARMONY_RES_1.2_CLUSTER_MARKERS.txt", sep = "\t", header = TRUE)
tab <- cluMarkers
# tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.75 & tab$pct.1 >= 0.75,]
tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.5,]
# tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.25,]
tab <- tab[c(7,6)]
tab$cluster <- as.factor(paste("Cluster_", sprintf("%02d", as.numeric(as.character(tab$cluster))), sep=""))
colnames(tab)=c("Gene","DEFINITION")
Genes=as.data.frame(table(tab$DEFINITION))


load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/A_JAN2022/05_SEURAT/MISC/Saunders2018DEG.RData")
# degSaunders
GeneSets <- degSaunders


for(i in 1:length(GeneSets)){
	colnames(GeneSets[[i]])[1] <- "Gene"
}

ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- length(unique(sort(cluMarkers$gene)))-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585 #13517 #nrow(tab) #8321 # length(unique(sort(cluMarkers$gene)))(genes in seurat object) 
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
	f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
	return(c(row,
			P_val = f$p.value,
			LogP = -log10(f$p.value), 
			OR = f$estimate[[1]],
			OR_Low = f$conf.int[1],
			OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat, TEMP, file= "NK_SEURAT_RES_1.2_FisherOutput_Saunders_Enrich.RData")


# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0


pdf("NK_SEURAT_RES_1.2_FisherOutput_Saunders_Enrich.pdf", width=12, height=16, pointsize=12)
par(mar=c(15, 7, 2, 2))
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5,
xLabelsAngle = 90)
dev.off()



FisherORt <- as.data.frame(t(FisherOR))
colnames(FisherORt) <- paste(colnames(FisherORt), "OR", sep = "_")

FisherAdjt <- as.data.frame(t(FisherAdj))
colnames(FisherAdjt) <- paste(colnames(FisherAdjt), "Pval", sep = "_")

FisherData <- merge(FisherORt, FisherAdjt, by = "row.names")
row.names(FisherData) <- FisherData$Row.names
FisherData$Row.names <- NULL

FisherData2 <- FisherData[,order(colnames(FisherData))]
write.table(FisherData2, "NK_SEURAT_RES_1.2_FisherOutput_Saunders_PlotData.txt", row.names = T, col.names = T, quote = F, sep = "\t")

```

### Update Seurat Object
```{R}

##------------------------------------
## Update Seurat Object
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(WGCNA)
library(scCustomize)

load("SEURAT_NK_P18_CKO_RES_HARMONY.RData")
# seuObj

Idents(seuObj) <- "integrated_snn_res.1.2"

metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp$CellBarcode <- row.names(metaTemp)

annotations <- read.table("NK_SEURAT_ANNOTATION_PER_CLUSTER_1.2.txt", header = TRUE, sep = "\t")

annotations2 <- merge(metaTemp, annotations, by = "RNA_snn_res.1.2")

annotations3 <- annotations2$CellType
names(annotations3) <- annotations2$CellBarcode

seuObj[["CellType"]] <- annotations3
seuObj[["CellType_Cluster"]] <- paste(seuObj@meta.data$CellType, seuObj@meta.data$RNA_snn_res.1.2, sep = "_")


pumap1a <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType", raster = TRUE, label = TRUE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$CellType))), palette = "varibow"), repel = TRUE) + NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_CELLTYPE.PDF", plot = pumap1a, width = 6, height = 6, units = "in", dpi = 150)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType", raster = TRUE, label = FALSE, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$CellType))), palette = "varibow")) + facet_wrap(~CellType, nrow = 4)
ggsave("NK_SEURAT_PLOT_UMAP_CELLTYPE_FACET.PDF", plot = pumap1b, width = 12, height = 12, units = "in", dpi = 150)

pumap1c <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType", raster = TRUE, label = FALSE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$CellType))), palette = "varibow")) #+ NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_CELLTYPE_NOLABEL.PDF", plot = pumap1c, width = 8, height = 6, units = "in", dpi = 150)

pumap5a <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType_Cluster", raster = TRUE, label = TRUE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$CellType_Cluster))), palette = "varibow"), repel = TRUE) + NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_CELLTYPE_CLUSTER.PDF", plot = pumap5a, width = 6, height = 6, units = "in", dpi = 150)

pumap5b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType_Cluster", raster = TRUE, label = FALSE, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(seuObj$CellType_Cluster))), palette = "varibow")) + facet_wrap(~CellType_Cluster, nrow = 5) + NoLegend()
ggsave("NK_SEURAT_PLOT_UMAP_CELLTYPE_CLUSTER_FACET.PDF", plot = pumap5b, width = 28, height = 16, units = "in", dpi = 150)

p6 <- Stacked_VlnPlot(seurat_object = seuObj, group.by = "CellType", features = c("nCount_RNA", "nFeature_RNA"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "NK_SEURAT_CELLTYPE_VIOLINPLOT_nUMI.pdf", plot = p6, width = 12, height = 8, units = "in", dpi = 300)


save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED.RData")

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

## load Seurat Harmony object
load("SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED.RData")
# seuObj


## WPRE-MCHERRY
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



## WPRE-MCHERRY by sample
seumeta <- as.data.frame(seuObj@meta.data)

df <- seumeta[,c("GenoVirAgeSample", "CellType_Cluster", "WPRE", "MCHERRY")]


for(myclu in unique(sort(df$CellType_Cluster)))
    {
    print(myclu)
    seumeta_sub <- seumeta[seumeta$CellType_Cluster == myclu,]
    print(nrow(seumeta_sub))
    wpresub <- as.data.frame.matrix(table(seumeta_sub$GenoVirAgeSample, seumeta_sub$WPRE))
    mcherrysub <- as.data.frame.matrix(table(seumeta_sub$GenoVirAgeSample, seumeta_sub$MCHERRY))
    write.table(wpresub, paste("NK_ALL_WPRE_BY_SAMPLE_BY_", myclu, ".txt", sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
    write.table(mcherrysub, paste("NK_ALL_MCHERRY_BY_SAMPLE_BY_", myclu, ".txt", sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
    }


## Updated Meta Data
seuMeta <- as.data.frame(seuObj@meta.data)
seuMeta$CellBarcode <- row.names(seuMeta)

annotation <- read.table("NK_SEURAT_ANNOTATION_PER_CLUSTER_1.2_2.txt", header = TRUE, sep = "\t")
seuMetaUpdated <- merge(seuMeta, annotation, by = "RNA_snn_res.1.2")

celltype2 <- seuMetaUpdated$CellType2
names(celltype2) <- seuMetaUpdated$CellBarcode

seuObj$CellType2 <- celltype2
seuObj$CellType2_Cluster <- paste(seuObj$CellType2, seuObj$RNA_snn_res.1.2, sep = "_")

pumap1a <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType2", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("SEURAT_NK_P18_CKO_RES_HARMONY_UMAP_CELLTYPE.PDF", plot = pumap1a, width = 6, height = 6, units = "in", dpi = 150)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj, pt.size = 0.01, group.by = "CellType2_Cluster", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("SEURAT_NK_P18_CKO_RES_HARMONY_UMAP_CELLTYPE_CLUSTER.PDF", plot = pumap1b, width = 6, height = 6, units = "in", dpi = 150)


Idents(seuObj) <- "CellType2"

seuObjSel <- subset(seuObj, idents = "Undetermined", invert = TRUE)

pumap2a <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.01, group.by = "CellType2", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("SEURAT_NK_P18_CKO_RES_HARMONY_UMAP_CELLTYPE_CLEANED.PDF", plot = pumap2a, width = 6, height = 6, units = "in", dpi = 150)

pumap2b <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.01, group.by = "CellType2_Cluster", raster = TRUE, label = TRUE, label.size = 2, ggplot_default_colors = TRUE) + NoLegend()
ggsave("SEURAT_NK_P18_CKO_RES_HARMONY_UMAP_CELLTYPE_CLUSTER_CLEANED.PDF", plot = pumap2b, width = 6, height = 6, units = "in", dpi = 150)


## Save Seurat Object
save(seuObj, seuObjSel, file = "SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED_ANNOTATED_CLEANED.RData")


```


### Intronic Reads Ratio
```{R}
## Load annotated seurat object
load("SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED.RData")
# seuObj

as.data.frame(table(seuObj$GenoVirAgeSample))
#               Var1  Freq
# 1 CKO_CTL_P18_3BF1 10834
# 2 CKO_CTL_P18_3BM3 11390
# 3 CKO_CTL_P18_3BM4  7243
# 4 CKO_RES_P18_3BM8 20586
# 5 CKO_RES_P18_3CF5 10169
# 6 CKO_RES_P18_3CF8 10073
# 7    WT_CTL_P18_C1  9507
# 8    WT_CTL_P18_C2 10996
# 9    WT_CTL_P18_C3  9206

metaTemp <- as.data.frame(seuObj@meta.data)

##------------------------------------
## load intronic reads ratios
load("NK_INTRONIC_READS_RATIO_P18.RData")
# WT_CTL_P18_C1_IR, WT_CTL_P18_C2_IR, WT_CTL_P18_C3_IR, CKO_CTL_P18_3BF1_IR, CKO_CTL_P18_3BM3_IR, CKO_CTL_P18_3BM4_IR, CKO_RES_P18_3BM8_IR, CKO_RES_P18_3CF5_IR, CKO_RES_P18_3CF8_IR

row.names(WT_CTL_P18_C1_IR) <- paste("P18_CTL_C1_", row.names(WT_CTL_P18_C1_IR), sep = "")
row.names(WT_CTL_P18_C2_IR) <- paste("P18_CTL_C2_", row.names(WT_CTL_P18_C2_IR), sep = "")
row.names(WT_CTL_P18_C3_IR) <- paste("P18_CTL_C3_", row.names(WT_CTL_P18_C3_IR), sep = "")
row.names(CKO_CTL_P18_3BF1_IR) <- paste("P18_CTL_3BF1_", row.names(CKO_CTL_P18_3BF1_IR), sep = "")
row.names(CKO_CTL_P18_3BM3_IR) <- paste("P18_CTL_3BM3_", row.names(CKO_CTL_P18_3BM3_IR), sep = "")
row.names(CKO_CTL_P18_3BM4_IR) <- paste("P18_CTL_3BM4_", row.names(CKO_CTL_P18_3BM4_IR), sep = "")
row.names(CKO_RES_P18_3BM8_IR) <- paste("P18_RES_3BM8_", row.names(CKO_RES_P18_3BM8_IR), sep = "")
row.names(CKO_RES_P18_3CF5_IR) <- paste("P18_RES_3CF5_", row.names(CKO_RES_P18_3CF5_IR), sep = "")
row.names(CKO_RES_P18_3CF8_IR) <- paste("P18_RES_3CF8_", row.names(CKO_RES_P18_3CF8_IR), sep = "")

intronic_ratio <- rbind(WT_CTL_P18_C1_IR, WT_CTL_P18_C2_IR, WT_CTL_P18_C3_IR, CKO_CTL_P18_3BF1_IR, CKO_CTL_P18_3BM3_IR, CKO_CTL_P18_3BM4_IR, CKO_RES_P18_3BM8_IR, CKO_RES_P18_3CF5_IR, CKO_RES_P18_3CF8_IR)

metaTemp2 <- merge(metaTemp, intronic_ratio, by = "row.names")
row.names(metaTemp2) <- metaTemp2$Row.names
metaTemp2$Row.names <- NULL

# Plot nuclear fraction per cell type
p1 <- ggboxplot(metaTemp2, x = 'CellType_Cluster', y = 'nuclear_fraction', color = 'CellType') + rotate_x_text(90)
ggsave("NuclearFraction_Per_CellType.pdf", p1, width = 20, height = 10, units = "in", dpi = 300)

pdf('NuclearFraction_Per_Cluster.pdf')
ggscatter(metaTemp2, x = 'CellType_Cluster', y = 'nuclear_fraction')
dev.off()

metaIRR <- metaTemp2$nuclear_fraction
names(metaIRR) <- row.names(metaTemp2)

seuObj$INTRONIC_READS_RATIO <- metaIRR

save(seuObj, file = "SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED_INTRONIC_RATIO.RData")

```
