# FIGURES FOR MANUSCRIPT

### Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

### All Cells
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

## Load cleaned, annotated and updated seurat object
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SEURAT_NK_P18_CKO_RES_HARMONY_UPDATED_INTRONIC_RATIO.RData")
# seuObj

## UPDATE CELL-TYPE LABELS
seuObj2 <- subset(seuObj, subset = CellType %in% c("Cortical", "Undetermined"), invert = TRUE)
seuObj2$Annotation <- gsub("^0$|^1$|^3$|^32$|^35$", "D1_SPNs", seuObj2$RNA_snn_res.1.2)
seuObj2$Annotation <- gsub("^4$|^9$", "D2_SPNs", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^10$|^18$|^25$", "eSPNs", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^2$|^12$|^14$|^21$", "Astrocytes", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^17$|^24$", "Endothelial_Cells", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^31$", "Interneurons", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^7$|^34$", "Microglia", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^22$|^41$", "Progenitors", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^8$|^23$", "Neurogenic_Progenitors", seuObj2$Annotation)
seuObj2$Annotation <- gsub("^5$|^6$|^11$|^13$|^16$|^19$|^30$|^33$|^38$", "Oligodendrocytes", seuObj2$Annotation)

## HEX VALUES FOR CELLTYPE COLORS
astrocytes <- "#077E97" # <--
d1spns <- "#FF0000" # <--
d2spns <- "#008000" # <--
endothelialcells <- "#C0C000" # <--
espns <- "#0000C0" # <--
interneurons <- "#808080" # <--
microglia <- "#F2B77C" # <--
neurogenic <- "#FFA040" # <--
oligodendrocytes <- "#F71480" # <--
progenitors <- "#800080" # <--

## COLOR PALETTE FOR CELLTYPE
mycolAll <- c(astrocytes, d1spns, d2spns, endothelialcells, espns, interneurons, microglia, neurogenic, oligodendrocytes, progenitors)

## UMAP PLOTS FOR CELLTYPE
pumap1a <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycolAll) #+ NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_A.PDF", plot = pumap1a, width = 6, height = 4, units = "in", dpi = 300)

pumap1b <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycolAll) + NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_B.PDF", plot = pumap1b, width = 4, height = 4, units = "in", dpi = 300)

pumap1c <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation", raster = TRUE, raster.dpi = c(1024, 1024), label = TRUE, label.size = 2, repel = TRUE, colors_use = mycolAll) + NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_C.PDF", plot = pumap1c, width = 4, height = 4, units = "in", dpi = 300)


## COLOR PALETTE FOR CELLTYPE CLUSTERS
seuObj2$Annotation_Cluster <- paste(seuObj2$Annotation, seuObj2$RNA_snn_res.1.2, sep = "_")
mycolAll2 <- c(astrocytes, astrocytes, astrocytes, astrocytes, d1spns, d1spns, d1spns, d1spns, d1spns, d2spns, d2spns, endothelialcells, endothelialcells, espns, espns, espns, interneurons, microglia, microglia, neurogenic, neurogenic, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, oligodendrocytes, progenitors, progenitors)

## UMAP PLOTS FOR CELLTYPE CLUSTERS
pumap1d <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation_Cluster", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycolAll2) #+ NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_D.PDF", plot = pumap1d, width = 12, height = 6, units = "in", dpi = 300)

pumap1e <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation_Cluster", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycolAll2) + NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_E.PDF", plot = pumap1e, width = 4, height = 4, units = "in", dpi = 300)

pumap1f <- DimPlot_scCustom(seurat_object = seuObj2, pt.size = 0.1, group.by = "Annotation_Cluster", raster = TRUE, raster.dpi = c(1024, 1024), label = TRUE, label.size = 2, repel = TRUE, colors_use = mycolAll2) + NoLegend()
ggsave("JUL2023_UMAP_ALL_CELLS_F.PDF", plot = pumap1f, width = 6, height = 6, units = "in", dpi = 300)


## VIOLIN PLOTS FOR CELLTYPE
vpl1 <- VlnPlot(object = seuObj2, group.by = "Annotation", features = c("nCount_RNA", "nFeature_RNA"), cols = mycolAll, pt.size = 0, ncol = 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_ALL_CELLS_A.pdf", plot = vpl1, width = 10, height = 8, units = "in", dpi = 300)

vpl2 <- VlnPlot(object = seuObj2, group.by = "Annotation", features = c("INTRONIC_READS_RATIO"), cols = mycolAll, pt.size = 0) + ylim(0, 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_ALL_CELLS_B.pdf", plot = vpl2, width = 10, height = 4, units = "in", dpi = 300)

## VIOLIN PLOTS FOR CELLTYPE CLUSTERS
vpl3 <- VlnPlot(object = seuObj2, group.by = "Annotation_Cluster", features = c("nCount_RNA", "nFeature_RNA"), cols = mycolAll2, pt.size = 0, ncol = 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_ALL_CELLS_C.pdf", plot = vpl3, width = 15, height = 8, units = "in", dpi = 300)

vpl4 <- VlnPlot(object = seuObj2, group.by = "Annotation_Cluster", features = c("INTRONIC_READS_RATIO"), cols = mycolAll2, pt.size = 0) + ylim(0, 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_ALL_CELLS_D.pdf", plot = vpl4, width = 15, height = 4, units = "in", dpi = 300)

## SAVE RDATA
save(seuObj2, mycolAll, mycolAll2, file = "JUL2023_SEURAT_ALL_CELLS.RData")



## DotPlot for cell-type markers
rm(list = ls())
load("JUL2023_SEURAT_ALL_CELLS.RData")
# seuObj2

seuObjSel <- seuObj2

seuObjSel$Annotation2 <- gsub("^Progenitors$", "A_Progenitors", seuObjSel$Annotation)
seuObjSel$Annotation2 <- gsub("^Neurogenic_Progenitors$", "B_Neurogenic_Progenitors", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^D1_SPNs$", "C_D1_SPNs", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^D2_SPNs$", "D_D2_SPNs", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^eSPNs$", "E_eSPNs", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^Astrocytes$", "F_Astrocytes", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^Oligodendrocytes$", "G_Oligodendrocytes", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^Microglia$", "H_Microglia", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^Endothelial_Cells$", "I_Endothelial_Cells", seuObjSel$Annotation2)
seuObjSel$Annotation2 <- gsub("^Interneurons$", "J_Interneurons", seuObjSel$Annotation2)


DefaultAssay(seuObjSel) <- "RNA"
seuObjSel <- NormalizeData(seuObjSel)

mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")

dpl1 <- DotPlot_scCustom(seurat_object = seuObjSel, features = mygenes, x_lab_rotate = TRUE, group.by = "Annotation2", col.min = 0, col.max = 3, dot.min = 0, dot.scale = 6)
ggsave(filename = "JUL2023_DOTPLOT_ALL_MARKER_GENES.pdf", plot = dpl1, width = 10, height = 4, units = "in", dpi = 300)

```


### SPNs
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
library(EnhancedVolcano)
library(pheatmap)
options(future.globals.maxSize= 50000 * 1024^2)

## Load SPN SUBCLUSTERING seurat object
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED_INTRONIC_RATIO.RData")
# seuObj

## UPDATE CELLTYPE LABELS
Idents(seuObj) <- "Annotation"
seuObjSel <- subset(seuObj, idents = "Undetermined", invert = TRUE)

seuObjSel$AnnotationA <- gsub("dSPN", "D1-SPNs", seuObjSel$Annotation)
seuObjSel$AnnotationA <- gsub("iSPN", "D2-SPNs", seuObjSel$AnnotationA)
seuObjSel$AnnotationA <- gsub("eSPN", "eSPNs", seuObjSel$AnnotationA)

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

seuObjSel$GenoVirA <- gsub("WT_CTL", "A_CTL", seuObjSel$GenoVir)
seuObjSel$GenoVirA <- gsub("CKO_CTL", "B_CKO", seuObjSel$GenoVirA)
seuObjSel$GenoVirA <- gsub("CKO_RES", "C_RES", seuObjSel$GenoVirA)

## HEX VALUES FOR CELLTYPE
d1spns <- "#FF0000" # <--
d2spns <- "#008000" # <--
espns <- "#0000C0" # <--

## HEX VALUES FOR GENOTYPES
ctl <- "#646464" # <--
cko <- "#FF6600" # <--
res <- "#02D4D4" # <--

## HEX VALUES FOR CELLTYPE CLUSTERS
d100 <- "#FF6464"
d101 <- "#FF2500"
d102 <- "#E12901"
d103 <- "#D31A38"
d104 <- "#800000"
d105 <- "#C04000"
d111 <- "#ED9121"
d117 <- "#EE5921"
d206 <- "#299617"
d207 <- "#00FA9A"
d210 <- "#29AB87"
d212 <- "#32CD32"
e08 <- "#0066FF"
e13 <- "#00BBEE"


## COLOR PALETTE FOR CELLTYPE
mycol1 <- c(d1spns, d2spns, espns)

## UMAP PLOTS FOR CELLTYPE
pumap1a <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationA", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol1) #+ NoLegend()
ggsave("JUL2023_UMAP_SPN_A.PDF", plot = pumap1a, width = 5, height = 4, units = "in", dpi = 300)

pumap1b <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationA", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol1) + NoLegend()
ggsave("JUL2023_UMAP_SPN_B.PDF", plot = pumap1b, width = 4, height = 4, units = "in", dpi = 300)

pumap1c <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationA", raster = TRUE, raster.dpi = c(1024, 1024), label = TRUE, label.size = 2, repel = TRUE, colors_use = mycol1) + NoLegend()
ggsave("JUL2023_UMAP_SPN_C.PDF", plot = pumap1c, width = 4, height = 4, units = "in", dpi = 300)

## VIOLIN PLOTS FOR CELLTYPE
vpl1 <- VlnPlot(object = seuObjSel, group.by = "AnnotationA", features = c("nCount_RNA", "nFeature_RNA"), cols = mycol1, pt.size = 0, ncol = 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_A.pdf", plot = vpl1, width = 4, height = 6, units = "in", dpi = 300)

vpl2 <- VlnPlot(object = seuObjSel, group.by = "AnnotationA", features = c("INTRONIC_READS_RATIO"), cols = mycol1, pt.size = 0) + ylim(0, 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_B.pdf", plot = vpl2, width = 4, height = 3, units = "in", dpi = 300)


## COLOR PALETTE FOR GENOTYPES
mycol2 <- c("#646464", "#FF6600", "#02D4D4")

## UMAP PLOTS FOR GENOTYPE
pumap1d <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "GenoVirA", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol2) #+ NoLegend()
ggsave("JUL2023_UMAP_SPN_D.PDF", plot = pumap1d, width = 5, height = 4, units = "in", dpi = 300)

pumap1e <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "GenoVirA", split.by = "GenoVirA", split_seurat = TRUE, raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol2) + NoLegend()
ggsave("JUL2023_UMAP_SPN_E.PDF", plot = pumap1e, width = 10, height = 4, units = "in", dpi = 300)

## VIOLIN PLOTS FOR GENOTYPE
vpl3 <- VlnPlot(object = seuObjSel, group.by = "GenoVirA", features = c("nCount_RNA", "nFeature_RNA"), cols = mycol2, pt.size = 0, ncol = 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_C.pdf", plot = vpl3, width = 4, height = 6, units = "in", dpi = 300)

vpl4 <- VlnPlot(object = seuObjSel, group.by = "GenoVirA", features = c("INTRONIC_READS_RATIO"), cols = mycol2, pt.size = 0) + ylim(0, 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_D.pdf", plot = vpl4, width = 4, height = 3, units = "in", dpi = 300)



## COLOR PALETTE FOR CELLTYPE CLUSTERS
mycol3 <- c(d100, d101, d102, d103, d104, d105, d111, d117, d206, d207, d210, d212, e08, e13)

## UMAP PLOTS FOR CELLTYPE CLUSTERS
pumap1a <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationB", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol3) #+ NoLegend()
ggsave("JUL2023_UMAP_SPN_F.PDF", plot = pumap1a, width = 5, height = 4, units = "in", dpi = 300)

pumap1b <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationB", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = mycol3) + NoLegend()
ggsave("JUL2023_UMAP_SPN_G.PDF", plot = pumap1b, width = 4, height = 4, units = "in", dpi = 300)

pumap1c <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 0.5, group.by = "AnnotationB", raster = TRUE, raster.dpi = c(1024, 1024), label = TRUE, label.size = 2, repel = TRUE, colors_use = mycol3) + NoLegend()
ggsave("JUL2023_UMAP_SPN_H.PDF", plot = pumap1c, width = 4, height = 4, units = "in", dpi = 300)

## VIOLIN PLOTS FOR CELLTYPE CLUSTERS
vpl1 <- VlnPlot(object = seuObjSel, group.by = "AnnotationB", features = c("nCount_RNA", "nFeature_RNA"), cols = mycol3, pt.size = 0, ncol = 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_E.pdf", plot = vpl1, width = 6, height = 6, units = "in", dpi = 300)

vpl2 <- VlnPlot(object = seuObjSel, group.by = "AnnotationB", features = c("INTRONIC_READS_RATIO"), cols = mycol3, pt.size = 0) + ylim(0, 1) + NoLegend()
ggsave(filename = "JUL2023_VIOLINPLOT_SPN_F.pdf", plot = vpl2, width = 6, height = 3, units = "in", dpi = 300)


## UMAP PLOTS FOR WPRE AND MCHERRY
pumap1d <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "WPRE_GTEQ_6", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = c("grey60", "blue")) #+ NoLegend()
ggsave("JUL2023_UMAP_SPN_I.PDF", plot = pumap1d, width = 4.5, height = 4, units = "in", dpi = 300)

pumap1e <- DimPlot_scCustom(seurat_object = seuObjSel, pt.size = 1, group.by = "MCHERRY_GTEQ_6", raster = TRUE, raster.dpi = c(1024, 1024), label = FALSE, colors_use = c("grey60", "red")) #+ NoLegend()
ggsave("JUL2023_UMAP_SPN_J.PDF", plot = pumap1e, width = 4.5, height = 4, units = "in", dpi = 300)


## SAVE RDATA
save(seuObjSel, mycol1, mycol2, mycol3, file = "JUL2023_SEURAT_SPN.RData")



## DGE VOLCANO PLOTS
rm(list = ls())

## GENES TO SHOW
genes2show <- c("Grin1", "Grin2a", "Kcnk2", "Gprin3", "Kcnj3", "Syngap1", "Adora2a", "Cnr1", "Syt6", "Dlg4", "Dlgap3", "Grm1", "Grm5", "Grm8", "Ppp1r1b")

## GENES_TO_REMOVE
chrmgenes <- scan("/GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n")
chrxgenes <- scan("GENCODE_vM17_MM10_ChrX_GENES.txt", what = "", sep = "\n")
chrygenes <- scan("GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n")
rpgenes <- scan("GENCODE_vM17_MM10_RP_GENES.txt", what = "", sep = "\n")
gmgenes <- scan("GENCODE_vM17_MM10_Gm_GENES.txt", what = "", sep = "\n")
genes2remove <- c(chrmgenes, chrxgenes, chrygenes, rpgenes, gmgenes)


## DGE TABLES
## All iSPNs, No Covariates
load("NK_DGE_PSEUDOBULK_EDGER_WT_CTL_x_CKO_CTL_DGE_AVG.RData")
# degavg, selLrt
degavgsel <- degavg[!row.names(degavg) %in% genes2remove,]
degavgsel$logFC <- -1 * degavgsel$logFC

fdrcutoff <- 0.05
fccutoff <- 1.3

## filter dge table
dge_sig_temp <- degavgsel[degavgsel$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])


## volcano plot
limit_x_max <- round(max(abs(degavgsel$logFC)) + 1.5)
limit_x_min <- -1 * limit_x_max
limit_y_max <- round(max(-log10(degavgsel$FDR)) + 10)
anno_x_max <- limit_x_max - 0.5
anno_x_min <- limit_x_min + 0.5
anno_y_max <- limit_y_max - 5

cluSelected <- "iSPN"
mygroup1 <- "CKO_CTL"
mygroup2 <- "WT_CTL"

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals.colour <- ifelse(
    degavgsel$logFC < -log2(fccutoff), 'green4',
        ifelse(degavgsel$logFC > log2(fccutoff), 'red3',
        'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'up'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'no change'
names(keyvals.colour)[keyvals.colour == 'green4'] <- 'down'

keyvalcol <- ifelse(degavgsel$logFC <= -log2(fccutoff) & degavgsel$FDR <= fdrcutoff, 'green4',
             ifelse(degavgsel$logFC >= log2(fccutoff) & degavgsel$FDR <= fdrcutoff, 'red3',
             ifelse(degavgsel$FDR <= fdrcutoff, 'black',
             'grey')))
keyvalcol[is.na(keyvalcol)] <- 'grey'
names(keyvalcol)[keyvalcol == 'red3'] <- 'up'
names(keyvalcol)[keyvalcol == 'black'] <- 'no_change'
names(keyvalcol)[keyvalcol == 'green4'] <- 'down'
names(keyvalcol)[keyvalcol == 'grey'] <- 'no_sig'


VOLplot <- EnhancedVolcano(degavgsel, 
                            lab = rownames(degavgsel), 
                            x = 'logFC', 
                            y = 'FDR', 
                            pointSize = 1,
                            xlim = c(limit_x_min, limit_x_max),
                            xlab = bquote(~log[2] ~ "fold change"), 
                            ylab = bquote(~-log[10] ~ italic("adj. p-val")), 
                            pCutoff = fdrcutoff, 
                            FCcutoff = log2(fccutoff), 
                            title = cluSelected, 
                            colCustom = keyvalcol,
                            colAlpha = 1,
                            subtitle = paste("Up-regulated Genes: ", length(dge_sig_up), " | ", "Down-regulated Genes: ", length(dge_sig_dn), sep = ""), caption = "",
                            legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(filename = paste("JUL2023_VOLCANO_ALL_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO_1.PDF", sep = ""), plot = VOLplot, width = 8, height = 9, units = "in", dpi = 300, useDingbats = FALSE)

options(ggrepel.max.overlaps = Inf)
VOLplot <- EnhancedVolcano(degavgsel, 
                            lab = rownames(degavgsel), 
                            x = 'logFC', 
                            y = 'FDR', 
                            pointSize = 1,
                            selectLab = genes2show,
                            drawConnectors = TRUE,
                            widthConnectors = 1,
                            maxoverlapsConnectors = Inf,
                            xlim = c(limit_x_min, limit_x_max),
                            xlab = bquote(~log[2] ~ "fold change"), 
                            ylab = bquote(~-log[10] ~ italic("adj. p-val")), 
                            pCutoff = fdrcutoff, 
                            FCcutoff = log2(fccutoff), 
                            title = cluSelected, 
                            colCustom = keyvalcol,
                            colAlpha = 1,
                            subtitle = paste("Up-regulated Genes: ", length(dge_sig_up), " | ", "Down-regulated Genes: ", length(dge_sig_dn), sep = ""), caption = "",
                            legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(filename = paste("JUL2023_VOLCANO_ALL_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO_2.PDF", sep = ""), plot = VOLplot, width = 8, height = 9, units = "in", dpi = 300, useDingbats = FALSE)


VOLplot <- EnhancedVolcano(degavgsel, 
                            lab = rownames(degavgsel), 
                            x = 'logFC', 
                            y = 'FDR', 
                            pointSize = 1,
                            labSize = 0,
                            xlim = c(limit_x_min, limit_x_max),
                            xlab = bquote(~log[2] ~ "fold change"), 
                            ylab = bquote(~-log[10] ~ italic("adj. p-val")), 
                            pCutoff = fdrcutoff, 
                            FCcutoff = log2(fccutoff), 
                            title = cluSelected, 
                            colCustom = keyvalcol,
                            colAlpha = 1,
                            subtitle = paste("Up-regulated Genes: ", length(dge_sig_up), " | ", "Down-regulated Genes: ", length(dge_sig_dn), sep = ""), caption = "",
                            legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(filename = paste("JUL2023_VOLCANO_ALL_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO_3.PDF", sep = ""), plot = VOLplot, width = 8, height = 9, units = "in", dpi = 300, useDingbats = FALSE)




## DGE DISEASE ENRICHMENT
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
library(tidyr)
library(RColorBrewer)
library(UpSetR)
library(ggpubr)
library(nichenetr) ## has function to convert human genes to mouse
library(SuperExactTest)
library(ggplot2)
library(ggnewscale)
library(rlist)
library(wesanderson)
library(edgeR)


## PREPARE REFERENCE GENE LISTS
## Load reference geneset
load("GeneSets_Disorders.RData")

asd13_human <- as.vector(GeneSets$`ASD_1-3`$Gene)
asd_human <- as.vector(GeneSets$ASD$Gene)
fmrp_human <- as.vector(GeneSets$FMRP$Gene)
id_human <- as.vector(GeneSets$ID$Gene)
syn_human <- as.vector(GeneSets$SYN$Gene)

asd13_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asd13_human))  # 467 genes
asd_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asd_human))  # 824 genes
fmrp_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = fmrp_human))  # 773 genes
id_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = id_human))  # 1010 genes
syn_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = syn_human))  # 1809 genes

asd13_mouse_df <- as.data.frame(asd13_mouse)
write.table(asd13_mouse_df, "FOR_NK_ASD_1-3_MOUSE.txt", row.names = F, col.names = F, quote = F, sep = "\t")
asd_mouse_df <- as.data.frame(asd_mouse)
write.table(asd_mouse_df, "FOR_NK_ASD_MOUSE.txt", row.names = F, col.names = F, quote = F, sep = "\t")
id_mouse_df <- as.data.frame(id_mouse)
write.table(id_mouse_df, "FOR_NK_ID_MOUSE.txt", row.names = F, col.names = F, quote = F, sep = "\t")


nk_syn <- read.table("NK_SynapticGenes.txt", sep = "\t", header = TRUE, fill = TRUE)

nk_syn_ex <- nk_syn$Excitatory_Synapse[grepl("^[A-Z]", nk_syn$Excitatory_Synapse)]
nk_syn_in <- nk_syn$Inhibitory_Synapse[grepl("^[A-Z]", nk_syn$Inhibitory_Synapse)]
nk_syn_pre <- nk_syn$Presynapse[grepl("^[A-Z]", nk_syn$Presynapse)]
nk_syn_post <- nk_syn$Postsynapse[grepl("^[A-Z]", nk_syn$Postsynapse)]
nk_syn_glut <- nk_syn$Glutamatergic_Synapse[grepl("^[A-Z]", nk_syn$Glutamatergic_Synapse)]
nk_syn_gaba <- nk_syn$Gabaergic_Synapse[grepl("^[A-Z]", nk_syn$Gabaergic_Synapse)]
nk_syn_dopa <- nk_syn$Dopaminergic_Synapse[grepl("^[A-Z]", nk_syn$Dopaminergic_Synapse)]
nk_syn_chol <- nk_syn$Cholinergic_Synapse[grepl("^[A-Z]", nk_syn$Cholinergic_Synapse)]
nk_syn_sero <- nk_syn$Serotonergic_synapse[grepl("^[A-Z]", nk_syn$Serotonergic_synapse)]
nk_syn_all <- nk_syn$All_Genes[grepl("^[A-Z]", nk_syn$All_Genes)]


nk_voltage <- scan("NK_VOLTAGE_GATED.txt", what = "", sep = "\n")
nk_ligand <- scan("NK_LIGAND_GATED.txt", what = "", sep = "\n")
# nk_othch <- scan("NK_OTHERS.txt", what = "", sep = "\n")
# nk_trans <- scan("NK_TRANSPORTERS.txt", what = "", sep = "\n")
# nk_catalityc <- scan("NK_CATALYTIC.txt", what = "", sep = "\n")
# nk_enz <- scan("NK_ENZYMES.txt", what = "", sep = "\n")
# nk_othprot <- scan("NK_OTHER_PROTEINS.txt", what = "", sep = "\n")
# nk_syn <- scan("NK_SYNAPSE.txt", what = "", sep = "\n")
# nk_channel <- scan("NK_CHANNELS.txt", what = "", sep = "\n")

## create a list for all reference gene lists
myref1 <- list("01_ASD_1-3" = asd13_mouse, 
               "02_ASD" = asd_mouse, 
               "03_FMRP" = fmrp_mouse, 
               "04_ID" = id_mouse)

# myref2 <- list("01_SYNAPSE_EXC" = nk_syn_ex,
#                "02_SYNAPSE_INH" = nk_syn_in)

# myref3 <- list("01_SYNAPSE_PRE" = nk_syn_pre,
#                "02_SYNAPSE_POST" = nk_syn_post)

myref4 <- list("01_SYNAPSE_GLUT" = nk_syn_glut,
               "02_SYNAPSE_GABA" = nk_syn_gaba,
               "03_SYNAPSE_DOPA" = nk_syn_dopa,
               "04_SYNAPSE_CHOL" = nk_syn_chol,
               "05_SYNAPSE_SERO" = nk_syn_sero)

myref5 <- list("01_VOLTAGE" = nk_voltage,
               "02_LIGAND" = nk_ligand)

# myref6 <- list("01_SYNAPSE_ALL" = nk_syn_all)



## READ & FILTER DEG TABLES
## GENES_TO_REMOVE
chrmgenes <- scan("GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n")
chrxgenes <- scan("GENCODE_vM17_MM10_ChrX_GENES.txt", what = "", sep = "\n")
chrygenes <- scan("GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n")
rpgenes <- scan("GENCODE_vM17_MM10_RP_GENES.txt", what = "", sep = "\n")
gmgenes <- scan("GENCODE_vM17_MM10_Gm_GENES.txt", what = "", sep = "\n")
genes2remove <- c(chrmgenes, chrxgenes, chrygenes, rpgenes, gmgenes)

## Selected iSPNs, wpre >= 6, mcherry >= 6, No Covariates
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/DGE_PSEUDOBULK_JUN2023_2/iSPNs_ALL_NOCOV/NK_DGE_PSEUDOBULK_EDGER_WT_CTL_x_CKO_CTL_DGE_AVG.RData")
# degavg, selLrt
degavgsel <- degavg[!row.names(degavg) %in% genes2remove,]
degavgsel$logFC <- -1 * degavgsel$logFC

fdrcutoff <- 0.05
fccutoff <- 1.3

## filter dge table
dge_sig_temp <- degavgsel[degavgsel$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])

nkdeg <- list("01_CKO_CTL_x_WT_CTL_UP" = dge_sig_up,
              "02_CKO_CTL_x_WT_CTL_DOWN" = dge_sig_dn)


## Perform Super Exact Test
setout <- list()
setoutAll <- list()
i <- 0

for(mod in names(nkdeg))
  {
  i <- i + 1
  newGenes <- list.append(myref5, mod = nkdeg[[mod]])
  names(newGenes)[[length(newGenes)]] <- mod
  setres <- supertest(newGenes, n = 7918, degree = 2)  # 7918 is background number of genes for CKO_CTL vs WT_CTL
  setresData <- as.data.frame(summary(setres)$Table)
  setoutAll[[i]] <- setresData
  setresDataSel <- setresData[grep(mod, setresData$Intersections),c(1, 3, 5, 6)]
  setout[[i]] <- setresDataSel
  print(paste(mod, i, sep = " "))
  }

names(setoutAll) <- names(nkdeg)
names(setout) <- names(nkdeg)

setoutData <- Reduce("rbind", setout)

save(setout, setoutData, setoutAll, file = "NK_DGE_ENRICHMENT_SET.RData")

## REORGANIZE DATA FOR PLOTS
setoutDataTemp <- as.data.frame(matrix(unlist(strsplit(setoutData$Intersections, " & ")), byrow = T, ncol = 2))
colnames(setoutDataTemp) <- c("REF_GENES", "CELLTYPE_DEG")
setoutDataTemp$Intersections <- setoutData$Intersections

setoutData2 <- merge(setoutData, setoutDataTemp, by = "Intersections")
setoutData2$NegLog10 <- -log10(setoutData2$P.value + 10^-10)

write.table(setoutData2, "NK_DGE_ENRICHMENT_SET_FOR_PLOT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

setoutData3 <- setoutData2

pdf("NK_DGE_ENRICHMENT_SET_FOR_PLOT_BUBBLE_PLOT.pdf", width = 8, height = 6)
ggscatter(setoutData3, 
          x = "CELLTYPE_DEG",
          y = "REF_GENES",
          size = "FE",
          color = "NegLog10",
          alpha = 1,
          xlab = "",
          ylab = "",) +
  theme_minimal() +
  rotate_x_text(angle = 45)+
  coord_flip()+
  scale_size(range = c(2, 12))+ 
  # gradient_color(c("grey95", "green4", "gold", "red3"))
  gradient_color(c("grey95", "blue2"))
dev.off()





## DGE HEATMAPS WITH DENDROGRAMS
rm(list = ls())
library(pheatmap)
library(RColorBrewer)

## GENES_TO_REMOVE
chrmgenes <- scan("./../GENES_TO_REMOVE/GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n")
chrxgenes <- scan("./../GENES_TO_REMOVE/GENCODE_vM17_MM10_ChrX_GENES.txt", what = "", sep = "\n")
chrygenes <- scan("./../GENES_TO_REMOVE/GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n")
rpgenes <- scan("./../GENES_TO_REMOVE/GENCODE_vM17_MM10_RP_GENES.txt", what = "", sep = "\n")
gmgenes <- scan("./../GENES_TO_REMOVE/GENCODE_vM17_MM10_Gm_GENES.txt", what = "", sep = "\n")
genes2remove <- c(chrmgenes, chrxgenes, chrygenes, rpgenes, gmgenes)

## Selected iSPNs, wpre >= 6, mcherry >= 6, No Covariates
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/DGE_PSEUDOBULK_JUN2023_2/iSPNs_WPRE_GTEQ_6_MCHERRY_GTEQ_6_NOCOV/NK_DGE_PSEUDOBULK_EDGER_WT_CTL_x_CKO_CTL_DGE_AVG.RData")
degavgsel1 <- degavg[!row.names(degavg) %in% genes2remove,]
rm(degavg)
rm(selLrt)

load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/DGE_PSEUDOBULK_JUN2023_2/iSPNs_WPRE_GTEQ_6_MCHERRY_GTEQ_6_NOCOV/NK_DGE_PSEUDOBULK_EDGER_CKO_CTL_x_CKO_RES_DGE_AVG.RData")
degavgsel2 <- degavg[!row.names(degavg) %in% genes2remove,]
rm(degavg)
rm(selLrt)

load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/DGE_PSEUDOBULK_JUN2023_2/iSPNs_WPRE_GTEQ_6_MCHERRY_GTEQ_6_NOCOV/NK_DGE_PSEUDOBULK_EDGER_WT_CTL_x_CKO_RES_DGE_AVG.RData")
degavgsel3 <- degavg[!row.names(degavg) %in% genes2remove,]
rm(degavg)
rm(selLrt)


## avg exp table
avg_exp <- merge(degavgsel1[,c("WT_CTL_P18_C1", "WT_CTL_P18_C2", "WT_CTL_P18_C3", "CKO_CTL_P18_3BF1", "CKO_CTL_P18_3BM3", "CKO_CTL_P18_3BM4")], degavgsel3[,c("CKO_RES_P18_3BM8", "CKO_RES_P18_3CF5", "CKO_RES_P18_3CF8")], by = "row.names")
row.names(avg_exp) <- avg_exp$Row.names
avg_exp$Row.names <- NULL
colnames(avg_exp) <- c("WT_CTL_1", "WT_CTL_2", "WT_CTL_3", "CKO_CTL_1", "CKO_CTL_2", "CKO_CTL_3", "CKO_RES_1", "CKO_RES_2", "CKO_RES_3")


## Filter DGEs
fdrcutoff <- 0.05
fccutoff <- 1.3

## filter dge table
degavgsel1$logFC <- -1 * degavgsel1$logFC
dge_sig_temp <- degavgsel1[degavgsel1$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])

## heatmap
avg_exp_up <- avg_exp[row.names(avg_exp) %in% dge_sig_up,]
avg_exp_dn <- avg_exp[row.names(avg_exp) %in% dge_sig_dn,]
avg_exp_all <- rbind(avg_exp_up, avg_exp_dn)

# prefix <- "CKO_CTL_x_WT_CTL"
# pheatmap(avg_exp_all, scale = "row", cluster_rows = T, cluster_cols = T, border_color = NA, show_rownames = F, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_1.pdf", sep = ""), width = 6, height = 6)


# ## heatmap2
nk_sel_up <- scan("NK_SEL_UP.txt", sep = "\n", what = "") # 59
nk_sel_dn <- scan("NK_SEL_DN.txt", sep = "\n", what = "") # 19

avg_exp_sel_up <- avg_exp[row.names(avg_exp) %in% nk_sel_up,]
avg_exp_sel_dn <- avg_exp[row.names(avg_exp) %in% nk_sel_dn,]
avg_exp_sel <- rbind(avg_exp_sel_up, avg_exp_sel_dn)


prefix <- "CKO_CTL_x_WT_CTL"
# pheatmap(avg_exp_sel, scale = "row", cluster_rows = T, cluster_cols = T, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_2.pdf", sep = ""), width = 6, height = 12)
# pheatmap(avg_exp_sel, scale = "row", cluster_rows = T, cluster_cols = F, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_3.pdf", sep = ""), width = 6, height = 12)

# avg_exp_sel_sorted <- avg_exp_sel[order(avg_exp_sel$WT_CTL_1, decreasing = TRUE),]
# pheatmap(avg_exp_sel, scale = "row", cluster_rows = F, cluster_cols = F, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_4.pdf", sep = ""), width = 6, height = 12)

# pheatmap(avg_exp_sel, scale = "row", cluster_rows = F, cluster_cols = T, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_5.pdf", sep = ""), width = 6, height = 12)


# pheatmap(avg_exp_sel_up, scale = "row", cluster_rows = T, cluster_cols = T, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_2_DOWN.pdf", sep = ""), width = 6, height = 15)
# pheatmap(avg_exp_sel_dn, scale = "row", cluster_rows = T, cluster_cols = T, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("JUL2023_HEATMAP_", prefix, "_2_UP.pdf", sep = ""), width = 6, height = 6)


nk_sel_order_up <- scan("NK_SEL_GENES_HEATMAP_UP2.txt", sep = "\n", what = "") # 19
nk_sel_order_dn <- scan("NK_SEL_GENES_HEATMAP_DN.txt", sep = "\n", what = "") # 59

avg_exp_sel_order_up <- avg_exp_sel_up[nk_sel_order_dn,]
avg_exp_sel_order_dn <- avg_exp_sel_dn[nk_sel_order_up,]
avg_exp_sel_order <- rbind(avg_exp_sel_order_up, avg_exp_sel_order_dn)

pheatmap(avg_exp_sel_order, scale = "row", cluster_rows = F, cluster_cols = T, border_color = NA, show_rownames = T, main = prefix, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(5000), filename = paste("AUG2023_HEATMAP_", prefix, "_6.pdf", sep = ""), width = 6, height = 12)

```
