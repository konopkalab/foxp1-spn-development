# DIFFERENTIAL GENE EXPRESSION ANALYSIS

### Load Modules
```{bash}
module purge && module load shared slurm python/3.7.x-anaconda
module load gcc/8.3.0
module load hdf5_18/1.8.17
module load R/4.1.1-gccmkl
```

### DGE Analysis
```{R}
## Load libraries
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(plyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(scran)
library(scater)
library(SingleCellExperiment)
library(edgeR)


##-------------------------------------------------------
## WT_CTL x CKO_CTL

rm(list = ls())
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED.RData")
# seuObj
seuMeta <- as.data.frame(seuObj@meta.data)

##-------------------------------------------------------
# ## filter for WPRE+/MCHERRY+ nuclei
# wpre_gteq_6_cells <- row.names(seuMeta[seuMeta$WPRE_GTEQ_6 == 1,])
# mcherry_gteq_6_cells <- row.names(seuMeta[seuMeta$MCHERRY_GTEQ_6 == 1,])
# cells2keep <- unique(sort(c(wpre_gteq_6_cells, mcherry_gteq_6_cells)))

# seuObjFilt <- subset(seuObj, cells = cells2keep, slot = "counts")
# seuObjFilt <- NormalizeData(seuObjFilt)
seuObjFilt <- seuObj

##-------------------------------------------------------
## variables
cluSelected <- "iSPN"
fdrcutoff <- 0.05
fccutoff <- 1.2

##-------------------------------------------------------
## subset for a cell-type iSPN
Idents(seuObjFilt) <- "Annotation"
print(table(seuObjFilt$Annotation))
#     dSPN         eSPN         iSPN Undetermined 
#    24230         3648         9484         4503

seuObjSelTemp <- subset(x = seuObjFilt, idents = c(cluSelected))
seuObjSel <- seuObjSelTemp

## create SCE object for all genes
strSCE <- as.SingleCellExperiment(seuObjSel)

strGroups <- colData(strSCE)[, c("GenoVirAgeSample", "GenoVir")]
strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
strAggMat <- strPseudoSCE@assays@data$sum
colnames(strAggMat) <- colData(strPseudoSCE)[['GenoVirAgeSample']]

## filter genes
wtctlExp = which(rowSums(strAggMat[, grepl("WT_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoctlExp = which(rowSums(strAggMat[, grepl("CKO_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoresExp = which(rowSums(strAggMat[, grepl("CKO_RES", colnames(strAggMat))] <= 10) == 0) %>% names
expgns = Reduce(union, list(wtctlExp, ckoctlExp, ckoresExp))

print(length(expgns))
# 9984

mygroup1 <- "WT_CTL"
mygroup2 <- "CKO_CTL"

seuObjSel <- subset(seuObjSel, subset = GenoVir %in% c("WT_CTL", "CKO_CTL"))
seuObjSel$GenoVir = droplevels(as.factor(seuObjSel$GenoVir))

## create SCE object for filtered genes
DefaultAssay(seuObjSel) <- "RNA"
selSCE <- as.SingleCellExperiment(seuObjSel)

## Pseudobulk SCE assay
selGroups <- colData(selSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
selGroups <- colData(selPseudoSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selGroups$GenoVirAgeSample <- factor(selGroups$GenoVirAgeSample)
selGroups$GenoVir <- factor(selGroups$GenoVir, levels = c("CKO_CTL", "WT_CTL"))
selGroups$Batch <- factor(selGroups$Batch)
selGroups$Sex <- factor(selGroups$Sex)

selAggMat <- selPseudoSCE@assays@data$sum
colnames(selAggMat) = colData(selPseudoSCE)[['GenoVirAgeSample']]
selAggMat <- selAggMat[expgns,]

## perform differential gene expression
selDGEL <- DGEList(counts = selAggMat)
selDGEL <- calcNormFactors(selDGEL)
# selDesign <- model.matrix(~GenoVir + Batch + Sex, data = selGroups)
# selDesign <- model.matrix(~GenoVir + Sex, data = selGroups)
# selDesign <- model.matrix(~Sex + GenoVir, data = selGroups)
selDesign <- model.matrix(~GenoVir, data = selGroups)
selDGEL <- estimateDisp(selDGEL, selDesign)

selFit <- glmFit(selDGEL, selDesign)
selLrt <- glmLRT(selFit, coef = 'GenoVirWT_CTL')
selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame

## average gene expression per sample
avg_exp <- as.data.frame(log2(cpm(selAggMat) + 1))

## merge dge table and avg gene expression
degavg <- merge(avg_exp, selRes, by = "row.names", all.y = TRUE)
row.names(degavg) <- degavg$Row.names
degavg$Row.names <- NULL

## write and save dge table
write.table(degavg, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_TABLE.tsv", sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
save(degavg, selLrt, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG.RData", sep = ""))

## filter dge table
dge_sig_temp <- degavg[degavg$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])

## volcano plot
limit_x_max <- round(max(abs(degavg$logFC)) + 1)
limit_x_min <- -1 * limit_x_max
limit_y_max <- round(max(-log10(degavg$FDR)) + 10)
anno_x_max <- limit_x_max - 0.5
anno_x_min <- limit_x_min + 0.5
anno_y_max <- limit_y_max - 5

pvol1 <- ggplot(data = degavg, aes(x = logFC, y = -log10(FDR))) +
            geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
            geom_point(data=dge_sig[dge_sig$logFC > 0,], colour="green3", size=1.2, shape=16, alpha=0.8) +
            geom_point(data=dge_sig[dge_sig$logFC < 0,], colour="royalblue", size=1.2, shape=16, alpha=0.8) +
            annotate(geom = "text", x = anno_x_max, y = anno_y_max, label = length(dge_sig_up)) +
            annotate(geom = "text", x = anno_x_min, y = anno_y_max, label = length(dge_sig_dn)) +
            theme_bw() +
            geom_vline(xintercept = log2(fccutoff), colour = "red", linetype = 2) + geom_vline(xintercept = -1 * log2(fccutoff), colour = "red", linetype = 2) +
            geom_hline(yintercept = -log10(fdrcutoff), colour = "red", linetype = 2) +
            xlim(limit_x_min, limit_x_max) +
            ylim(0, limit_y_max) +
            labs(x="Avg log Fold Change", y="- log10 FDR")
ggsave(filename = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO.PDF", sep = ""), plot = pvol1, width = 4, height = 4, units = "in", dpi = 300, useDingbats = FALSE)






##-------------------------------------------------------
## CKO_CTL x CKO_RES
rm(list = ls())
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED.RData")
# seuObj
seuMeta <- as.data.frame(seuObj@meta.data)

##-------------------------------------------------------
# ## filter for WPRE+/MCHERRY+ nuclei
# wpre_gteq_6_cells <- row.names(seuMeta[seuMeta$WPRE_GTEQ_6 == 1,])
# mcherry_gteq_6_cells <- row.names(seuMeta[seuMeta$MCHERRY_GTEQ_6 == 1,])
# cells2keep <- unique(sort(c(wpre_gteq_6_cells, mcherry_gteq_6_cells)))

# seuObjFilt <- subset(seuObj, cells = cells2keep, slot = "counts")
# seuObjFilt <- NormalizeData(seuObjFilt)
seuObjFilt <- seuObj

##-------------------------------------------------------
## variables
cluSelected <- "iSPN"
fdrcutoff <- 0.05
fccutoff <- 1.2

##-------------------------------------------------------
## subset for a cell-type iSPN
Idents(seuObjFilt) <- "Annotation"
print(table(seuObjFilt$Annotation))
#     dSPN         eSPN         iSPN Undetermined 
#    24230         3648         9484         4503

seuObjSelTemp <- subset(x = seuObjFilt, idents = c(cluSelected))
seuObjSel <- seuObjSelTemp

## create SCE object for all genes
strSCE <- as.SingleCellExperiment(seuObjSel)

strGroups <- colData(strSCE)[, c("GenoVirAgeSample", "GenoVir")]
strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
strAggMat <- strPseudoSCE@assays@data$sum
colnames(strAggMat) <- colData(strPseudoSCE)[['GenoVirAgeSample']]

## filter genes
wtctlExp = which(rowSums(strAggMat[, grepl("WT_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoctlExp = which(rowSums(strAggMat[, grepl("CKO_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoresExp = which(rowSums(strAggMat[, grepl("CKO_RES", colnames(strAggMat))] <= 10) == 0) %>% names
expgns = Reduce(union, list(wtctlExp, ckoctlExp, ckoresExp))

print(length(expgns))
# 9984

mygroup1 <- "CKO_CTL"
mygroup2 <- "CKO_RES"

seuObjSel <- subset(seuObjSel, subset = GenoVir %in% c("CKO_CTL", "CKO_RES"))
seuObjSel$GenoVir = droplevels(as.factor(seuObjSel$GenoVir))

## create SCE object for filtered genes
DefaultAssay(seuObjSel) <- "RNA"
selSCE <- as.SingleCellExperiment(seuObjSel)

## Pseudobulk SCE assay
selGroups <- colData(selSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
selGroups <- colData(selPseudoSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selGroups$GenoVirAgeSample <- factor(selGroups$GenoVirAgeSample)
selGroups$GenoVir <- factor(selGroups$GenoVir, levels = c("CKO_RES", "CKO_CTL"))
selGroups$Batch <- factor(selGroups$Batch)
selGroups$Sex <- factor(selGroups$Sex)

selAggMat <- selPseudoSCE@assays@data$sum
colnames(selAggMat) = colData(selPseudoSCE)[['GenoVirAgeSample']]
selAggMat <- selAggMat[expgns,]

## perform differential gene expression
selDGEL <- DGEList(counts = selAggMat)
selDGEL <- calcNormFactors(selDGEL)
# selDesign <- model.matrix(~GenoVir + Batch + Sex, data = selGroups)
# selDesign <- model.matrix(~GenoVir + Sex, data = selGroups)
# selDesign <- model.matrix(~Sex + GenoVir, data = selGroups)
selDesign <- model.matrix(~GenoVir, data = selGroups)
selDGEL <- estimateDisp(selDGEL, selDesign)

selFit <- glmFit(selDGEL, selDesign)
selLrt <- glmLRT(selFit, coef = 'GenoVirCKO_CTL')
selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame

## average gene expression per sample
avg_exp <- as.data.frame(log2(cpm(selAggMat) + 1))

## merge dge table and avg gene expression
degavg <- merge(avg_exp, selRes, by = "row.names", all.y = TRUE)
row.names(degavg) <- degavg$Row.names
degavg$Row.names <- NULL

## write and save dge table
write.table(degavg, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_TABLE.tsv", sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
save(degavg, selLrt, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG.RData", sep = ""))

## filter dge table
dge_sig_temp <- degavg[degavg$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])

## volcano plot
limit_x_max <- round(max(abs(degavg$logFC)) + 1)
limit_x_min <- -1 * limit_x_max
limit_y_max <- round(max(-log10(degavg$FDR)) + 10)
anno_x_max <- limit_x_max - 0.5
anno_x_min <- limit_x_min + 0.5
anno_y_max <- limit_y_max - 5

pvol1 <- ggplot(data = degavg, aes(x = logFC, y = -log10(FDR))) +
            geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
            geom_point(data=dge_sig[dge_sig$logFC > 0,], colour="green3", size=1.2, shape=16, alpha=0.8) +
            geom_point(data=dge_sig[dge_sig$logFC < 0,], colour="royalblue", size=1.2, shape=16, alpha=0.8) +
            annotate(geom = "text", x = anno_x_max, y = anno_y_max, label = length(dge_sig_up)) +
            annotate(geom = "text", x = anno_x_min, y = anno_y_max, label = length(dge_sig_dn)) +
            theme_bw() +
            geom_vline(xintercept = log2(fccutoff), colour = "red", linetype = 2) + geom_vline(xintercept = -1 * log2(fccutoff), colour = "red", linetype = 2) +
            geom_hline(yintercept = -log10(fdrcutoff), colour = "red", linetype = 2) +
            xlim(limit_x_min, limit_x_max) +
            ylim(0, limit_y_max) +
            labs(x="Avg log Fold Change", y="- log10 FDR")
ggsave(filename = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO.PDF", sep = ""), plot = pvol1, width = 4, height = 4, units = "in", dpi = 300, useDingbats = FALSE)







##-------------------------------------------------------
## WT_CTL x CKO_RES
rm(list = ls())
load("/work/Neuroinformatics_Core/akulk1/CELLRANGER/NITIN_DATA/E_INTEGRATION_HARMONY/SUBCLUSTERING_SPN_01/SEURAT_NK_P18_CKO_RES_HARMONY_WPRE_MCHERRY_ANNOTATED.RData")
# seuObj
seuMeta <- as.data.frame(seuObj@meta.data)

##-------------------------------------------------------
# ## filter for WPRE+/MCHERRY+ nuclei
# wpre_gteq_6_cells <- row.names(seuMeta[seuMeta$WPRE_GTEQ_6 == 1,])
# mcherry_gteq_6_cells <- row.names(seuMeta[seuMeta$MCHERRY_GTEQ_6 == 1,])
# cells2keep <- unique(sort(c(wpre_gteq_6_cells, mcherry_gteq_6_cells)))

# seuObjFilt <- subset(seuObj, cells = cells2keep, slot = "counts")
# seuObjFilt <- NormalizeData(seuObjFilt)
seuObjFilt <- seuObj

##-------------------------------------------------------
## variables
cluSelected <- "iSPN"
fdrcutoff <- 0.05
fccutoff <- 1.2

##-------------------------------------------------------
## subset for a cell-type iSPN
Idents(seuObjFilt) <- "Annotation"
print(table(seuObjFilt$Annotation))
#     dSPN         eSPN         iSPN Undetermined 
#    10589         1099         4018         1796

seuObjSelTemp <- subset(x = seuObjFilt, idents = c(cluSelected))
seuObjSel <- seuObjSelTemp

## create SCE object for all genes
strSCE <- as.SingleCellExperiment(seuObjSel)

strGroups <- colData(strSCE)[, c("GenoVirAgeSample", "GenoVir")]
strPseudoSCE <- sumCountsAcrossCells(strSCE, strGroups)
strAggMat <- strPseudoSCE@assays@data$sum
colnames(strAggMat) <- colData(strPseudoSCE)[['GenoVirAgeSample']]

## filter genes
wtctlExp = which(rowSums(strAggMat[, grepl("WT_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoctlExp = which(rowSums(strAggMat[, grepl("CKO_CTL", colnames(strAggMat))] <= 10) == 0) %>% names
ckoresExp = which(rowSums(strAggMat[, grepl("CKO_RES", colnames(strAggMat))] <= 10) == 0) %>% names
expgns = Reduce(union, list(wtctlExp, ckoctlExp, ckoresExp))

print(length(expgns))
# 9984

mygroup1 <- "WT_CTL"
mygroup2 <- "CKO_RES"

seuObjSel <- subset(seuObjSel, subset = GenoVir %in% c("WT_CTL", "CKO_RES"))
seuObjSel$GenoVir = droplevels(as.factor(seuObjSel$GenoVir))

## create SCE object for filtered genes
DefaultAssay(seuObjSel) <- "RNA"
selSCE <- as.SingleCellExperiment(seuObjSel)

## Pseudobulk SCE assay
selGroups <- colData(selSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
selGroups <- colData(selPseudoSCE)[, c("GenoVirAgeSample", "GenoVir", "Batch", "Sex")]
selGroups$GenoVirAgeSample <- factor(selGroups$GenoVirAgeSample)
selGroups$GenoVir <- factor(selGroups$GenoVir, levels = c("CKO_RES", "WT_CTL"))
selGroups$Batch <- factor(selGroups$Batch)
selGroups$Sex <- factor(selGroups$Sex)

selAggMat <- selPseudoSCE@assays@data$sum
colnames(selAggMat) = colData(selPseudoSCE)[['GenoVirAgeSample']]
selAggMat <- selAggMat[expgns,]

## perform differential gene expression
selDGEL <- DGEList(counts = selAggMat)
selDGEL <- calcNormFactors(selDGEL)
# selDesign <- model.matrix(~GenoVir + Batch + Sex, data = selGroups)
# selDesign <- model.matrix(~GenoVir + Sex, data = selGroups)
# selDesign <- model.matrix(~Sex + GenoVir, data = selGroups)
selDesign <- model.matrix(~GenoVir, data = selGroups)
selDGEL <- estimateDisp(selDGEL, selDesign)

selFit <- glmFit(selDGEL, selDesign)
selLrt <- glmLRT(selFit, coef = 'GenoVirWT_CTL')
selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame

## average gene expression per sample
avg_exp <- as.data.frame(log2(cpm(selAggMat) + 1))

## merge dge table and avg gene expression
degavg <- merge(avg_exp, selRes, by = "row.names", all.y = TRUE)
row.names(degavg) <- degavg$Row.names
degavg$Row.names <- NULL

## write and save dge table
write.table(degavg, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_TABLE.tsv", sep = ""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
save(degavg, selLrt, file = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG.RData", sep = ""))

## filter dge table
dge_sig_temp <- degavg[degavg$FDR <= fdrcutoff,]
dge_sig <- dge_sig_temp[abs(dge_sig_temp$logFC) >= sprintf("%.2f", log2(fccutoff)),]
dge_sig_up <- row.names(dge_sig[dge_sig$logFC > 0,])
dge_sig_dn <- row.names(dge_sig[dge_sig$logFC < 0,])

## volcano plot
limit_x_max <- round(max(abs(degavg$logFC)) + 1)
limit_x_min <- -1 * limit_x_max
limit_y_max <- round(max(-log10(degavg$FDR)) + 10)
anno_x_max <- limit_x_max - 0.5
anno_x_min <- limit_x_min + 0.5
anno_y_max <- limit_y_max - 5

pvol1 <- ggplot(data = degavg, aes(x = logFC, y = -log10(FDR))) +
            geom_point(colour="darkgrey", size=1.2, shape=16, alpha=0.6) +
            geom_point(data=dge_sig[dge_sig$logFC > 0,], colour="green3", size=1.2, shape=16, alpha=0.8) +
            geom_point(data=dge_sig[dge_sig$logFC < 0,], colour="royalblue", size=1.2, shape=16, alpha=0.8) +
            annotate(geom = "text", x = anno_x_max, y = anno_y_max, label = length(dge_sig_up)) +
            annotate(geom = "text", x = anno_x_min, y = anno_y_max, label = length(dge_sig_dn)) +
            theme_bw() +
            geom_vline(xintercept = log2(fccutoff), colour = "red", linetype = 2) + geom_vline(xintercept = -1 * log2(fccutoff), colour = "red", linetype = 2) +
            geom_hline(yintercept = -log10(fdrcutoff), colour = "red", linetype = 2) +
            xlim(limit_x_min, limit_x_max) +
            ylim(0, limit_y_max) +
            labs(x="Avg log Fold Change", y="- log10 FDR")
ggsave(filename = paste("NK_DGE_PSEUDOBULK_EDGER_", mygroup1, "_x_", mygroup2, "_DGE_AVG_VOLCANO.PDF", sep = ""), plot = pvol1, width = 4, height = 4, units = "in", dpi = 300, useDingbats = FALSE)

```

