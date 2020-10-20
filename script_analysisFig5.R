library('Seurat')
library('dplyr')
library('Matrix')
library("scater")
library("scran")
library('gridExtra')

################################################################################
########## LOAD DATA
################################################################################
getwd()

dataAggr <- Read10X("aggr_nerveStStD1D5/filtered_gene_bc_matrices_mex/mm10/")
rawData<-as.matrix(dataAggr)

rawDataExport<-as.data.frame(dataAggr)
rawDataExport$gene<-rownames(rawDataExport)
otherCols<-setdiff(colnames(rawDataExport),'gene')
rawDataNew<-rawDataExport[,c('gene',otherCols)]
colnames(rawDataNew)<-gsub('-1','-StSt',colnames(rawDataNew))
colnames(rawDataNew)<-gsub('-2','-D1',colnames(rawDataNew))
colnames(rawDataNew)<-gsub('-3','-D5',colnames(rawDataNew))

write.table(rawDataNew, file="countTable_aggrNerveStStD1D5.txt", sep="\t", row.names = F, col.names = T)

########################################
##### Characteristics of data
########################################
rawData[1:5,1:5]

dim(rawData)
#27998  5586

sum(rawData==0)
#151520524
# Total dataset: 27998*5586 = 156396828
# % zero: 151520524/156396828*100 = 96.8821

nrCells<-apply(rawData,1,function (x){sum(x>0)})
min(nrCells)
max(nrCells)
##min: 0
##max: 5582
length(nrCells[nrCells<3])
##15584
##27998-15584 = 12414 genes to keep
length(nrCells[nrCells==0])
##13735


nrGenes<-apply(rawData,2,function (x){sum(x>0)})
min(nrGenes)
max(nrGenes)
##min: 27
##max: 2294
length(nrGenes[nrGenes<200])
##6
##5586-6= 5580 cells to keep

rm(nrGenes)
rm(nrCells)
rm(dataAggr)

################################################################################
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
library("scater")
sce<-SingleCellExperiment(list(counts=rawData))
dim(sce)
# 27998  5586

##### Get spike inns #####
is.spike <- grepl("^ERCC", rownames(sce))
sum(is.spike)
##0

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)
##13

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))
# 5586   28
colnames(colData(sce))

##### Create metaData matrix #####
metaData<-data.frame("staticNr"=colnames(rawData),"nGene"=sce$total_features,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

########################################
########## Get outliers
########################################

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop <- isOutlier(sce$total_features, nmads=4, type="lower", log=TRUE)
sum(feature.drop)
#2

##same as UMI in Seurat pipeline
libsize.drop <- isOutlier(sce$total_counts, nmads=4, type="lower", log=TRUE)
sum(libsize.drop)
#0

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=4, type="higher")
sum(mito.drop)
#219

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)
# 27998  5365
###5586-5365 = 221 cells removed


########################################
########## Create violinPlots
########################################

toPlot<-metaData
toPlot<-metaData[! metaData$final.drop,]
toPlot<-metaData[! metaData$pca.drop,]

p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="aggr_nerveStStD1D5/results/plots/1a_beforeFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="aggr_nerveStStD1D5/results/plots/2_afterFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="aggr_nerveStStD1D5/results/plots/3b_beforePcaFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="aggr_nerveStStD1D5/results/plots/4_afterPcaFiltering.png")


########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

##nGene
png(file="aggr_nerveStStD1D5/results/plots/1b_nGene.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file="aggr_nerveStStD1D5/results/plots/1v_nUMI.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

##percent.mito
png(file="aggr_nerveStStD1D5/results/plots/1d_percMito.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
dev.off()


########################################
########## Create PCA
########################################
colnames(colData(sce))
##on position 4, 22 and 18
colnames(colData(sce))[colnames(colData(sce))=="log10_total_counts"]<-"log10_total_counts_feature_controls"
colnames(colData(sce))[colnames(colData(sce))=="pct_counts_feature_control"]<-"pct_counts_feature_controls"
colnames(colData(sce))[colnames(colData(sce))=="total_features_feature_control"]<-"total_features_feature_controls"

png(file="aggr_nerveStStD1D5/results/plots/3a_pca.png",  width = 850, height = 642)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata", detect_outliers=TRUE) + fontsize
dev.off()

##### Detect bad cells #####
sceNew<-runPCA(sce, pca_data_input = "pdata", detect_outliers=TRUE)
table(sceNew$outlier)
# FALSE 
# 5365
outs<-colnames(sceNew)[sceNew$outlier]
### Add to metaData
metaData$pca.drop<-metaData$final.drop
metaData[outs,which(colnames(metaData)=="pca.drop")]<-TRUE

##### Color bad cells on PCA plot #####
colorDF<-as.data.frame(cbind(colnames(sceNew),"1"), stringsAsFactors=F)
rownames(colorDF)<-colorDF[,1]
colorDF[outs,2]<-"2"
colorDF[,2]<-as.factor(colorDF[,2])
tmp<-colorDF[,2,drop=F]

scater::plotPCA(sceNew, pca_data_input="pdata", return_SCESet = FALSE, colour_by=tmp, shape_by=tmp, theme_size=14)

##### Create violin plots again ####
##change color of geom_jitter() into 'pca.drop'

##### Remove the bad cells based on the PCA plot ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
sum(pca.drop)
##0

# sce <- sce[,!(pca.drop)]
dim(sce)
## 27998  5365


rm(sceNew)

################################################################################
########## QC: GENES
################################################################################

ave.counts <- calcAverage(sce)
png(file="aggr_nerveStStD1D5/results/plots/5_QCgenes.png", width = 850, height = 642)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=-3.1,col="blue",lwd=2,lty=2)
dev.off()

cutoff<-10^-3.1
demo.keep <- ave.counts >= cutoff
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
# Mode   FALSE    TRUE 
# logical   15903   12095
sce<-filtered.sce

rm(filtered.sce)
rm(ave.counts)
rm(demo.keep)

########## Keep cells with at least 200 genes
to.keep<-which(sce@colData$total_counts>=200)
length(to.keep)
# 5365 = all the cells


################################################################################
########## NORMALIZATION
################################################################################
library('limSolve')
dim(sce)
# 12095  5365

##### Method 2 #####
# While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero. 
# This manifests as negative size factors, which are obviously nonsensical. To avoid this, the computeSumFactors function will automatically remove 
# low-abundance genes prior to the calculation of size factors. Genes with an average count below a specified threshold (min.mean) are ignored. 
# For read count data, the default value of 1 is usually satisfactory. For UMI data, counts are lower so a threshold of 0.1 is recommended.
sce <- computeSumFactors(sce, min.mean=0.1)
summary(sizeFactors(sce))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04821 0.51456 0.87511 1.00000 1.36615 4.60299

png(file="aggr_nerveStStD1D5/results/plots/6_normalization.png", width = 850, height = 642)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
dev.off()


sce <- normalize(sce)
mat<-exprs(sce)

saveRDS(sce, file="aggr_nerveStStD1D5/Robjects/sce.rds")
saveRDS(mat, file="aggr_nerveStStD1D5/Robjects/mat.rds")


################################################################################
########## CREATE SEURAT OBJECT
################################################################################

dim(mat)
# 12095  5365

##### Create object #####
seuratObj <- CreateSeuratObject(raw.data = mat, min.cells = 3, min.genes = 0, project = "seuratObj")
dim(seuratObj@raw.data)
# 11957  5365
# 138 extra genes removed

# ##### Fill var.genes (HVG genes) #####
# seuratObj@var.genes<-rownames(hvg.out)
# length(seuratObj@var.genes)
# ##992

##### Fill scaled.data (normalized data) #####
seuratObj@data<-mat

################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="aggr_nerveStStD1D5/results/plots/7_afterNormalization.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

################################################################################
########## GET HVG
################################################################################
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.7)

length(seuratObj@var.genes)
# 1081

################################################################################
########## SCALE DATA
################################################################################
seuratObj <- ScaleData(object = seuratObj)


################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 10, pcs.compute = 40)

########################################
########## PCA PLOT
########################################
pdf(file="aggr_nerveStStD1D5/results/plots/8a_PCA.pdf", width = 10)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 3)
dev.off()

########################################
########## HEATMAP OF PCs
########################################
pdf(file="aggr_nerveStStD1D5/results/plots/9a_selectPC.pdf")
PCHeatmap(seuratObj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 37:40, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

PCHeatmap(seuratObj, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

seuratObj <- JackStraw(seuratObj, num.replicate = 100, num.pc = 40)
pdf(file="aggr_nerveStStD1D5/results/plots/9c_jackStrawPlot.pdf")
JackStrawPlot(seuratObj, PCs = 1:40)
dev.off()

png(file="aggr_nerveStStD1D5/results/plots/9b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################

##try PCs 5,6,8,9,10,11,12,15,17,20,25



seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:12, resolution = 0.8, print.output = 0, save.SNN = F)
##### Create tSNE plot
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:12, do.fast = TRUE)
png(file="aggr_nerveStStD1D5/results/plots/10_tSNE_1-12_res08.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)
dev.off()

##### Save object
saveRDS(seuratObj, file="aggr_nerveStStD1D5/Robjects/seuratObj.rds")


################################################################################
########## 2. CREATE SEURAT OBJECT
################################################################################

library('Seurat')
library('dplyr')

getwd()

########## Load first seuratObj ##########
seuratObj<-readRDS(file="aggr_nerveStStD1D5/Robjects/seuratObj.rds")

logTable_v1<-as.matrix(seuratObj@data)
clusterMatrix_v1<-seuratObj@meta.data
tsneTable_v1<-as.data.frame(seuratObj@dr$tsne@cell.embeddings, stringsAsFactors = F)

##### tSNE plot
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)


########## Remove bad cells (clusters 9) ##########
clusterMatrix<-clusterMatrix_v1[! clusterMatrix_v1$res.0.8 == 9,]
dim(clusterMatrix)
# 5284    5

table(clusterMatrix_v1$res.0.8)

logTable<-logTable_v1[,rownames(clusterMatrix)]
dim(logTable)
# 12095  5284

########## Create new SeuratObj ##########
seuratObj <- CreateSeuratObject(raw.data = logTable, min.cells = 3, min.genes = 200, project = "seuratObj")
dim(seuratObj@raw.data)
# 11930  5282

##### Fill scaled.data (normalized data) #####
seuratObj@data<-seuratObj@raw.data

################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/1_afterNormalization.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

################################################################################
########## GET HIGH VARIABLE GENES
################################################################################
png(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/2_HVG.png", width = 850, height = 642)
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.7)
dev.off()

length(seuratObj@var.genes)
# 1106

################################################################################
########## SCALE DATA
################################################################################
seuratObj <- ScaleData(object = seuratObj)



################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 10, pcs.compute = 40)

########################################
########## PCA PLOT
########################################
pdf(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/3_pca.pdf", width=10)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 3)
dev.off()

########################################
########## HEATMAP OF PCs
########################################
pdf(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/4a_selectPC.pdf")
PCHeatmap(seuratObj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 37:40, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

PCHeatmap(seuratObj, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

seuratObj <- JackStraw(seuratObj, num.replicate = 100, num.pc = 40)
pdf(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/4c_jackStrawPlot.pdf")
JackStrawPlot(seuratObj, PCs = 1:40)
dev.off()

png(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/4b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 40)
dev.off()

################################################################################
########## CLUSTER THE CELLS
################################################################################

#try PCs 5,7,8,9,10,11,12,14,16,20,25

seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:11, resolution = 0.6, print.output = 0, 
                          save.SNN = F)

##### Create tSNE plot
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:11, do.fast = TRUE)
png(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/5_tSNE_1-11_res06_FINAL.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)
dev.off()

##### Save object
saveRDS(seuratObj, file="aggr_nerveStStD1D5/Robjects/seuratObj_final.rds")

########## Create UMAP ##########
seuratObj <- RunUMAP(seuratObj, reduction.use = "pca", dims.use = 1:11)

# png(file="aggr_nerveStStD1D5/results/plots/clean_tSNE/5_UMAP_1-11_FINAL.png", width = 850, height = 642)
DimPlot(seuratObj, reduction.use = "umap", do.label = T, label.size = 8)
# dev.off()

umapTable<-as.data.frame(seuratObj@dr$umap@cell.embeddings, stringsAsFactors = F)
saveRDS(seuratObj, file="aggr_nerveStStD1D5/Robjects/seuratObj_final_withUMAP.rds")


################################################################################
########## DE GENES
################################################################################

########################################
##### Focus on new groups
########################################
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)

markers_recruitedMFsD1<-FindMarkers(seuratObj, ident.1 = c(1,2), min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)
markers_recruitedMFsD5<-FindMarkers(seuratObj, ident.1 = c(0,3), min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)
markers_recruitedMFsD1and5<-FindMarkers(seuratObj, ident.1 = 6, min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)
markers_insideD1<-FindMarkers(seuratObj, ident.1 = 7, min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)
markers_insideStStandD5<-FindMarkers(seuratObj, ident.1 = 4, min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)
markers_outside<-FindMarkers(seuratObj, ident.1 = 5, min.pct = 0.10, logfc.threshold = 0.50, only.pos = T)

listDEgenes_groups<-list("recruitedMFsD1"=markers_recruitedMFsD1, "recruitedMFsD5"=markers_recruitedMFsD5, 
                         "recruitedMFsD1andD5"=markers_recruitedMFsD1and5, "insideD1"=markers_insideD1,
                         "insideStStandD5"=markers_insideStStandD5, "outside"=markers_outside)

##Add geneSymbol in column (for the export)
listDEgenes_groups<-lapply(listDEgenes_groups,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenes_groups<-lapply(listDEgenes_groups, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Calculate score
listDEgenes_groups<-lapply(listDEgenes_groups, function(x){dplyr::mutate(x, score=pct.1/pct.2*avg_logFC)})
##Sort on logFC
listDEgenes_groups<-lapply(listDEgenes_groups,function(x){x<-x[order(x$avg_logFC, decreasing=T),]})

saveRDS(listDEgenes_groups, file="aggr_nerveStStD1D5/Robjects/markersGroups.rds")

### Diagnostics
lapply(listDEgenes_groups, function(x){dim(x)})
# $recruitedMFsD1
# [1] 59  7
# 
# $recruitedMFsD5
# [1] 18  7
# 
# $recruitedMFsD1andD5
# [1] 57  7
# 
# $insideD1
# [1] 52  7
# 
# $insideStStandD5
# [1] 165   7
# 
# $outside
# [1] 90  7


################################################################################
########## HEATMAP FINAL GROUPS
################################################################################

myColorPalette<-c("#1c6ab0", "#1f6db1", "#2372b4", "#2675b6", "#2a79b8", "#2e7ebb", "#3181bc", "#3685bf", "#3888c1", "#3d8dc3", 
                  "#4090c5", "#4794c7", "#4c98c9", "#539dcb", "#5ba1ce", "#60a5d0", "#67aad2", "#6cadd4", "#74b2d6", "#79b5d8", 
                  "#81badb", "#8bbfde", "#99c7e2", "#a8d0e6", "#b2d5e9", "#c1dded", "#cbe3f0", "#daebf4", "#e4f0f7", "#f3f8fb", 
                  "#fffefd", "#fff8f5", "#feefe9", "#fee9e1", "#fee0d4", "#fedacc", "#fdd1c0", "#fdcbb8", "#fdc2ac", "#fcb9a0", 
                  "#fcb398", "#fcab8f", "#fca68a", "#fc9e82", "#fc997c", "#fc9174", "#fc8c6e", "#fb8466", "#fb7d5d", "#fb7758", 
                  "#fb7050", "#f96a4d", "#f66348", "#f45e45", "#f15640", "#ee4e3b", "#ec4937", "#ea4133", "#e83c2f", "#e5342a")

##### Needed genes #####
## 3. top 20 DEgenes
###Remove ActD genes. From Van Hove et al. Nature Neuroscience. 2019
dissociationGenes<-read.table(file="neededFiles/dissociationGenes.txt")[,1]
length(dissociationGenes)
# 224

neededGenesTable<-rbind(cbind(listDEgenes_groups$outside,'source'='group6'), 
                        cbind(listDEgenes_groups$insideStStandD5,'source'='group5'), 
                        cbind(listDEgenes_groups$recruitedMFsD1,'source'='group1'), 
                        cbind(listDEgenes_groups$insideD1,'source'='group4'), 
                        cbind(listDEgenes_groups$recruitedMFsD5,'source'='group2'),  
                        cbind(listDEgenes_groups$recruitedMFsD1andD5,'source'='group3'))
neededGenesTable<-as.data.frame(neededGenesTable, stringsAsFactors=F)
neededGenesTable$gene<-as.character(neededGenesTable$gene)
dim(neededGenesTable)
# 441 8
rownames(neededGenesTable)<-1:nrow(neededGenesTable)

### Remove duplicated genes
sum(duplicated(neededGenesTable$gene))
duplicatedGenes<-neededGenesTable$gene[duplicated(neededGenesTable$gene)]
neededGenesTable[neededGenesTable$gene %in% duplicatedGenes,] %>% .[order(.$gene),]
##Keep the gene with highest logFC
rowIDs_toRemove<-c()
for(geneSymbol in unique(duplicatedGenes)){
 tmp<-neededGenesTable[neededGenesTable$gene==geneSymbol,]
 rowIDtoKeep<-rownames(tmp[tmp$avg_logFC == max(tmp$avg_logFC),])
 rowIDs_toRemove<-c(rowIDs_toRemove, setdiff(rownames(tmp),rowIDtoKeep))
 if(length(rowIDtoKeep)>1){
   print(geneSymbol)
 }
}
neededGenesTable<-neededGenesTable[setdiff(rownames(neededGenesTable), rowIDs_toRemove),]
dim(neededGenesTable)
# 380 8

### Get top 20 per group
toPlot<-c()
for(toSearch in c('group6','group5','group1','group4','group2','group3')){
  tmp<-neededGenesTable[neededGenesTable$source==toSearch,]
  myGenes<-setdiff(tmp$gene, dissociationGenes)[1:20]
  print(length(myGenes[!is.na(myGenes)]))
  toPlot<-rbind(toPlot, tmp[tmp$gene %in% myGenes,])
  print(dim(toPlot))
}
dim(toPlot)
# 111   8
toPlot$gene<-as.character(toPlot$gene)
rownames(toPlot)<-1:nrow(toPlot)
sum(duplicated(toPlot$gene))


neededGenes<-toPlot$gene
length(neededGenes)
# 111

tableRes<-table(toPlot$source)
gapsRow<-c(sum(tableRes[1]), sum(tableRes[1:2]), sum(tableRes[1:3]), sum(tableRes[1:4]), sum(tableRes[1:5]))

##### Needed cells #####
cellsRecruitedMFsD1<-rownames(clusterMatrix[clusterMatrix$res.0.6 %in% c(1,2),]) #group1
cellsRecruitedMFsD5<-rownames(clusterMatrix[clusterMatrix$res.0.6 %in% c(0,3),]) #group2
cellsRecruitedMFsD1andD5<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "6",]) #group3
cellsInsideD1<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "7",]) #group4
cellsInsideStStandD5<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "4",]) #group5
cellsOutside<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "5",]) #group6

##order of the cells
# group6 (cl5), group5 (cl4) ,group4 (cl7), group1 (cl1+2),  group2 (cl0+3) en group3 (cl6)
neededCells<-c(cellsOutside, cellsInsideStStandD5, cellsInsideD1, cellsRecruitedMFsD1, cellsRecruitedMFsD5, cellsRecruitedMFsD1andD5)

length(cellsOutside)
# 247
length(cellsInsideStStandD5)
# 518
length(cellsInsideD1)
# 124
length(cellsRecruitedMFsD1)
# 2318
length(cellsRecruitedMFsD5)
# 1857
gapCols<-c(247,765,889,3207,5064)
### gaps: 247 / 247+518= 765 / 765+124= 889 / 889+2318= 3207 / 3207+1857= 5064


##### Prepare heatmap #####
library('pheatmap')
heatmapValues<-seuratObj@scale.data[neededGenes,neededCells]
dim(heatmapValues)
# 111 5282 (top 20 DE genes)

##### Plot heatmap #####
##Test heatmap
pheatmap(t(SCORPIUS::scale_quantile(t(heatmapValues))), cluster_cols=F, cluster_rows=F, gaps_col = gapCols, show_colnames=F,
         fontsize = 5, cellheight=4.9)


