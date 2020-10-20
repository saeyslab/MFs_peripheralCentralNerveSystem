library('Seurat') # Seurat_2.3.2  
library('dplyr')
library('Matrix')
library("scater")
library("scran")
library('gridExtra')

################################################################################
########## LOAD DATA
################################################################################
getwd()

dataStSt <- Read10X("sample_nerveStSt/filtered_gene_bc_matrices/mm10/")
rawData<-as.matrix(dataStSt)

rawDataExport<-as.data.frame(dataStSt)
rawDataExport$gene<-rownames(rawDataExport)
otherCols<-setdiff(colnames(rawDataExport),'gene')
rawDataNew<-rawDataExport[,c('gene',otherCols)]

write.table(rawDataNew, file="countTable_nerveStSt.txt", sep="\t", row.names = F, col.names = T)

########################################
##### Characteristics of data
########################################
rawData[1:5,1:5]

dim(rawData)
#27998  731

sum(rawData==0)
#19950120
# Total dataset: 27998*731 = 20466538
# % zero: 19950120/20466538*100 = 97.47677

nrCells<-apply(rawData,1,function (x){sum(x>0)})
min(nrCells)
max(nrCells)
##min: 0
##max: 730
length(nrCells[nrCells<3])
##18517
##27998-18517 = 9481 genes to keep
length(nrCells[nrCells==0])
##16580


nrGenes<-apply(rawData,2,function (x){sum(x>0)})
min(nrGenes)
max(nrGenes)
##min: 176
##max: 1624
length(nrGenes[nrGenes<200])
##3
##731-3= 728 cells to keep

rm(nrGenes)
rm(nrCells)

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
# 27998  731

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
# 731   28
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
#11

##same as UMI in Seurat pipeline
libsize.drop <- isOutlier(sce$total_counts, nmads=4, type="lower", log=TRUE)
sum(libsize.drop)
#19

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=4, type="higher")
sum(mito.drop)
#20

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
# 27998  691
###731-691 = 40 cells removed


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
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="sample_nerveStSt/results/plots/1a_beforeFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="sample_nerveStSt/results/plots/2_afterFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="sample_nerveStSt/results/plots/3b_beforePcaFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="sample_nerveStSt/results/plots/4_afterPcaFiltering.png")


########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

##nGene
png(file="sample_nerveStSt/results/plots/1b_nGene.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file="sample_nerveStSt/results/plots/1v_nUMI.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

##percent.mito
png(file="sample_nerveStSt/results/plots/1d_percMito.png", width=850)
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

png(file="sample_nerveStSt/results/plots/3a_pca.png",  width = 850, height = 642)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata", detect_outliers=TRUE) + fontsize
dev.off()

##### Detect bad cells #####
sceNew<-runPCA(sce, pca_data_input = "pdata", detect_outliers=TRUE)
table(sceNew$outlier)
# FALSE 
# 691
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
## 27998  691


rm(sceNew)

################################################################################
########## QC: GENES
################################################################################

ave.counts <- calcAverage(sce)
png(file="sample_nerveStSt/results/plots/5_QCgenes.png", width = 850, height = 642)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=-2.7,col="blue",lwd=2,lty=2)
dev.off()

cutoff<-10^-2.7
demo.keep <- ave.counts >= cutoff
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
# Mode   FALSE    TRUE 
# logical   17549   10449
sce<-filtered.sce

rm(filtered.sce)
rm(ave.counts)
rm(demo.keep)

########## Keep cells with at least 200 genes
to.keep<-which(sce@colData$total_counts>=200)
length(to.keep)
# 691 = all the cells


################################################################################
########## NORMALIZATION
################################################################################
library('limSolve')
dim(sce)
# 10449   691

##### Method 2 #####
# While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero. 
# This manifests as negative size factors, which are obviously nonsensical. To avoid this, the computeSumFactors function will automatically remove 
# low-abundance genes prior to the calculation of size factors. Genes with an average count below a specified threshold (min.mean) are ignored. 
# For read count data, the default value of 1 is usually satisfactory. For UMI data, counts are lower so a threshold of 0.1 is recommended.
sce <- computeSumFactors(sce, min.mean=0.1)
summary(sizeFactors(sce))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2130  0.7292  0.9380  1.0000  1.1749  3.4940

png(file="sample_nerveStSt/results/plots/6_normalization.png", width = 850, height = 642)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
dev.off()


sce <- normalize(sce)
mat<-exprs(sce)

saveRDS(sce, file="sample_nerveStSt/Robjects/sce.rds")
saveRDS(mat, file="sample_nerveStSt/Robjects/mat.rds")

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

dim(mat)
# 10449   691

##### Create object #####
seuratObj <- CreateSeuratObject(raw.data = mat, min.cells = 3, min.genes = 0, project = "seuratObj")
dim(seuratObj@raw.data)
# 9437  691
# 1012 extra genes removed

##### Fill scaled.data (normalized data) #####
seuratObj@data<-mat

################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="sample_nerveStSt/results/plots/7_afterNormalization.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

################################################################################
########## GET HVG
################################################################################
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.7)

length(seuratObj@var.genes)
# 1084

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
pdf(file="sample_nerveStSt/results/plots/8a_PCA.pdf", width = 10)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 3)
dev.off()

########################################
########## HEATMAP OF PCs
########################################
pdf(file="sample_nerveStSt/results/plots/9a_selectPC.pdf")
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
pdf(file="sample_nerveStSt/results/plots/9c_jackStrawPlot.pdf")
JackStrawPlot(seuratObj, PCs = 1:40)
dev.off()

png(file="sample_nerveStSt/results/plots/9b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################

##try PCs 3,5,6,8,10

seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = F)
##### Create tSNE plot
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:8, do.fast = TRUE)
png(file="sample_nerveStSt/results/plots/10_tSNE_1-8_res06.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)
dev.off()


##### Save object
saveRDS(seuratObj, file="sample_nerveStSt/Robjects/seuratObj.rds")

png(file="sample_nerveStSt/results/plots/10_tSNE_1-8_res06_FINALFINAL.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 4)
dev.off()




################################################################################
########## DE GENES
################################################################################
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj@dr$tsne@cell.embeddings, stringsAsFactors = F)
logTable<-seuratObj@raw.data

########################################
##### all clusters vs all clusters
########################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10, logfc.threshold = 0.40, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
# 0   1 
# 150 137

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$diff<-tmp$pct.1/tmp$pct.2
  tmp$score<-tmp$diff*tmp$avg_logFC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Remove dissociation genes
cl0_cleaned<-markersList$cluster0[setdiff(markersList$cluster0$gene, dissociationGenes),]
cl1_cleaned<-markersList$cluster1[setdiff(markersList$cluster1$gene, dissociationGenes),]

markersList$cl0_cleaned<-cl0_cleaned
markersList$cl1_cleaned<-cl1_cleaned


################################################################################
########## HEATMAP
################################################################################

###Remove ActD genes. From Van Hove et al. Nature Neuroscience. 2019
dissociationGenes<-read.table(file="neededFiles/dissociationGenes.txt")[,1]
length(dissociationGenes)
# 224


##### Needed genes #####
### 1. DE genes
markersList<-lapply(markersList,function(x){x<-x[order(x$avg_logFC, decreasing=T),]})
neededGenes<-as.character(unlist(lapply(markersList,function(x){x %>% head(.,20) %>% .$gene})))

##remove dissocation genes ##
upGenesCl0<-markersList$cluster0[markersList$cluster0$avg_logFC > 0,]
upGenesCl1<-markersList$cluster1[markersList$cluster1$avg_logFC > 0,]
wantedGenes<-c('Hexb',"Cx3cr1",'Ccl12','Gpr183','Gpr34')
which(upGenesCl0$gene %in% wantedGenes)
# 5  29  69  86 128

upGenesCl0_cleaned<-upGenesCl0[setdiff(upGenesCl0$gene,dissociationGenes),]
upGenesCl1_cleaned<-upGenesCl1[setdiff(upGenesCl1$gene,dissociationGenes),]
for(i in wantedGenes){
  print(paste0(i," on position ",which(upGenesCl0_cleaned$gene==i)))
}


neededGenes<-c(head(upGenesCl0_cleaned$gene,30), head(upGenesCl1_cleaned$gene,30))
cellHeight<-9.5
fontSize<-8
neededGenes<-c(head(upGenesCl0_cleaned$gene,45), head(upGenesCl1_cleaned$gene,45))
cellHeight<-6.3
fontSize<-6


### 2. Genes from bulk heatmap
neededGenes<-c("Mgl2","Fcrls","Gas6","Cbr2",
               "Ccl12","Cx3cr1","Gpr34","Gpr183","Hexb","Mef2c","Serpinf1","St3gal6","Tagap","Trem2",
               "Ccr2","Cxcl1","Selm","Pla2g2d","Qpct","Tnfsf9","Xist","Il1rl1",
               "Bin1","Cd209d","Fxyd2","Stxbp6","Ccl8","Tslp","Cd209a","Mmp9")

gapsRow<-c(which(neededGenes=="Cbr2"), which(neededGenes=="Trem2"), which(neededGenes=="Il1rl1"))
cellHeight<-15
fontSize<-13


### 3. Dissociation genes
neededGenes<-intersect(dissociationGenes, rownames(logTable))
length(neededGenes)
# 219

tmp<-logTable[neededGenes,]
myRes<-cbind(apply(tmp[,cellsInside],1,function(x){sum(x>0)/length(cellsInside)}),
             apply(tmp[,cellsOutside],1,function(x){sum(x>0)/length(cellsOutside)}))
myRes_v2<-myRes>0.10
myRes_v2<-cbind(myRes_v2,'final'=myRes_v2[,1] | myRes_v2[,2])
neededGenes<-rownames(myRes_v2[myRes_v2[,3]==TRUE,])
length(neededGenes)
# 145

cellHeight<-3.7
fontSize<-4.3



##### Needed cells #####
cellsInside<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "0",])
cellsOutside<-rownames(clusterMatrix[clusterMatrix$res.0.6 == "1",])

neededCells<-c(cellsInside,cellsOutside)

##### Create heatmap #####
library('pheatmap')
heatmapValues<-seuratObj@data[neededGenes,neededCells]
dim(heatmapValues)
# 90 691
# 29 691
# 219 691

length(cellsInside)
# 502
### gaps: 502



##Create heatmap
pheatmap(t(SCORPIUS::scale_quantile(t(heatmapValues))), cluster_cols=F, cluster_rows=F, gaps_col = c(502), gaps_row = gapsRow,
         show_colnames=F, cellheight=cellHeight, fontsize_row = fontSize)


### For heatmap with the genes from the bulk heatmap
##sort on expression in WT group
normValues<-t(SCORPIUS::scale_quantile(t(heatmapValues)))
dim(normValues)

valuesInside<-apply(normValues[,cellsInside],1,mean)
valuesOutside<-apply(normValues[,cellsOutside],1,mean)
normValuesMean<-cbind(valuesInside, valuesOutside)
normValuesMean<-normValuesMean[order(normValuesMean[,1], decreasing = T),]
neededGenes<-rownames(normValuesMean)
heatmapValues<-seuratObj@scale.data[neededGenes,neededCells]
dim(heatmapValues)


