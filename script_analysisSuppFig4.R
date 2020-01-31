library('Seurat')
library('ggplot2')
library('gridExtra')

########## Function #########
addDashes<-function(myVar){
  myVar[grep('K10_',myVar)]<-paste0(myVar[grep('K10_',myVar)],"-1")
  myVar[grep('K11_',myVar)]<-paste0(myVar[grep('K11_',myVar)],"-2")
  myVar[grep('K22_',myVar)]<-paste0(myVar[grep('K22_',myVar)],"-3")
  myVar[grep('K23_',myVar)]<-paste0(myVar[grep('K23_',myVar)],"-4")
  return(myVar)
}


################################################################################
########## LOAD SEURAT OBJECT
################################################################################

### Load seuratObj Van Hove
load(file='neededData/seuratObjMF.Robj')
seuratObjKia<-UpdateSeuratObject(seuratObjMFAnew)
rm(seuratObjMFAnew)
DimPlot(seuratObjKia, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
allCellsKia<-colnames(seuratObjKia)
length(allCellsKia)
# 10863

### Load seuratObj
seuratObjMartin<-readRDS(file='neededData/seuratObj.rds')
seuratObjMartin<-UpdateSeuratObject(seuratObjMartin)
DimPlot(seuratObjMartin, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
allCellsMartin<-colnames(seuratObjMartin)
length(allCellsMartin)
# 691

################################################################################
########## MERGE DATA
################################################################################

normDataKia<-seuratObjKia[['RNA']]@data
dim(normDataKia)
# 13518 10863

normDataMartin<-seuratObjMartin[['RNA']]@data
dim(normDataMartin)
# 10449   691
colnames(normDataMartin)<-paste0('Martin_',colnames(normDataMartin))

########## Perpare for merge ##########
### Get overlap
overlapGenes<-intersect(rownames(normDataKia), rownames(normDataMartin))
length(overlapGenes)
# 10340

normDataOverlap<-cbind(normDataKia[overlapGenes,], normDataMartin[overlapGenes,])
dim(normDataOverlap)
# 10340 11554

### Get genes of sample1
genesOnlyKia<-setdiff(rownames(normDataKia), rownames(normDataMartin))
length(genesOnlyKia)
# 3178

tmp<-matrix(c(0),nrow=length(genesOnlyKia), ncol=ncol(normDataMartin))
rownames(tmp)<-genesOnlyKia
colnames(tmp)<-colnames(normDataMartin)
normDataOnlyKia<-cbind(as.matrix(normDataKia[genesOnlyKia,]), tmp)
dim(normDataOnlyKia)
# 3178 11554

### Get genes of sample2
genesOnlyMartin<-setdiff(rownames(normDataMartin), rownames(normDataKia))
length(genesOnlyMartin)
# 109

tmp<-matrix(c(0),nrow=length(genesOnlyMartin), ncol=ncol(normDataKia))
rownames(tmp)<-genesOnlyMartin
colnames(tmp)<-colnames(normDataKia)
normDataOnlyMartin<-cbind(tmp, as.matrix(normDataMartin[genesOnlyMartin,]))
dim(normDataOnlyMartin)
# 109 11554

########## Do merge ##########
mergedNormData<-rbind(as.matrix(normDataOverlap), normDataOnlyKia, normDataOnlyMartin)
dim(mergedNormData)
# 13627 11554


########## Quantile normalization ##########
library('preprocessCore')

normData<-normalize.quantiles(mergedNormData)
dim(normData)
# 13627 11554

rownames(normData)<-rownames(mergedNormData)
colnames(normData)<-colnames(mergedNormData)


###################################################################
########## CREATE SEURAT OBJECT
################################################################################
listLabels<-list('K10','K11','K22','K23','Martin')

##### Create object #####
seuratObj <- CreateSeuratObject(counts = normData, project = "seuratObj", min.cells = 3, min.features = 200)
dim(seuratObj)
# 13544 11554

seuratObj[['RNA']]@counts[1:5,1:5]


################################################################################
########## FILTER DATA
################################################################################

seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0("results_quantile_final/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
for(i in 1:length(listLabels)){
  toSearch<-paste0('-',i)
  metaDataTable[grep(toSearch,rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-listLabels[[i]]
}
metaDataTable$source<-"Kia"
metaDataTable[metaDataTable$orig.ident=="Martin","source"]<-"Martin"

seuratObj@meta.data<-metaDataTable
head(seuratObj@meta.data)
table(seuratObj@meta.data$orig.ident)
table(seuratObj@meta.data$source)


################################################################################
########## NORMALIZE
################################################################################
seuratObj[['RNA']]@data<-mergedNormData

### Get normalised values
seuratObj[['RNA']]@data[1:5,1:5]


##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})
mitoGenes<-rownames(seuratObj)[grep('^mt-',rownames(seuratObj))]
metaDataTable$percent.mito<-colSums(as.matrix(seuratObj[['RNA']]@data[mitoGenes,]))/metaDataTable$nUMI*100

drawVlnPlotSeurat_split(metaDataTable, paste0("results_quantile_final/5_afterNorm_splitted.png"))

################################################################################
########## GET HVG
################################################################################

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj))

### Get more info about HVGs (mean, dispersion and dispersion scaled)
head(HVFInfo(seuratObj))


### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# png(file=paste0(sampleFolder,"results_quantile_final/clean/6_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
# dev.off()


################################################################################
########## SCALE DATA
################################################################################
# Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
# Scaling data => Shifts the expression of each gene, so that the mean expression across cells is 0
# Scaling data => Scales the expression of each gene, so that the variance across cells is 1 (so that highly-expressed genes do not dominate)

seuratObj <- ScaleData(object = seuratObj)

### Get scaled values
seuratObj[['RNA']]@scale.data[1:5,1:5]

##### Check per group #####
head(metaDataTable)
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@scale.data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@scale.data),2,function(x){sum(x>0)})
tmp<-intersect(mitoGenes, rownames(seuratObj[['RNA']]@scale.data))
metaDataTable$percent.mito<-colSums(as.matrix(seuratObj[['RNA']]@data[tmp,]))/metaDataTable$nUMI*100

drawVlnPlotSeurat_split(metaDataTable, paste0("results_quantile_final/7_afterScale_splitted.png"))

################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

names(seuratObj)
# "RNA" "pca"
seuratObj[['pca']]@cell.embeddings[1:5,1:5]

########################################
########## PCA PLOT
########################################
pdf(file=paste0("results_quantile_final/8a_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3), group.by = "orig.ident")
dev.off()

########################################
########## HEATMAP OF PCs
########################################

### Create heatmap of PC 1-40
pdf(file=paste0("results_quantile_final/9a_selectPC.pdf"))
PCHeatmap(seuratObj, dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

### Create PCElbowplot
png(file=paste0("results_quantile_final/9b_selectPC.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################
# dimsToTry<-c(15,20,25,35)
# resToUse<-0.8

### Final
dimsToTry<-c(20)
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse)
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0("results_quantile_final/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30)
  umapPlot<-DimPlot(seuratObj, reduction.use = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction.use = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0("results_quantile_final/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
  
  ##### Rename clusters
  ### Rename clusters
  seuratObjTmp<-seuratObj
  neededCells<-WhichCells(seuratObjKia,idents = 0)
  Idents(seuratObjTmp, cells=neededCells)<-"microglia"
  
  neededCells<-WhichCells(seuratObjKia,idents = 1)
  Idents(seuratObjTmp, cells=neededCells)<-"cpHiBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 4)
  Idents(seuratObjTmp, cells=neededCells)<-"SdBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 5)
  Idents(seuratObjTmp, cells=neededCells)<-"DloBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 6)
  Idents(seuratObjTmp, cells=neededCells)<-"DhiBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 7)
  Idents(seuratObjTmp, cells=neededCells)<-"cpLoBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 8)
  Idents(seuratObjTmp, cells=neededCells)<-"cpEpiBAM"
  
  neededCells<-paste0('Martin_',WhichCells(seuratObjMartin,idents = 0))
  Idents(seuratObjTmp, cells=neededCells)<-"endo"
  
  neededCells<-paste0('Martin_',WhichCells(seuratObjMartin,idents = 1))
  Idents(seuratObjTmp, cells=neededCells)<-"epi"
  
  p1<-DimPlot(seuratObjTmp, reduction.use = "umap", label = T, label.size = 8)
  p2<-DimPlot(seuratObjTmp, reduction.use = "umap", label = T, label.size = 8, group.by = 'orig.ident')
  
  ggsave(grid.arrange(p1, p2, ncol=2),
         file=paste0("results_quantile_final/10c_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
}

saveRDS(seuratObj, file="Robjects/seuratObjFirst.rds")

seuratObjFirst<-seuratObj


################################################################################
########## DETECT THE CELLS AFFECTED BY DISSOCIATION
################################################################################
###Remove ActD genes. From Van Hove et al. Nature Neuroscience. 2019
DEgenes_exp1<-readRDS(file="/home/liesbetm/Documents/Projects/extern_collab/KiavashMovahedi/singleCell_exp2836_revision/JP30-JP31/Robjects/listDEgenes.rds")
DEgenes_exp2<-readRDS(file="/home/liesbetm/Documents/Projects/extern_collab/KiavashMovahedi/singleCell_exp2836_revision/JP32-JP33/Robjects/listDEgenes.rds")
tmp<-DEgenes_exp1$JP31_vs_JP30
DEgenes_exp1<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])
tmp<-DEgenes_exp2$JP33_vs_JP32
DEgenes_exp2<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])
dissociationGenes<-intersect(DEgenes_exp1,DEgenes_exp2)
length(dissociationGenes)
# 224

dissocationGenes<-dissociationGenes
normData<-seuratObj[['RNA']]@data

setdiff(dissocationGenes, rownames(normData))
dissocationGenes<-intersect(dissocationGenes, rownames(normData))
setdiff(dissocationGenes, rownames(normData))

########## Start the code from the paper ##########
Selection<-normData[dissocationGenes,]
Selection[is.na(Selection)]<-0

#Generate  a  new dataframe  (DataPercentages)  containing  for  each  cell:  t-SNE  dimensions  ("V1"  and "V2"), 
#sum of reads from all dissociation-affected genes ("Sums") and percentage of transcriptome of that cell that maps to 
#dissociation affected reads ("Percentage"; this equals "Sums"-column divided by total read count of that cell):
tSNECoordinates <-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)
head(cbind(rownames(tSNECoordinates),colnames(normData)))
tail(cbind(rownames(tSNECoordinates),colnames(normData)))

Sums <-colSums(Selection)
TotalSums <-colSums(normData)
DataPercentages <-merge(tSNECoordinates, Sums, by="row.names", all=T)
row.names(DataPercentages) <-DataPercentages$Row.names
DataPercentages <-DataPercentages[,-1]
colnames(DataPercentages)[3] <-"Sums"
DataPercentages$Percentage <-DataPercentages$Sums*100/TotalSums

#Make a histogram showing the distribution of the metric "Percentage of dissociation-affected reads per cell":
png(file="results_quantile_final/13a_histogramDissociationAffectedCells.png")
hist(DataPercentages$Percentage, breaks = 100, col = "lightgrey", main = "Expression level dissociation-affected genes", 
     xlab = "Sum expression level of dissociation-affected genes", ylab = "Number of cells")
abline(v=5.8,col="blue",lwd=2,lty=2)
dev.off()

#Based  on  this  histogram,  choose  an  appropriate  cutoff  value.  All  cells  with  a  percentage  equal  to  or above 
#this value are affected by the dissociation procedure, will be labelled as "Dissociation-affected" in the "DataPercentages" dataframe:
SetCutoff <-5.8
DataPercentages$Dissociation_affected <-ifelse(DataPercentages$Percentage >= SetCutoff,1,0)

#Calculate  what percentage of the cells  will be annotated as dissociation-affected cells  with this cutoff value:
PercentageDissociationAffected   <-round((nrow(DataPercentages[DataPercentages$Percentage   >= SetCutoff, ]))/(nrow(DataPercentages))*100, digits = 2)
print(c("Percentage of cells annotated as dissociation-affected cells is", PercentageDissociationAffected))
##87.67

cellsAffected<-rownames(DataPercentages[DataPercentages$Dissociation_affected == "1",])
length(cellsAffected)
# 10129
# saveRDS(cellsAffected, file="Robjects/cellsAffected_quantile_final.rds")

### Adjust cellsAffected
tmpCells<-setdiff(colnames(seuratObj)[grep('Martin',colnames(seuratObj))], cellsAffected[grep('Martin',cellsAffected)])
cellsAffected<-c(cellsAffected, tmpCells)
length(cellsAffected)
# 10200
saveRDS(cellsAffected, file="Robjects/cellsAffected_quantile_final.rds")

########## Color the affected cells on the tSNE ##########
###F8766D=red
###00BFC4=cyan
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

tmp<-clusterMatrix
tmp$affected<-"FALSE"
tmp[cellsAffected,which(colnames(tmp)=="affected")]<-"TRUE"

p <- ggplot()+
  geom_point(aes(x=UMAP_1,y=UMAP_2, colour=tmp$affected), data=umapTable, size=2, shape=20) +
  scale_color_manual(values=c("FALSE"="#00BFC4", "TRUE"="#F8766D")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
print(p)
ggsave(p, file="results_quantile_final/13b_colorDissociationAffectedCells_adjusted.png")



FeaturePlot(object = seuratObj, features = "Fos", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
FeaturePlot(object = seuratObj, features = "Dusp1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)




################################################################################
########## 2. CREATE SEURAT OBJECT
################################################################################
listLabels<-list('K10','K11','K22','K23','Martin')

##### Create object #####
seuratObj <- CreateSeuratObject(counts = normData, project = "seuratObj", min.cells = 3, min.features = 200)
dim(seuratObj)
# 13544 11554

seuratObj[['RNA']]@counts[1:5,1:5]


################################################################################
########## FILTER DATA
################################################################################

seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0("results_quantile_final/clean/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
for(i in 1:length(listLabels)){
  toSearch<-paste0('-',i)
  metaDataTable[grep(toSearch,rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-listLabels[[i]]
}
##Add source
metaDataTable$source<-"Kia"
metaDataTable[metaDataTable$orig.ident=="Martin","source"]<-"Martin"
##Add affected
metaDataTable$affected<-'no'
metaDataTable[cellsAffected, 'affected']<-'yes'

seuratObj@meta.data<-metaDataTable
head(seuratObj@meta.data)
table(seuratObj@meta.data$orig.ident)
table(seuratObj@meta.data$source)
table(seuratObj@meta.data$affected)

################################################################################
########## NORMALIZE
################################################################################
seuratObj[['RNA']]@data<-mergedNormData

### Get normalised values
seuratObj[['RNA']]@data[1:5,1:5]


##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})
mitoGenes<-rownames(seuratObj)[grep('^mt-',rownames(seuratObj))]
metaDataTable$percent.mito<-colSums(as.matrix(seuratObj[['RNA']]@data[mitoGenes,]))/metaDataTable$nUMI*100

drawVlnPlotSeurat_split(metaDataTable, paste0("results_quantile_final/clean/5_afterNorm_splitted.png"))

################################################################################
########## GET HVG
################################################################################

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj))

### Get more info about HVGs (mean, dispersion and dispersion scaled)
head(HVFInfo(seuratObj))

HvgGenes<-setdiff(VariableFeatures(seuratObj), dissocationGenes)
VariableFeatures(seuratObj)<-HvgGenes
length(VariableFeatures(seuratObj))
# 1830


### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# png(file=paste0(sampleFolder,"results_quantile_final/clean/6_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
# dev.off()


################################################################################
########## SCALE DATA
################################################################################
# Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
# Scaling data => Shifts the expression of each gene, so that the mean expression across cells is 0
# Scaling data => Scales the expression of each gene, so that the variance across cells is 1 (so that highly-expressed genes do not dominate)

# seuratObj <- ScaleData(object = seuratObj, vars.to.regress = 'affected')
seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c('source','affected'))

### Get scaled values
seuratObj[['RNA']]@scale.data[1:5,1:5]

##### Check per group #####
head(metaDataTable)
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@scale.data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@scale.data),2,function(x){sum(x>0)})
tmp<-intersect(mitoGenes, rownames(seuratObj[['RNA']]@scale.data))
metaDataTable$percent.mito<-colSums(as.matrix(seuratObj[['RNA']]@data[tmp,]))/metaDataTable$nUMI*100

drawVlnPlotSeurat_split(metaDataTable, paste0("results_quantile_final/clean/7_afterScale_splitted.png"))

################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

names(seuratObj)
# "RNA" "pca"
seuratObj[['pca']]@cell.embeddings[1:5,1:5]

########################################
########## PCA PLOT
########################################
pdf(file=paste0("results_quantile_final/clean/8a_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3), group.by = "orig.ident")
dev.off()

########################################
########## HEATMAP OF PCs
########################################

### Create heatmap of PC 1-40
pdf(file=paste0("results_quantile_final/clean/9a_selectPC.pdf"))
PCHeatmap(seuratObj, dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

### Create PCElbowplot
png(file=paste0("results_quantile_final/clean/9b_selectPC.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################
dimsToTry<-c(15,20,25,30,35)
resToUse<-0.8

### Final
dimsToTry<-c(30)
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse)
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0("results_quantile_final/clean/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30)
  umapPlot<-DimPlot(seuratObj, reduction.use = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction.use = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0("results_quantile_final/clean/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
  
  ##### Rename clusters
  ### Rename clusters
  seuratObjTmp<-seuratObj
  neededCells<-WhichCells(seuratObjKia,idents = 0)
  Idents(seuratObjTmp, cells=neededCells)<-"microglia"
  
  neededCells<-WhichCells(seuratObjKia,idents = 1)
  Idents(seuratObjTmp, cells=neededCells)<-"cpHiBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 4)
  Idents(seuratObjTmp, cells=neededCells)<-"SdBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 5)
  Idents(seuratObjTmp, cells=neededCells)<-"DloBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 6)
  Idents(seuratObjTmp, cells=neededCells)<-"DhiBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 7)
  Idents(seuratObjTmp, cells=neededCells)<-"cpLoBAM"
  
  neededCells<-WhichCells(seuratObjKia,idents = 8)
  Idents(seuratObjTmp, cells=neededCells)<-"cpEpiBAM"
  
  neededCells<-paste0('Martin_',WhichCells(seuratObjMartin,idents = 0))
  Idents(seuratObjTmp, cells=neededCells)<-"endo"
  
  neededCells<-paste0('Martin_',WhichCells(seuratObjMartin,idents = 1))
  Idents(seuratObjTmp, cells=neededCells)<-"epi"
  
  p1<-DimPlot(seuratObjTmp, reduction.use = "umap", label = T, label.size = 8)
  p2<-DimPlot(seuratObjTmp, reduction.use = "umap", label = T, label.size = 8, group.by = 'orig.ident')
  
  ggsave(grid.arrange(p1, p2, ncol=2),
         file=paste0("results_quantile_final/clean/10c_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20)
  
}

##### Save object
saveRDS(seuratObj, file=paste0("Robjects/seuratObj_quantile_clean_finalFinal.rds"))


FeaturePlot(object = seuratObj, features = "Dusp1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
FeaturePlot(object = seuratObj, features = "Sall1", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
FeaturePlot(object = seuratObj, features = "Clec10a", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)


################################################################################
########## GET DE GENES
################################################################################
DimPlot(seuratObj, reduction.use = "umap", label = T, label.size = 8)

########################################
##### some clusters vs some clusters
########################################

endo_vs_all <- FindMarkers(seuratObj, ident.1 = c(7), min.pct = 0.10, 
                              min.diff.pct=0.20, logfc.threshold = 0.25, only.pos = TRUE)

epi_vs_all <- FindMarkers(seuratObj, ident.1 = c(8), min.pct = 0.10, 
                              min.diff.pct=0.20, logfc.threshold = 0.30, only.pos = TRUE)

endoEpi_vs_all <- FindMarkers(seuratObj, ident.1 = c(7,8), min.pct = 0.10, 
                                            min.diff.pct=0.25, logfc.threshold = 0.25, only.pos = TRUE)

epi_vs_endo <- FindMarkers(seuratObj, ident.1 = 7, ident.2 = 8, min.pct = 0.10, 
                           min.diff.pct=0.20, logfc.threshold = 0.25, only.pos = TRUE)


endo_vs_cpBam <- FindMarkers(seuratObj, ident.1 = 7, ident.2 = c(1,5), min.pct = 0.10, 
                           min.diff.pct=0.20, logfc.threshold = 0.25, only.pos = TRUE)


endo_vs_dBam <- FindMarkers(seuratObj, ident.1 = 7, ident.2 = c(4,3), min.pct = 0.10, 
                             min.diff.pct=0.20, logfc.threshold = 0.25, only.pos = TRUE)

endo_vs_cpBamDbam <- FindMarkers(seuratObj, ident.1 = 7, ident.2 = c(1,5,4,3), min.pct = 0.10, 
                            min.diff.pct=0.20, logfc.threshold = 0.25, only.pos = TRUE)


##### Create list
listDEgenesGroups<-tibble::lst(endo_vs_all, epi_vs_all, endoEpi_vs_all, epi_vs_endo, endo_vs_cpBam, endo_vs_dBam, endo_vs_cpBamDbam)

##Add geneSymbol in column (for the export)
listDEgenesGroups<-lapply(listDEgenesGroups,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesGroups<-lapply(listDEgenesGroups, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesGroups<-lapply(listDEgenesGroups, function(x){dplyr::mutate(x,'score'=pct.1/(pct.2+0.01)*avg_logFC)})
##Sort on logFC
listDEgenesGroups<-lapply(listDEgenesGroups,function(x){x<-x[order(x$score, decreasing=T),]})


##write to Excel
library('openxlsx')
write.xlsx(listDEgenesGroups, paste0("results_quantile_final/clean_final/markersClusters.xlsx"))


