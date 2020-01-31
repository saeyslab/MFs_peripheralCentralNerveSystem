library("limma")
library("edgeR")
library("ggplot2")

# devtools::install_github("saeyslab/triwise")
library("triwise")
library("htmlwidgets")

##### SN = sciatic nerve macrophages (4 replicates)
##### ON = optic nerve microglia (4 replicates)
##### SPF = brain microglia (5 replicates)

########################################
##### Functions
########################################

###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}

################################################################################
######### LOAD DATA
################################################################################

getwd()

##### Load raw counts
countData<-read.table(file="data_Duitsland/rawData/counts_combined_SN_ON_MG.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
dim(countData)
# 48679    15

##Clean table
rownames(countData)<-countData[,1]
countData<-countData[, - c(1,2)]
dim(countData)
# 48679    13

##Remove ERCC genes
length(grep('ERCC-',rownames(countData)))
# 92
countData<-countData[- grep('ERCC-',rownames(countData)),]
dim(countData)
# 48587    13

##Rename columns
newColnames<-c("SN_1","SN_2","SN_3","SN_4","ON_1","ON_2","ON_3","ON_4","SPF_1","SPF_2","SPF_3","SPF_4","SPF_5")
cbind(colnames(countData),newColnames)
colnames(countData)<-newColnames

### Load meta data
colData<-read.table(file="data_Duitsland/rawData/metadata.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE)
rownames(colData)<-colData$fileName
colData<-colData[,-1]
dim(colData)
# 13  2

### Reorder
countData<-countData[,rownames(colData)]
dim(countData)
# 48587    13

countData_tmp<-countData
countData_tmp$Gene<-rownames(countData_tmp)
otherCols<-setdiff(colnames(countData_tmp),"Gene")
countData_tmp<-countData_tmp[,c('Gene',otherCols)]
write.table(countData_tmp, file="data_Duitsland/rawData/rawData_forUploadNCBI.txt", sep="\t", row.names = F, col.names = T)

################################################################################
########## CREATE OBJECT
################################################################################

y <- DGEList(counts = countData)

################################################################################
########## FILTER DATA
################################################################################

##### Filter low count genes
## always work with count-per-million (CPM) instead of raw counts
## Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library
## Imagine that the lowest lib size is around 6 million reads => threshold is set on CPM>1
## But what if lib size is 20 million? Here CPM of 1 means 20 counts. Do we still use 5 counts and then the CPM cut off
## will be 0.25 or is the threshold of 5 counts a cut off for lib sizes of around 5-6 million? Then we need to put the
## cut off on 20 for lib sizes around 20 million and so use a CPM of 1.
## Threshold needs to be true in at least x samples. x is always the lowest number of replicates.
## for example: 3 samples with each 2 replicates => x set on 2
## => This ensures that a gene will be retained if it is only expressed in both replicates of a certain group

## Do filtering
yNoFilter<-y
myCpm<-cpm(y)

keep = rowSums(cpm(y)>2) >= 4
y = y[keep,]
dim(y)
##11185
dim(yNoFilter)
##48587

##### Reset lib sizes
y$samples$lib.size = colSums(y$counts)
y$samples

################################################################################
########## NORMALISATION
################################################################################

##### Scale normalisation
yNoNormalisation<-y
y <- calcNormFactors(y)

##### MDS-plot
theColors<-c(rep("royalblue",4),rep("darkgreen",4),rep("darkorange",5))
cbind(colData, theColors)


plotMDS(y, cex=1.2, col=theColors)


################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
TS <- paste(colData$cell, colData$condition, sep=".")
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
png(file="data_Duitsland/results/meanVariancePlot.png",width = 1515, height = 1138)
v <- voom(y, design, plot = TRUE)
dev.off()

expTable<-v$E

##### Normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)


################################################################################
########## PCA
################################################################################
library("rgl")

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

PoV <- pca$sdev^2/sum(pca$sdev^2)
summary(pca)

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.5)

### Save as image
rgl.viewpoint(0, 0)
rgl.viewpoint(35, 0)


################################################################################
########## GET DE GENES
################################################################################

#### Fit linear model on data
fit <- lmFit(v, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=mac.sciaticNerve-mac.opticNerve, 
                             group2=mac.sciaticNerve-mac.brain, 
                             group3=mac.opticNerve-mac.brain, 
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Quick list of DE genes
summa.fit <- decideTests(fit.eb, p.value = 0.01, lfc = 1)
summary(summa.fit)
# group1 group2 group3
# Down      704   1674   1108
# NotSig   9541   7461   8709
# Up        940   2050   1368


#### Quick list of DE genes via treat
# When there is a lot of differential expression, sometimes we may want to cut-off on a fold change threshold as well as a p-value threshold 
# so that we follow up on the most biologically significant genes. However, it is not recommended to simply rank by p-value and then discard 
# genes with small logFC’s, as this has been shown to increase the false discovery rate. In other words, you are not controlling the false 
# discovery rate at 5% any more. There is a function called treat in the limma package that performs this style of analysis correctly 
# (McCarthy and Smyth 2009). treat will simply take our fit.cont object, as well as a user-specified log fold change cut-off, and recalculate 
# the moderated t-statistics and p-values with the new information about logFC.
tfit <- treat(fit2, lfc=1)
dt <- decideTests(tfit, p.value = 0.01, lfc = 1)
summary(dt)
# group1 group2 group3
# Down      284    600    129
# NotSig  10455   9420  10443
# Up        446   1165    613


volcanoplot(fit.eb,coef=1,highlight=100)
volcanoplot(fit.eb,coef=2,highlight=100)
volcanoplot(fit.eb,coef=3,highlight=100)

########################################
##### 1. Sciatic vs optic
########################################
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.01,1)
dim(DEgenesGroup1)
##1644

allGenesGroup1_treat<-topTreat(tfit,number=Inf, coef=1)
DEgenesGroup1_treat<-getDEgenes(allGenesGroup1_treat,0.01,1)
dim(DEgenesGroup1_treat)
##730


########################################
##### 2. Sciatic vs microglia
########################################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.01,1)
dim(DEgenesGroup2)
##3724

allGenesGroup2_treat<-topTreat(tfit,number=Inf, coef=2)
DEgenesGroup2_treat<-getDEgenes(allGenesGroup2_treat,0.01,1)
dim(DEgenesGroup2_treat)
##1765

########################################
##### 3. Optic vs microglia
########################################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.01,1)
dim(DEgenesGroup3)
##2476

allGenesGroup3_treat<-topTreat(tfit,number=Inf, coef=3)
DEgenesGroup3_treat<-getDEgenes(allGenesGroup3_treat,0.01,1)
dim(DEgenesGroup3_treat)
##742


######################################################################
########## SQLDF
######################################################################

##################################################
########## Preparation
##################################################

########## Create mean expTable ##########
colsSN<-grep('SN_',colnames(expTable))
colsON<-grep('ON_',colnames(expTable))
colsSPF<-grep('SPF_',colnames(expTable))

sn_mean<-apply(expTable[,colsSN],1,mean)
on_mean<-apply(expTable[,colsON],1,mean)
spf_mean<-apply(expTable[,colsSPF],1,mean)
spf_sn_mean<-apply(expTable[,c(colsSPF,colsSN)],1,mean)
spf_on_mean<-apply(expTable[,c(colsSPF,colsON)],1,mean)
sn_on_mean<-apply(expTable[,c(colsSN,colsON)],1,mean)

expTableMean<-cbind(spf_mean,on_mean,sn_mean,spf_sn_mean,spf_on_mean,sn_on_mean)

########## Create data frame for sqldf ##########
library("sqldf")

tmp<-as.data.frame(expTableMean)
tmp$GeneSymbol<-rownames(tmp)
expTableMeanSQLDF<-tmp

##### Only search in the DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1_treat),rownames(DEgenesGroup2_treat),rownames(DEgenesGroup3_treat)))
length(allDEgenes)
##1983


expTableMeanSQLDF<-expTableMeanSQLDF[allDEgenes,]
dim(expTableMeanSQLDF)
# 1983    4

##################################################
########## Searching genes
##################################################

### Looking for geneProfile 1: genes up only in SPF
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where 
                spf_mean > on_mean+2.5 and
                spf_mean > sn_mean+2.5")
genesCluster1<-as.character(theQuery$GeneSymbol)
length(genesCluster1)
# 101


### Looking for geneProfile 2: genes up only in ON
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where 
                on_mean > spf_mean+2.5 and
                on_mean > sn_mean+2.5")
genesCluster2<-as.character(theQuery$GeneSymbol)
length(genesCluster2)
# 42


### Looking for geneProfile 3: genes up only in SN
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where 
                sn_mean > spf_mean+2.5 and
                sn_mean > on_mean+2.5")
genesCluster3<-as.character(theQuery$GeneSymbol)
length(genesCluster3)
# 450


### Looking for geneProfile 4: genes up in SPF and ON
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where
                spf_mean < on_mean+0.5 and
                spf_mean > on_mean-0.5 and
                spf_mean > sn_mean+2.5")
genesCluster4<-as.character(theQuery$GeneSymbol)
length(genesCluster4)
# 79
genesCluster4[genesCluster4=="Sema4g"]


theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean,spf_on_mean from expTableMeanSQLDF where 
                spf_mean < on_mean+4 and
                spf_mean > on_mean-4 and
                spf_on_mean > sn_mean+5")
genesCluster4_bis<-as.character(theQuery$GeneSymbol)
length(genesCluster4_bis)
# 161
intersect(genesCluster4, genesCluster4_bis)
genesCluster4_bis<-setdiff(genesCluster4_bis, genesCluster4)
length(genesCluster4_bis)
# 128

### Looking for geneProfile 5: genes up in ON and SN
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where 
                sn_mean < on_mean+0.5 and
                sn_mean > on_mean-0.5 and
                sn_mean > spf_mean+2.5")
genesCluster5<-as.character(theQuery$GeneSymbol)
length(genesCluster5)
# 137


### Looking for geneProfile 6: genes up in SPF and SN
theQuery<-sqldf("select GeneSymbol,spf_mean,on_mean,sn_mean from expTableMeanSQLDF where 
                spf_mean < sn_mean+0.5 and
                spf_mean > sn_mean-0.5 and
                spf_mean > on_mean+2.5")
genesCluster6<-as.character(theQuery$GeneSymbol)
length(genesCluster6)
# 4



################################################################################
########## CREATE HEATMAP
################################################################################

colsSN<-grep('SN_',colnames(expTable))
colsON<-grep('ON_',colnames(expTable))
colsSPF<-grep('SPF_',colnames(expTable))

##################################################
##### Preparation
##################################################
wantedGenes<-c(genesCluster1, genesCluster2, genesCluster3, genesCluster4, genesCluster4_bis, genesCluster5, genesCluster6)
length(wantedGenes)
# 941

### Remove ActD genes. From Van Hove et al. Nature Neuroscience. 2019
DEgenes_exp1<-readRDS(file="/path/to/JP30-JP31/Robjects/listDEgenes.rds")
DEgenes_exp2<-readRDS(file="/path/to/JP32-JP33/Robjects/listDEgenes.rds")
tmp<-DEgenes_exp1$JP31_vs_JP30
DEgenes_exp1<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])
tmp<-DEgenes_exp2$JP33_vs_JP32
DEgenes_exp2<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])
dissociationGenes<-intersect(DEgenes_exp1,DEgenes_exp2)
length(dissociationGenes)
# 224

wantedGenes<-setdiff(wantedGenes, dissociationGenes)
length(wantedGenes)
# 856

wantedSamples<-c(colsSPF, colsON, colsSN)
expProfiles<-expTable[wantedGenes,wantedSamples]

##### Normalize #####
expProfilesNorm<-normalizePerGene(expProfiles)

gem<-apply(expProfiles,1,mean)
expProfilesNorm2<-expProfiles-gem

toPlot<-expProfilesNorm2
gapsCol<-c(length(colsSPF), length(colsSPF)+length(colsON))

##################################################
##### Filter genes
##################################################
#### Only keep the genes with perfect red-blue profile
colsSPF<-grep('SPF_',colnames(toPlot))
colsON<-grep('ON_',colnames(toPlot))
colsSN<-grep('SN_',colnames(toPlot))

##### geneProfile1 #####
goodGenes<-vector()
length(genesCluster1)
# 101
genesCluster1<-setdiff(genesCluster1, dissociationGenes)
length(genesCluster1)
# 101

tmp<-apply(toPlot[genesCluster1,colsSPF],1,function(x){sum(x>0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster1,colsON],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster1,colsSN],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile1<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile1)
# 19

toPlotTmp<-toPlot[wantedGenesProfile1,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsSPF[1]], decreasing = T),]
wantedGenesProfile1<-rownames(toPlotTmpSort)
length(wantedGenesProfile1)
# 19

##### geneProfile2 #####
goodGenes<-vector()
length(genesCluster2)
# 42
genesCluster2<-setdiff(genesCluster2, dissociationGenes)
length(genesCluster2)
# 41

tmp<-apply(toPlot[genesCluster2,colsSPF],1,function(x){sum(x<0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster2,colsON],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster2,colsSN],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile2<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile2)
# 6

toPlotTmp<-toPlot[wantedGenesProfile2,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsON[1]], decreasing = T),]
wantedGenesProfile2<-rownames(toPlotTmpSort)
length(wantedGenesProfile2)
# 6

##### geneProfile3 #####
goodGenes<-vector()
length(genesCluster3)
# 450
genesCluster3<-setdiff(genesCluster3, dissociationGenes)
length(genesCluster3)
# 436

tmp<-apply(toPlot[genesCluster3,colsSPF],1,function(x){sum(x<0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster3,colsON],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster3,colsSN],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile3<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile3)
# 79

toPlotTmp<-toPlot[wantedGenesProfile3,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsSN[1]], decreasing = T),]
wantedGenesProfile3<-rownames(toPlotTmpSort)
length(wantedGenesProfile3)
# 79

##### geneProfile4 #####
goodGenes<-vector()
length(genesCluster4)
# 79
genesCluster4<-setdiff(genesCluster4, dissociationGenes)
length(genesCluster4)
# 78

tmp<-apply(toPlot[genesCluster4,colsSPF],1,function(x){sum(x>0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster4,colsON],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster4,colsSN],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile4<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile4)
# 75

toPlotTmp<-toPlot[wantedGenesProfile4,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsSPF[1]], decreasing = T),]
wantedGenesProfile4<-rownames(toPlotTmpSort)
length(wantedGenesProfile4)
# 75

##### geneProfile4 bis #####
goodGenes<-vector()
length(genesCluster4_bis)
# 128
genesCluster4_bis<-setdiff(genesCluster4_bis, dissociationGenes)
length(genesCluster4_bis)
# 128

tmp<-apply(toPlot[genesCluster4_bis,colsSPF],1,function(x){sum(x>0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster4_bis,colsON],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster4_bis,colsSN],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile4_bis<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile4_bis)
# 91

toPlotTmp<-toPlot[wantedGenesProfile4_bis,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsSPF[1]], decreasing = T),]
wantedGenesProfile4_bis<-rownames(toPlotTmpSort)
length(wantedGenesProfile4_bis)
# 91

##### geneProfile5 #####
goodGenes<-vector()
length(genesCluster5)
# 137
genesCluster5<-setdiff(genesCluster5, dissociationGenes)
length(genesCluster5)
# 85

tmp<-apply(toPlot[genesCluster5,colsSPF],1,function(x){sum(x<0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster5,colsON],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster5,colsSN],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile5<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile5)
# 82

toPlotTmp<-toPlot[wantedGenesProfile5,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsON[1]], decreasing = T),]
wantedGenesProfile5<-rownames(toPlotTmpSort)
length(wantedGenesProfile5)
# 82


##### geneProfile6 #####
goodGenes<-vector()
length(genesCluster6)
# 4
genesCluster6<-setdiff(genesCluster6, dissociationGenes)
length(genesCluster6)
# 4

tmp<-apply(toPlot[genesCluster6,colsSPF],1,function(x){sum(x>0)})
goodGenes<-tmp==length(colsSPF)

tmp<-apply(toPlot[genesCluster6,colsON],1,function(x){sum(x<0)})
goodGenes<-cbind(goodGenes, tmp==length(colsON))

tmp<-apply(toPlot[genesCluster6,colsSN],1,function(x){sum(x>0)})
goodGenes<-cbind(goodGenes, tmp==length(colsSN))

goodGenes<-cbind(goodGenes, apply(goodGenes,1,sum))
wantedGenesProfile6<-rownames(goodGenes[goodGenes[,4] == 3,])
length(wantedGenesProfile6)
# 4

toPlotTmp<-toPlot[wantedGenesProfile6,]
toPlotTmpSort<-toPlotTmp[order(toPlotTmp[,colsSN[1]], decreasing = T),]
wantedGenesProfile6<-rownames(toPlotTmpSort)
length(wantedGenesProfile6)
# 4

### Big heatmap with profile4_bis (version 2)
toPlot<-toPlot[c(head(wantedGenesProfile1, 20), head(wantedGenesProfile2, 20), head(wantedGenesProfile3, 20),
                 head(wantedGenesProfile4, 20), head(wantedGenesProfile4_bis, 20), head(wantedGenesProfile5, 20), 
                 head(wantedGenesProfile6, 20)),]
dim(toPlot)
# 109  13

gap1<-length(head(wantedGenesProfile1, 20))
gap2<-gap1+length(head(wantedGenesProfile2, 20))
gap3<-gap2+length(head(wantedGenesProfile3, 20))
gap4<-gap3+length(head(wantedGenesProfile4, 20))
gap5<-gap4+length(head(wantedGenesProfile4_bis, 20))
gap6<-gap5+length(head(wantedGenesProfile5, 20))
gapsRow<-c(gap1, gap2, gap3, gap4, gap5, gap6)


##################################################
##### Plot heatmap
##################################################
library('pheatmap')
library('grid')

myColorPalette<-c("#08306b", "#08326e", "#083573", "#083876", "#083b7b", "#083d7f", "#084083", "#084387", "#08468c", "#084990", 
                  "#084b94", "#0a4f97", "#0b5199", "#0e559d", "#0f579f", "#125ba2", "#135da4", "#1661a7", "#1864aa", "#1967ad", 
                  "#1c6ab0", "#1f6db1", "#2372b4", "#2675b6", "#2a79b8", "#2e7ebb", "#3181bc", "#3685bf", "#3888c1", "#3d8dc3", 
                  "#4090c5", "#4794c7", "#4c98c9", "#539dcb", "#5ba1ce", "#60a5d0", "#67aad2", "#6cadd4", "#74b2d6", "#79b5d8", 
                  "#81badb", "#8bbfde", "#99c7e2", "#a8d0e6", "#b2d5e9", "#c1dded", "#cbe3f0", "#daebf4", "#e4f0f7", "#f3f8fb", 
                  "#fffefd", "#fff8f5", "#feefe9", "#fee9e1", "#fee0d4", "#fedacc", "#fdd1c0", "#fdcbb8", "#fdc2ac", "#fcb9a0", 
                  "#fcb398", "#fcab8f", "#fca68a", "#fc9e82", "#fc997c", "#fc9174", "#fc8c6e", "#fb8466", "#fb7d5d", "#fb7758", 
                  "#fb7050", "#f96a4d", "#f66348", "#f45e45", "#f15640", "#ee4e3b", "#ec4937", "#ea4133", "#e83c2f", "#e5342a", 
                  "#e32f27", "#dd2c25", "#d92924", "#d32622", "#cd2220", "#c9201f", "#c31d1d", "#bf1a1c", "#b9171a", "#b51419", 
                  "#ae1117", "#a91016", "#a00e15", "#980c14", "#920a13", "#890812", "#840711", "#7b0510", "#75030f", "#6d010e")

paletteLength<-100


########## Create heatmap ##########
myBreaks <- c(seq(min(toPlot), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(toPlot)/paletteLength, max(toPlot), length.out=floor(paletteLength/2)))


pheatmap(as.matrix(toPlot),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F, 
         border_color = "gray25", show_rownames = T, show_colnames = T, breaks = myBreaks,
         fontsize=8, cellwidth = 8.5, cellheight = 8.5, gaps_col = gapsCol, gaps_row = gapsRow, 
         filename = "heatmap.pdf", width = 11.48, height = 15.82)
dev.off()




#######################################################################################################
############################## ------------- TRIWISE PLOT -------------- ##############################
#######################################################################################################

allDEgenes<-unique(c(rownames(DEgenesGroup1_treat),rownames(DEgenesGroup2_treat),rownames(DEgenesGroup3_treat)))
length(allDEgenes)
##1983


######################################################################
########## PREPARATION
######################################################################
wantedColors7<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red",
                 nodiffSet2="limegreen",diffSet2="darkgreen", nodiffSet3="lightblue", diffSet3="darkblue",
                 nodiffSet4="darkorchid1", diffSet4="darkmagenta", nodiffSet5="lightsalmon", diffSet5="darkorange",
                 nodiffSet6="lightskyblue", diffSet6="turquoise4", nodiffSet7="bisque", diffSet7="burlywood4")
wantedColors<-c(nodiffall="gray55",diffall="indianred1")

##### Triwise plots
colsSN<-grep("SN_",colnames(expTable))
colsON<-grep("ON_",colnames(expTable))
colsSPF<-grep("SPF",colnames(expTable))

SN_mean<-apply(expTable[,colsSN],1,mean)
ON_mean<-apply(expTable[,colsON],1,mean)
SPF_mean<-apply(expTable[,colsSPF],1,mean)

expTable_mean<-cbind(SPF_mean,SN_mean,ON_mean)
colnames(expTable_mean)<-c('SPF','SN','ON')

######################################################################
########## TRIWISE PLOTS
######################################################################
theBarycoords<-transformBarycentric(expTable_mean)

###Triwise plot
p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors, showlabels = F)
print(p)

###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)

###Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE)
print(p)



