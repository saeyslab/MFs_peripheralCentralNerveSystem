library("limma")
library("oligo")
library("mogene10stv1cdf")
library("mogene10sttranscriptcluster.db")
library("pd.mogene.1.0.st.v1")
library("genefilter")
library("annotate")
library("ggplot2")
library("rgl")
library("reshape")
library("pheatmap")
library("sva")

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

celFilesOverview <- read.table(file="documentation/neededMicroarrays_nerve.txt", sep="\t", 
                               stringsAsFactors=FALSE, comment.char="", header=TRUE)

### Load raw data
pathMA<-"/path/to/Exp2227_RawData/"
celFiles<-paste(pathMA, celFilesOverview[,1],sep="")

rawData<-read.celfiles(celFiles)

### Create ph
ph<-pData(rawData)
check<-cbind(rownames(ph),celFilesOverview[,1])
ph<-cbind(ph, celFilesOverview[,-1])

########################################
##### NORMALIZE DATA
########################################
dataRMA<-rma(rawData)
dataRMAmatrix<-exprs(dataRMA)


########################################
##### FILTER DATA
########################################
##Create mogene1map
websiteAnnot<-read.table("/path/to/annotationMoGene1-0.txt",stringsAsFactors=FALSE)
colnames(websiteAnnot)<-c('probe_id','symbol')
mogene1map <- sapply(websiteAnnot$symbol,as.character)
names(mogene1map)<-websiteAnnot$probe_id

setdiff(rownames(dataRMAmatrix),names(mogene1map))

##Filter 
geneNames<-unlist(lapply(rownames(dataRMAmatrix),function(x) mogene1map[[x]]))
t=aggregate(dataRMAmatrix,list(Gene = geneNames),max)
newexp<-t[-1,]
rownames(newexp)<-newexp[,1]
newexp<-newexp[,-1]
check<-cbind(colnames(newexp),rownames(ph), ph$fullName)
colnames(newexp)<-ph$fullName

##Reorder
expTable<-newexp
dim(expTable)
# 22320    25

################################################################################
########## PCA
################################################################################
library("rgl")

toPlot<-expTable

### Calculate variance
variance<-apply(toPlot, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-toPlot[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-ph$PCAcolor


### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab="pc1", ylab="pc2", zlab="pc3")
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-3), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.3)
grid3d(c("x+", "y-", "z-"))

### Save
rgl.viewpoint(0, 5)
rgl.viewpoint(-30, 5)



################################################################################
########## GET DE GENES
################################################################################

#### Create design matrix
TS<-ph$subgroup
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

#### Fit linear model on data
fit <- lmFit(expTable, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=MHCII_D1_5-ResMF_StSt, 
                             group2=MHCII_D5-ResMF_StSt, 
                             group3=MHCII_D10-ResMF_StSt, 
                             group4=ResMF_D37-ResMF_StSt, 
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Quick list of DE genes
summa.fit <- decideTests(fit.eb, p.value = 0.01, lfc = 1)
summary(summa.fit)
# group1 group2 group3 group4
# Down      622    303    144     55
# NotSig  21034  21516  21941  22116
# Up        664    501    235    149

summa.fit <- decideTests(fit.eb, p.value = 0.01, lfc = 2)
summary(summa.fit)
# group1 group2 group3 group4
# Down      158     49     24     10
# NotSig  22017  22185  22252  22276
# Up        145     86     44     34



########################################
##### 1. ResMF_StSt vs MHCII_D1.5
########################################
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.01,2)
dim(DEgenesGroup1)
##303


########################################
##### 2. ResMF_StSt vs MHCII_D5
########################################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.01,2)
dim(DEgenesGroup2)
##135


########################################
##### 3. ResMF_StSt vs MHCII_D10
########################################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.01,2)
dim(DEgenesGroup3)
##68


########################################
##### 4. ResMF_StSt vs ResMF_D37
########################################
allGenesGroup4<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=4)
DEgenesGroup4<-getDEgenes(allGenesGroup4,0.01,2)
dim(DEgenesGroup4)
##44



########################################
##### Write results
########################################
allDEgenes<-c(rownames(DEgenesGroup1), rownames(DEgenesGroup2), rownames(DEgenesGroup3), rownames(DEgenesGroup4))
length(allDEgenes)
# 550


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

DEgenesGroup1<-DEgenesGroup1[setdiff(rownames(DEgenesGroup1), dissociationGenes),]
rownames(DEgenesGroup1)[grep('Gm|Rik',rownames(DEgenesGroup1))]
DEgenesGroup1<-DEgenesGroup1[-grep('Gm|Rik',rownames(DEgenesGroup1)),]

DEgenesGroup2<-DEgenesGroup2[setdiff(rownames(DEgenesGroup2), dissociationGenes),]
rownames(DEgenesGroup2)[grep('Gm|Rik',rownames(DEgenesGroup2))]
DEgenesGroup2<-DEgenesGroup2[-grep('Gm|Rik',rownames(DEgenesGroup2)),]

DEgenesGroup3<-DEgenesGroup3[setdiff(rownames(DEgenesGroup3), dissociationGenes),]
rownames(DEgenesGroup3)[grep('Gm|Rik',rownames(DEgenesGroup3))]
DEgenesGroup3<-DEgenesGroup3[-grep('Gm|Rik',rownames(DEgenesGroup3)),]

DEgenesGroup4<-DEgenesGroup4[setdiff(rownames(DEgenesGroup4), dissociationGenes),]
rownames(DEgenesGroup4)[grep('Gm|Rik',rownames(DEgenesGroup4))]
DEgenesGroup4<-DEgenesGroup4[-grep('Gm|Rik',rownames(DEgenesGroup4)),]



listDEgenes<-list('StSt_vs_D1.5'=DEgenesGroup1,'StSt_vs_D5'=DEgenesGroup2,'StSt_vs_D10'=DEgenesGroup3,'StSt_vs_D37'=DEgenesGroup4)
lapply(listDEgenes,function(x){dim(x)})


################################################################################
########## CREATE HEATMAP
################################################################################

##################################################
##### Preparation
##################################################

##### Get genes #####
upGenes<-unlist(lapply(listDEgenes, function(x){c(head(rownames(x), 10))}))
names(upGenes)<-NULL
sum(duplicated(upGenes))


downGenes<-unlist(lapply(listDEgenes, function(x){c(tail(rownames(x), 10))}))
names(downGenes)<-NULL
sum(duplicated(downGenes))

wantedGenes<-c(upGenes[!duplicated(upGenes)], downGenes[! duplicated(downGenes)])
length(wantedGenes)


##### Get samples #####
wantedSamples<-colnames(expTable)
expProfiles<-expTable[wantedGenes,wantedSamples]

##### Normalize #####
expProfilesNorm<-normalizePerGene(expProfiles)

gem<-apply(expProfiles,1,mean)
expProfilesNorm2<-expProfiles-gem

toPlot<-expProfilesNorm2
gapsCol<-c(3,9,15,22)


#### HEATMAP TO FIND NICE ORDER ####
library("gplots")
myColors<-colorRampPalette(c("royalblue","skyblue3", "gray88", "indianred1","indianred3"))
tmp<-heatmap.2(as.matrix(toPlot), col=myColors, scale="none", density.info="none", trace="none", key=TRUE, Colv=FALSE, dendrogram="row", cexCol=0.5, cexRow=0.6)

wantedGenes<-allDEgenes[tmp$rowInd]



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

##test heatmap
pheatmap(as.matrix(toPlot),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 15, cellheight = 15, gaps_col = gapsCol, fontsize=15,
         show_rownames = T, show_colnames = T, breaks = myBreaks)


##small heatmap
pheatmap(as.matrix(toPlot),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 15, cellheight = 15, gaps_col = gapsCol, fontsize=15,
         show_rownames = T, show_colnames = T, breaks = myBreaks, 
         filename = "heatmap.pdf", width = 11.48, height = 15.82)












