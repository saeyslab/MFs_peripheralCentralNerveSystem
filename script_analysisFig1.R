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

celFilesOverview <- read.table(file="documentation/neededMicroarrays_compareImmgen.txt", sep="\t", 
                               stringsAsFactors=FALSE, comment.char="", header=TRUE)

pathMAImmgen<-"/path/to/IMMGEN_MICROARRAYS/"
pathMA1<-"/path/to/Exp1833_and_1890_RawData/"
pathMA1<-"/path/to/Exp2227_RawData/"
pathMA3<-"/path/to/Exp2290_RawData/"
celFilesImmgen<-paste(pathMAImmgen, celFilesOverview[celFilesOverview$source=="1",1],sep="")
celFiles1<-paste(pathMA1, celFilesOverview[celFilesOverview$source=="2",1],sep="")
celFiles2<-paste(pathMA2, celFilesOverview[celFilesOverview$source=="3",1],sep="")
celFiles3<-paste(pathMA3, celFilesOverview[celFilesOverview$source=="4",1],sep="")
celFiles<-c(celFilesImmgen,celFiles1,celFiles2,celFiles3)

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

########################################
##### COMBAT
########################################
combatM<-ph[,c(grep('batch',colnames(ph)), grep('subgroups',colnames(ph)))]
rownames(combatM)<-colnames(newexp)
colnames(combatM)<-c("batch","subgroups")
combatM

modcombat<-model.matrix(~subgroups, data=combatM)
expTable_combat<-ComBat(dat=newexp, batch=as.numeric(combatM$batch), mod=modcombat, par.prior=TRUE)

### To check combat results
expTableAll<-expTable_combat
phAll<-ph


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
PCAcolors<-phAll$PCAcolor


### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab="pc1", ylab="pc2", zlab="pc3")
# text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-3), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.3)
grid3d(c("x+", "y-", "z-"))


rgl.viewpoint(0, 5)
rgl.viewpoint(-30, 5)



################################################################################
########## GET DE GENES
################################################################################

#### Create design matrix: version 1
TS<-c(rep("microgliaResMF",3),rep("theRest",9),rep("microgliaResMF",3),rep("theRest",3))
TS <- factor(TS, levels=unique(TS))
design1 <- model.matrix(~0+TS)
colnames(design1) <- levels(TS)
design1

#### Create design matrix: version 2
TS<-c(rep("theRest",12),rep("resMF",3),rep("theRest",3))
TS <- factor(TS, levels=unique(TS))
design2 <- model.matrix(~0+TS)
colnames(design2) <- levels(TS)
design2

########################################
##### 1. (resMF+microglia) vs theRest
########################################

#### Fit linear model on data
fit <- lmFit(expTable, design1)
cont.matrix <- makeContrasts(microgliaResMFvsTheRest=microgliaResMF-theRest,levels=design1)
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
options(digits=3)
tTall_part1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf)

#### Identification upregulated and downregulated genes
DEgenes_part1_lowS<-getDEgenes(tTall_part1,0.05,1)
DEgenes_part1_midS<-getDEgenes(tTall_part1,0.01,1)
DEgenes_part1_highS<-getDEgenes(tTall_part1,0.01,2)

dim(DEgenes_part1_lowS)
##728
dim(DEgenes_part1_midS)
##325
dim(DEgenes_part1_highS)
##103



########################################
##### 2. resMF vs (theRest+microglia)
########################################

#### Fit linear model on data
fit <- lmFit(expTable, design2)
cont.matrix <- makeContrasts(resMFvsTheRest=resMF-theRest,levels=design2)
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
options(digits=3)
tTall_part2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf)

#### Identification upregulated and downregulated genes
DEgenes_part2_lowS<-getDEgenes(tTall_part2,0.05,1)
DEgenes_part2_midS<-getDEgenes(tTall_part2,0.01,1)
DEgenes_part2_highS<-getDEgenes(tTall_part2,0.01,2)


dim(DEgenes_part2_lowS)
##367
dim(DEgenes_part2_midS)
##204
dim(DEgenes_part2_highS)
##84

########################################
##### resMF vs microglia (for triwise)
########################################

#### Fit linear model on data
fit <- lmFit(expTable, design)
cont.matrix <- makeContrasts(resMFvsMicroglia=resMF-microglia,levels=design)
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
options(digits=3)
tTall_part1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf)

#### Identification upregulated and downregulated genes
DEgenes_part1_lowS<-getDEgenes(tTall_part1,0.05,1)
DEgenes_part1_midS<-getDEgenes(tTall_part1,0.01,1)
DEgenes_part1_highS<-getDEgenes(tTall_part1,0.01,2)

dim(DEgenes_part1_midS)
##859


########################################
##### resMF vs theRest (for triwise)
########################################

#### Fit linear model on data
fit <- lmFit(expTable, design)
cont.matrix <- makeContrasts(resMFvsTheRest=resMF-theRest,levels=design)
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
options(digits=3)
tTall_part2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf)

#### Identification upregulated and downregulated genes
DEgenes_part2_lowS<-getDEgenes(tTall_part2,0.05,1)
DEgenes_part2_midS<-getDEgenes(tTall_part2,0.01,1)
DEgenes_part2_highS<-getDEgenes(tTall_part2,0.01,2)

dim(DEgenes_part2_midS)
##360

########################################
##### microglia vs theRest (for triwise)
########################################

#### Fit linear model on data
fit <- lmFit(expTable, design)
cont.matrix <- makeContrasts(microgliavsTheRest=microglia-theRest,levels=design)
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
options(digits=3)
tTall_part3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf)

#### Identification upregulated and downregulated genes
DEgenes_part3_lowS<-getDEgenes(tTall_part3,0.05,1)
DEgenes_part3_midS<-getDEgenes(tTall_part3,0.01,1)
DEgenes_part3_highS<-getDEgenes(tTall_part3,0.01,2)

dim(DEgenes_part3_midS)
##1447



######################################################################
########## SQLDF
######################################################################

##################################################
########## Preparation
##################################################

########## Create mean expTable ##########
colsMicroglia<-grep('Microglia',colnames(expTable))
colsLung<-grep('MF.Lu',colnames(expTable))
colsSpleen<-grep('RP.Sp',colnames(expTable))
colsPC<-grep('480hi.PC',colnames(expTable))
colsResMF<-grep('ResMF',colnames(expTable))
colsKC<-grep('KC_DTR',colnames(expTable))

microglia_mean<-apply(expTable[,colsMicroglia],1,mean)
lung_mean<-apply(expTable[,colsLung],1,mean)
spleen_mean<-apply(expTable[,colsSpleen],1,mean)
pc_mean<-apply(expTable[,colsPC],1,mean)
resMF_mean<-apply(expTable[,colsResMF],1,mean)
KC_mean<-apply(expTable[,colsKC],1,mean)
microgliaResMF_mean<-apply(expTable[,c(colsMicroglia,colsResMF)],1,median)

expTableMean<-cbind(microglia_mean,lung_mean,spleen_mean,pc_mean,resMF_mean,KC_mean, microgliaResMF_mean)

########## Create data frame for sqldf ##########
library("sqldf")

tmp<-as.data.frame(expTableMean)
tmp$GeneSymbol<-rownames(tmp)
expTableMeanSQLDF<-tmp

##################################################
########## Searching genes
##################################################

### Looking for geneProfile 1: genes up in microglia and resMF
theQuery<-sqldf("select GeneSymbol,microglia_mean,lung_mean,spleen_mean,pc_mean,resMF_mean,KC_mean from expTableMeanSQLDF where 
                microglia_mean < resMF_mean+3.2 and
                microglia_mean > resMF_mean-3.2 and
                lung_mean < microgliaResMF_mean-1.5 and
                spleen_mean < microgliaResMF_mean-1.5 and
                pc_mean < microgliaResMF_mean-1.5 and
                KC_mean < microgliaResMF_mean-1.5")
genesCluster1<-as.character(theQuery$GeneSymbol)
genesCluster1

### Looking for geneProfile 2: genes up only in resMF
theQuery<-sqldf("select GeneSymbol,microglia_mean,lung_mean,spleen_mean,pc_mean,resMF_mean,KC_mean from expTableMeanSQLDF where 
                microglia_mean < resMF_mean-1.5 and
                lung_mean < resMF_mean-1.5 and
                spleen_mean < resMF_mean-1.5 and
                pc_mean < resMF_mean-1.5 and
                KC_mean < resMF_mean-1.5")
genesCluster2<-as.character(theQuery$GeneSymbol)
genesCluster2


### Looking for geneProfile 4: genes up only in microglia
theQuery<-sqldf("select GeneSymbol,microglia_mean,lung_mean,spleen_mean,pc_mean,resMF_mean,KC_mean from expTableMeanSQLDF where 
                resMF_mean < microglia_mean-2 and
                lung_mean < microglia_mean-2 and
                spleen_mean < microglia_mean-2 and
                pc_mean < microglia_mean-2 and
                KC_mean < microglia_mean-2")
genesCluster4<-as.character(theQuery$GeneSymbol)
genesCluster4

##Check overlap: all genes better be with genesCluster1, except for Rasal3 and Rtn4rl1
intersect(genesCluster4, c(genesCluster2))
intersect(genesCluster4, c(genesCluster1))
genesCluster1<-genesCluster1[- grep('Rasal3|Rtn4rl1',genesCluster1)]
genesCluster4<-setdiff(genesCluster4, genesCluster1)


################################################################################
########## CREATE HEATMAP
################################################################################

##################################################
##### Preparation
##################################################

##### Needed genes #####
wantedGenes<-c(genesCluster1 , genesCluster2, genesCluster4)
##Change order
#Hexb in in genesCluster4, but move to genesCluster1 and put together with Sema4d, Cxxc5, Trem2
genesCluster1<-genesCluster1[- which(genesCluster1 %in% c('Sema4d','Cxxc5','Trem2'))]
genesCluster1<-c(genesCluster1,c('Sema4d','Cxxc5','Trem2','Hexb'))
genesCluster4<-genesCluster4[-which(genesCluster4=="Hexb")]
wantedGenes<-c(genesCluster1 , genesCluster2, genesCluster4)

sum(duplicated(wantedGenes))

### Remove ActD genes. From Van Hove et al. Nature Neuroscience. 2019
dissociationGenes<-read.table(file="neededFiles/dissociationGenes.txt")[,1]
length(dissociationGenes)
# 224

wantedGenes<-setdiff(wantedGenes, dissociationGenes)
length(wantedGenes)
# 131

### Remove Gm... and ...Rik genes
badIDs<-grep('Gm|Rik',wantedGenes)
wantedGenes[badIDs]
wantedGenes<-wantedGenes[- badIDs]
length(wantedGenes)
# 121
sum(duplicated(wantedGenes))

### Make sure genes are in correct order
tmp<-c(intersect(genesCluster1, wantedGenes), intersect(genesCluster2, wantedGenes), intersect(genesCluster4, wantedGenes))
wantedGenes<-tmp
length(wantedGenes)
# 121

length(intersect(wantedGenes,genesCluster1))
# 29
length(intersect(wantedGenes,genesCluster2))
# 52
length(intersect(wantedGenes,genesCluster4))
# 40

gapsRow<-c(length(intersect(wantedGenes, genesCluster1)),
           length(intersect(wantedGenes, genesCluster1))+length(intersect(wantedGenes, genesCluster2)))
cellSize<-4.5
fontSize<-5
fileName<-"heatmap_geneProfile1and2and4"


##### Needed samples #####
wantedSamples<-colnames(expTable)[c(colsMicroglia, colsResMF, colsLung, colsSpleen, colsPC, colsKC)]
expProfiles<-expTable[wantedGenes,wantedSamples]

##### Normalize #####
expProfilesNorm<-normalizePerGene(expProfiles)

gem<-apply(expProfiles,1,mean)
expProfilesNorm2<-expProfiles-gem


### Prepare for heatmap
toPlot<-expProfilesNorm2



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
pheatmap(as.matrix(toPlot),color=myColorPalette,scale="none", cluster_rows=T, cluster_cols=F,
         border_color = "gray25", cellwidth = cellSize, cellheight = cellSize, gaps_col = c(6), gaps_row = gapsRow, fontsize=fontSize,
         show_rownames = T, show_colnames = T, breaks = myBreaks)

##nice heatmap
pheatmap(as.matrix(toPlot),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = cellSize, cellheight = cellSize, gaps_col = c(6), gaps_row = gapsRow, fontsize=fontSize,
         show_rownames = T, show_colnames = T, breaks = myBreaks,
         filename = paste0(fileName,".pdf"))




################################################################################
########## TRIWISE PLOT
################################################################################

allDEgenes<-unique(c(rownames(DEgenes_part1_midS),rownames(DEgenes_part2_midS),rownames(DEgenes_part3_midS)))
length(allDEgenes)
##1918

######################################################################
########## PREPARATION
######################################################################
wantedColors<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red")

##### Triwise plots 
microglia_mean<-apply(expTable[,1:3],1,mean)
resMF_mean<-apply(expTable[,13:15],1,mean)
theRest_mean<-apply(expTable[,c(4:12,16:18)],1,mean)

expTable_mean<-cbind(microglia_mean,resMF_mean,theRest_mean)
colnames(expTable_mean)<-c('microglia','ResMfStSt','theRest')


######################################################################
########## TRIWISE PLOTS
######################################################################
theBarycoords<-transformBarycentric(expTable_mean)

p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors)
print(p)

p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes)
print(p)

###Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE)
print(p)









