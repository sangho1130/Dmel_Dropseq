library('SCENIC'); packageVersion('SCENIC')
library('RcisTarget'); packageVersion('RcisTarget')

expr <- read.delim('merged.3tps.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
label <- read.delim('merged.3tps.clustering_labels.txt', header = T, sep = '\t', check.names = F, row.names = 1)
colnames(label) <- c('Timepoint', 'Celltype', 'Subcluster')

label <- subset(label, Subcluster == 'PSC' | Subcluster == 'PH1' | Subcluster == 'PH2')
head(label); nrow(label) #345

expr <- expr[, rownames(label)]
identical(colnames(expr), rownames(label)); dim(expr)

org <- 'dmel'
dbDir <- 'databases' # RcisTarget databases location
myDatasetTitle='PSC_PH1-2'
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=20)

label$nCount_RNA <- as.numeric(label$nCount_RNA)
label$nFeature_RNA <- as.numeric(label$nFeature_RNA)
label$percent.mt <- as.numeric(label$percent.mt)

dir.create("int")
saveRDS(label, file="int/cellInfo.Rds")

scenicOptions@inputDatasetInfo$label <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network ###
nCellsPerGene <- apply(expr, 1, function(x) sum(x>0))
nCountsPerGene <- apply(expr, 1, sum)
summary(nCellsPerGene)
summary(nCountsPerGene)
max(expr) # 1001
sum(expr>0) / sum(expr==0) # 0.08411523

# 345*0.025 = 8.625
# 345*0.05 = 17.25
# First filter:
minReads <- 3*0.025*ncol(expr) # 25.875
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads) # 4677 genes
# Second filter:
minSamples <- ncol(expr)*0.025 #8.625
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells) # 4588 genes

motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases) # 
interestingGenes <- c("Dl")
interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]


genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
expr_filtered <- expr[genesKept, ]
expr_filtered <- as.matrix(expr_filtered)

### Correlation ###
corrMat <- cor(t(expr_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
save.image('tmp1.Rdata')


### GENIE3 ###
#load('tmp1.Rdata')
expr_filtered <- log2(expr_filtered + 1)
runGenie3(as.matrix(expr_filtered), scenicOptions)
save.image('tmp2.Rdata')

### runSCENIC ###
#load('tmp2.Rdata')
expr <- as.matrix(expr)
expr <- log2(expr+1)
dim(expr) # 
#scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose #<- TRUE
scenicOptions@settings$seed <- 1907151
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, expr)
save.image('tmp3.Rdata')


#load('tmp3.Rdata')
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(15, 20, 25, 30, 35, 40, 45, 50), perpl=c(10, 20, 30, 40, 50))#, seed = 19070611)
pdf('int/tSNE_compare.pdf')
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="Subclustering_v2", cex=.5) # saving in pdf files
dev.off()
pdf('int/tSNE_compare.oneshot.pdf', width = 30, height = 22)
par(mfcol=c(5, 8))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="Subclustering_v2", cex=1)
dev.off()

# Using only "high-confidence" regulons (normally similar)
fileNames.highconf <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(15, 20, 25, 30, 35, 40, 45, 50), perpl=c(10, 20, 30, 40, 50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
#, seed = 19070651)
pdf('int/tSNE_compare.oHC.pdf')
fileNames.highconf <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames.highconf, scenicOptions, showLegend=T, varName="Subclustering_v2", cex=.5)
dev.off()
pdf('int/tSNE_compare.oHC.oneshot.pdf', width = 30, height = 22)
par(mfcol=c(5, 8))
fileNames.highconf <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames.highconf, scenicOptions, showLegend=T, varName="Subclustering_v2", cex=1)
dev.off()
save.image('tmp4.post_tSNE.Rdata')


scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 45
scenicOptions@settings$defaultTsne$perpl <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

runSCENIC_3_scoreCells(scenicOptions, expr) 
save.image('tmp4.post_tSNE.scoreCells.Rdata')


######################
### Local computer ###
######################
library('SCENIC'); packageVersion('SCENIC')
library('RcisTarget'); packageVersion('RcisTarget')

load('tmp4.post_tSNE.scoreCells.Rdata')
scenicOptions@settings$nCores <- 3
exprMat <- expr
aucellApp <- plotTsne_AUCellApp(scenicOptions, expr)
savedSelections <- shiny::runApp(aucellApp)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
#scenicOptions@settings$seed <- 123
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#save.image('tmp4.post_tSNE.scoreCells.newthresholds.Rdata')

# scenicOptions@settings$devType = "png"
load('tmp4.post_tSNE.scoreCells.newthresholds.Rdata')
head(label)

#
addlab <- read.delim('../../Subclustering_post/PSC+PH1-2/PSC+PH1-2.txt', header = T, row.names = 1)
head(addlab)
addlab <- subset(addlab, PSC_SCs_res_0.7 != 'Neurons')
addlab <- droplevels(addlab)
addlab <- addlab[rownames(label), ]
identical(rownames(addlab), rownames(label))

addlab_psc <- subset(addlab, Subclustering == 'PSC')
addlab_psc <- droplevels(addlab_psc)

label$Subclustering_v3 <- as.character(label$Subclustering_v2)
label[rownames(addlab_psc), 'Subclustering_v3'] <- as.character(addlab_psc$PSC_SCs_res_0.7)
unique(label$Subclustering_v3)
label$Subclustering_v3 <- factor(label$Subclustering_v3,
                                 levels = c('PH_1', 'PH_2', 'PSC_1', 'PSC_2', 'PSC_3'))
label$Subclustering_v2 <- NULL
#

label <- droplevels(label)
head(label)
saveRDS(label, file="int/cellInfo.Rds")

runSCENIC_4_aucell_binarize(scenicOptions)
save.image('tmp5.threshold_updates.ver1.Rdata')


