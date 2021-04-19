# Radiogenomics analysis
# Date: 2021-04-19
# Authors: Qiuchang Sun and Zhi-Cheng Li


# set path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('./radiogenomics_function.R')
import_library()
# ==================== load rna data and clincal data ====================
# read radioGenTrain data and test datasetwd()
options(stringsAsFactors = FALSE);

radioGenTrainData = read.csv(file = "./input/radiogenomics/radioGenTrainData.csv",header = TRUE)
externalTestData = read.csv(file = "./input/radiogenomics/externalTestData.csv",header = TRUE)
tcgaTestData = read.csv(file = "./input/radiogenomics/tcgaTestData.csv",header = TRUE)

radioGenTrainData = log(radioGenTrainData+1)
externalTestData = log(externalTestData+1)
tcgaTestData = log(tcgaTestData+1)

# Os, status, radiomics signature, features
# radioGenTrainData rgtOs rgtSig rgtStatus rgtRadFeature
# externalTestData  etOs etSig etStatus etRadFeature
# tcgaTestData  tcgaOs tcgaStatus

radioGenTrainRad = read.csv("./input/radiogenomics/radioGenTrainRadiomics.csv")
externalTestRad = read.csv("./input/radiogenomics/externalTestRadiomics.csv")
tcgaTestclin = read.csv("./input/radiogenomics/tcgaTestclinical.csv")
rgtOs = radioGenTrainRad$os
rgtStatus = radioGenTrainRad$stauts
rgtSig = radioGenTrainRad$signature
rgtRadFeature = radioGenTrainRad[,4:16]

etOs = externalTestRad$os
etStatus = externalTestRad$stauts
etSig = externalTestRad$signature
etRadFeature = externalTestRad[,4:16]

tcgaOs = tcgaTestclin$os
tcgaStatus = tcgaTestclin$status

# ==================== load geneset ====================
# read radioGenTrain data and test datasetwd()

gmtfile_hallmark = './input/radiogenomics/gmt/h.all.v7.0.symbols.gmt'
gmtfile_biocarta = './input/radiogenomics/gmt/c2.cp.biocarta.v7.0.symbols.gmt'
gmtfile_CPG = './input/radiogenomics/gmt/c2.cgp.v7.0.symbols.gmt'
gmtfile_CP = './input/radiogenomics/gmt/c2.cp.v7.0.symbols.gmt'
gmtfile_PID = './input/radiogenomics/gmt/c2.cp.pid.v7.0.symbols.gmt'
gmtfile_GO = './input/radiogenomics/gmt/c5.bp.v7.0.symbols.gmt'
gmtfile_KEGG = "./input/radiogenomics/gmt/c2.cp.kegg.v7.0.symbols.gmt"
gmtfile_REACTOME = "./input/radiogenomics/gmt/c2.cp.reactome.v7.0.symbols.gmt"
gmtfile_C2ALL = "./input/radiogenomics/gmt/c2.all.v7.0.symbols.gmt"

geneset_hallmark <- read.gmt(gmtfile_hallmark)
geneset_biocarta <- read.gmt(gmtfile_biocarta)
geneset_CPG <- read.gmt(gmtfile_CPG)
geneset_CP <- read.gmt(gmtfile_CP)
geneset_PID <- read.gmt(gmtfile_PID)
geneset_GO <- read.gmt(gmtfile_GO)
geneset_KEGG <- read.gmt(gmtfile_KEGG)
geneset_REACTOME <- read.gmt(gmtfile_REACTOME)

gsva_C2ALL = getGmt(gmtfile_C2ALL)
gsva_GO = getGmt(gmtfile_GO)
gsva_HALL = getGmt(gmtfile_hallmark)

# ==================== radioGenTrain ====================
# wgcna
# for test
load("data.Rdata")
resultDir =  file.path(getwd(),"output/radiogenomics")
if(!(file.exists(resultDir))){
  dir.create(resultDir)}

if(!(file.exists(file.path(resultDir,name = "rgt")))){
  dir.create(file.path(resultDir,name = "rgt"))}
rgtWgcna = wgcnaData(rnaData = G4[1:95,], setPower = TRUE, cutHeight = 70,resultDir = resultDir,name = "rgt")
rgtClust = rgtWgcna[[2]]
table(rgtClust)
rgtData = rgtWgcna[[3]]
rgtPower = rgtWgcna[[4]]
rgtNet = rgtWgcna[[5]]
rgtData = rgtData[rgtNet[["goodGenes"]]]
rgtMEs = rgtWgcna[[6]]
rgtMM = rgtWgcna[[7]]
rgtModuleColor = rgtWgcna[[8]]
rgtModuleColor = rgtModuleColor[rgtNet[["goodGenes"]]]
rgtColorUnique = unique(rgtModuleColor)

# make gsva set
name = "rgt"
rgtGeneName = colnames(rgtData)
rgtGsvaSet = makeGsvaSet(color = rgtColorUnique, moduleColors = rgtModuleColor, geneName = rgtGeneName)
rgtList = rgtGsvaSet[[1]]
rgtGeneSet <- rbind.fill(rgtList)
rgt_gene_set = data.frame(rgtColorUnique,DESCRIPTION = "no",rgtGeneSet)
write.table(file = file.path(resultDir,name,"rgt_gene_set.gmt"), rgt_gene_set, sep = "\t", col.names = F, row.names = F, quote = F)

# calculate gsva score
rgtGsvaSetName = file.path(resultDir,name,"rgt_gene_set.gmt")
rgtGsvaSetGmt = getGmt(rgtGsvaSetName)
rgtDataGsva= as.matrix(t(exp(rgtData)-1))
rgtGsvaResult <- gsva(rgtDataGsva, rgtGsvaSetGmt, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=2)

keepSamples = (rgtClust==1)
rgtOs = rgtOs[keepSamples]
rgtSig = rgtSig[keepSamples]
rgtStatus = rgtStatus[keepSamples]
rgtRadFeature = rgtRadFeature[keepSamples,]
# calculate correaltion between gsva and radiomics siganture
rgtGsvaTraitCor = cor(as.data.frame(t(rgtGsvaResult)), rgtSig, method= 'pearson',use = "p");
rgtGsvaTraitPvalue = corPvalueStudent(rgtGsvaTraitCor, length(rgtSig));
rgtPathGsva = rownames(rgtGsvaTraitPvalue)[rgtGsvaTraitPvalue<0.1]

# enrich correaltion module
length(rgtPathGsva)
name = "rgt"
resultDir =  file.path(getwd(),"output/radiogenomics/gsvaPathway/rgt")
if(!(file.exists(file.path(getwd(),"output/radiogenomics/gsvaPathway")))){
  dir.create(file.path(getwd(),"output/radiogenomics/gsvaPathway"))}
if(!(file.exists(resultDir))){
  dir.create(resultDir)}
rgtGeneSetPathway = genesetPathwayModule(color = rgtPathGsva, moduleColors = rgtModuleColor, geneName = rgtGeneName, MEs = rgtMEs, MM = rgtMM, name = name, resultDir = resultDir)

# calculate gsva score

rgtGsvaC2all <- gsva(rgtDataGsva, gsva_C2ALL, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=2)
rgtGsvaGo <- gsva(rgtDataGsva, gsva_GO, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                      verbose=FALSE, parallel.sz=2)
rgtGsvaHall <- gsva(rgtDataGsva, gsva_HALL, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                      verbose=FALSE, parallel.sz=2)
rgtGsvaAll = rbind(rgtGsvaC2all, rgtGsvaGo, rgtGsvaHall)
rgtGsvaAllPathway = data.frame(rgtGsvaAll)

# found key gene
module = "red"
rgtGeneRed = keyGene(data = rgtData, sig = rgtSig, module = module, moduleColor = rgtModuleColor, MM = rgtMM)
colnames(rgtGeneRed) = c("keyGeneName","keyGeneMM","keyGeneGS","keyGeneGSPvalue")
module = "blue"
rgtGeneBlue = keyGene(data = rgtData, sig = rgtSig, module = module, moduleColor = rgtModuleColor, MM = rgtMM)
colnames(rgtGeneBlue) = c("keyGeneName","keyGeneMM","keyGeneGS","keyGeneGSPvalue")


# ==================== externalTestData ====================
# wgcna
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
setwd(resultDir)
if(!(file.exists(file.path(resultDir,name = "et")))){
  dir.create(file.path(resultDir,name = "et"))}
etWgcna = wgcnaData(rnaData = G4[1:78,], setPower = rgtPower, cutHeight = 85, resultDir = resultDir, name = "et")
etClust = etWgcna[[2]]
table(etClust)
etData = etWgcna[[3]]
etNet = etWgcna[[5]]
etMEs = etWgcna[[6]]
etMM = etWgcna[[7]]
etModuleColor = etWgcna[[8]]
etData = etData[etNet[["goodGenes"]]]
etModuleColor = etModuleColor[etNet[["goodGenes"]]]
etColorUnique = unique(etModuleColor)

# zsummary
etPreservatModule = zsummaryData(dataTrain = rgtData, dataTest = etData, colorTrain = rgtModuleColor, colorTest = etModuleColor,resultDir = resultDir, name = "et")
etPreservatTable = etPreservatModule[[1]]
etPreservatMp = etPreservatModule[[2]]
etZsummary = rownames(etPreservatTable)[which(etPreservatTable$Zsummary.pres>10)]
name="et"
# make gsva set
etGeneName = colnames(etData)
etGsvaSet = makeGsvaSet(color = etColorUnique, moduleColors = etModuleColor, geneName = etGeneName)
etList = etGsvaSet[[1]]
etGeneSet <- rbind.fill(etList)
et_gene_set = data.frame(etColorUnique,DESCRIPTION = "no",etGeneSet)
write.table(file = file.path(resultDir,name,"et_gene_set.gmt"), et_gene_set, sep = "\t", col.names = F, row.names = F, quote = F)

# calculate gsva score
etGsvaSetName = file.path(resultDir,name,"et_gene_set.gmt")
etGsvaSetGmt = getGmt(etGsvaSetName)
etGsvaData= as.matrix(t(etData)-1)
etGsvaResult <- gsva(etGsvaData, etGsvaSetGmt, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=2)

# calculate correaltion between gsva and radiomics siganture
etGsvaTraitCor = cor(as.data.frame(t(etGsvaResult)), etSig, method= 'pearson',use = "p");
etGsvaTraitPvalue = corPvalueStudent(etGsvaTraitCor, length(etSig));
etPathGsva = rownames(etGsvaTraitPvalue)[etGsvaTraitPvalue<0.1]

# enrich correaltion module
length(etPathGsva)
name = "et"
resultDir =  file.path(getwd(),"output/radiogenomics/gsvaPathway/et")
if(!(file.exists(resultDir))){
  dir.create(resultDir)}
etGeneSetPathway = genesetPathwayModule(color = etPathGsva, moduleColors = etModuleColor, geneName = etGeneName, MEs = etMEs, MM = etMM, name = name, resultDir = resultDir)

# calculate gsva score

etGsvaC2all <- gsva(etGsvaData, gsva_C2ALL, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=2)
etGsvaGo <- gsva(etGsvaData, gsva_GO, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                      verbose=FALSE, parallel.sz=2)
etGsvaHall <- gsva(etGsvaData, gsva_HALL, min.sz=1, max.sz=1000, mx.diff=TRUE, 
                      verbose=FALSE, parallel.sz=2)
etGsvaAll = rbind(etGsvaC2all, etGsvaGo, etGsvaHall)
etGsvaAllPathway = data.frame(etGsvaAll)

# found key gene
module = "red"
etGeneRed = keyGene(data = etData, sig = etSig, module = module, moduleColor = etModuleColor, MM = etMM)
colnames(etGeneRed) = c("keyGeneName","keyGeneMM","keyGeneGS","keyGeneGSPvalue")
module = "blue"
etGeneBlue = keyGene(data = etData, sig = etSig, module = module, moduleColor = etModuleColor, MM = etMM)
colnames(etGeneBlue) = c("keyGeneName","keyGeneMM","keyGeneGS","keyGeneGSPvalue")

# ==================== TCGA test set ====================
# read radioGenTrain data and test datasetwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
setwd(resultDir)
if(!(file.exists(file.path(resultDir,name = "tcga")))){
  dir.create(file.path(resultDir,name = "tcga"))}
tcgaWgcna = wgcnaData(rnaData = G4, setPower = TRUE, cutHeight = 90, resultDir = resultDir, name = "tcga")
tcgaClust = tcgaWgcna[[2]]
table(tcgaClust)
tcgaData = tcgaWgcna[[3]]
tcgaNet = tcgaWgcna[[5]]
tcgaMEs = tcgaWgcna[[6]]
tcgaMM = tcgaWgcna[[7]]
tcgaModuleColor = tcgaWgcna[[8]]
tcgaData = tcgaData[tcgaNet[["goodGenes"]]]
tcgaModuleColor = tcgaModuleColor[tcgaNet[["goodGenes"]]]
tcgaColorUnique = unique(tcgaModuleColor)

# zsummary
tcgaPreservatModule = zsummaryData(dataTrain = rgtData, dataTest = tcgaData, colorTrain = rgtModuleColor, colorTest = tcgaModuleColor,resultDir = resultDir, name = "tcga")
tcgaPreservatTable = tcgaPreservatModule[[1]]
tcgaPreservatMp = tcgaPreservatModule[[2]]

# ==================== Heatmap of z-score-normalized gene expressions ====================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
setwd(resultDir)
rgtModuleRed = colnames(rgtData)[rgtModuleColor=="red"]
etModuleRed = colnames(etData)[etModuleColor=="red"]
rgt_et_red = intersect(rgtModuleRed,etModuleRed)
rgtRed = rgtData[,rgt_et_red]
etRed = etData[,rgt_et_red]

rgtModuleBlue = colnames(rgtData)[rgtModuleColor=="blue"]
etModuleBlue = colnames(etData)[etModuleColor=="blue"]
rgt_et_blue = intersect(rgtModuleBlue,etModuleBlue)
rgtBlue = rgtData[,rgt_et_blue]
etBlue = etData[,rgt_et_blue]

rgtGeneExpression = cbind(rgtSig, rgtOs, rgtStatus, rgtRed, rgtBlue)
etGeneExpression = cbind(etSig, etOs, etStatus, etRed, etBlue)
rgtGeneExpression = rgtGeneExpression[order(rgtGeneExpression[,1],decreasing=T),]
etGeneExpression = etGeneExpression[order(etGeneExpression[,1],decreasing=T),]

rgt_gene_heat = heatmapGene(radioGenTrain_heatmap1_gene_expressions = rgtGeneExpression, externalTest_heatmap1_gene_expressions = etGeneExpression, resultDir = resultDir)

# ==================== Heatmap of z-score-normalized gene expressions ====================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
rgtRedPath = file.path(resultDir,'gsvaPathway/rgt/red')
rgtRedPathway = makePathFile(rgtRedPath)
rgtBluePath = file.path(resultDir,'gsvaPathway/rgt/blue')
rgtBluePathway = makePathFile(rgtBluePath)

etRedPath = file.path(resultDir,'gsvaPathway/et/red')
etRedPathway = makePathFile(etRedPath)
etBluePath = file.path(resultDir,'gsvaPathway/et/red')
etBluePathway = makePathFile(etBluePath)

rgt_rt_RedPathway = inner_join(rgtRedPathway,etRedPathway,by=c("ID"="ID"))
rgt_rt_BluePathway = inner_join(rgtBluePathway,etBluePathway,by=c("ID"="ID"))

write.csv(rgt_rt_RedPathway,file="rgt_rt_RedPathway.csv")
write.csv(rgt_rt_BluePathway,file="rgt_rt_BluePathway.csv")

module = read.csv("Module.csv",header = TRUE)
class_1 = module1[which(module1[,6]==1),]
class_2 = module1[which(module1[,6]==2),]
class_3 = module1[which(module1[,6]==3),]
class_4 = module1[which(module1[,6]==4),]

rgtGsvaAllPathway = data.frame(rgtGsvaAll)
rgtClass1Pathway = rgtGsvaAllPathway[match(class_1[,1],rownames(rgtGsvaAllPathway)),]
rgtClass1Pathway = na.omit(rgtClass1Pathway)

rgtClass2Pathway = rgtGsvaAllPathway[match(class_2[,1],rownames(rgtGsvaAllPathway)),]
rgtClass2Pathway = na.omit(rgtClass2Pathway)

rgtClass3Pathway = rgtGsvaAllPathway[match(class_3[,1],rownames(rgtGsvaAllPathway)),]
rgtClass3Pathway = na.omit(rgtClass3Pathway)

rgtClass4Pathway = rgtGsvaAllPathway[match(class_4[,1],rownames(rgtGsvaAllPathway)),]
rgtClass4Pathway = na.omit(rgtClass4Pathway)


rgtClass1RadFeature = as.data.frame(cor(t(rgtClass1Pathway),rgtRadFeature, method= 'pearson'))
rgtClass1RadFeaturePvalue = as.data.frame(corPvalueStudent(as.matrix(rgtClass1RadFeature), nrow(rgtRadFeature)))

rgtClass2RadFeature = as.data.frame(cor(t(rgtClass2Pathway),rgtRadFeature, method= 'pearson'))
rgtClass2RadFeaturePvalue = as.data.frame(corPvalueStudent(as.matrix(rgtClass2RadFeature), nrow(rgtRadFeature)))

rgtClass3RadFeature = as.data.frame(cor(t(rgtClass3Pathway),rgtRadFeature, method= 'pearson'))
rgtClass3RadFeaturePvalue = as.data.frame(corPvalueStudent(as.matrix(rgtClass3RadFeature), nrow(rgtRadFeature)))

rgtClass4RadFeature = as.data.frame(cor(t(rgtClass4Pathway),rgtRadFeature, method= 'pearson'))
rgtClass4RadFeaturePvalue = as.data.frame(corPvalueStudent(as.matrix(rgtClass4RadFeature), nrow(rgtRadFeature)))

list_class_1=list()
list_class = data.frame(feature7 = rownames(rgtClass1RadFeature)[order(abs(rgtClass1RadFeature[,2]),decreasing=TRUE)[1:10]])
list_class_1[[1]] = list_class[c(1,2,3,4,5),]
rgtRadFeature2 = rgtRadFeature[,2] 
rgtClass1Feature2 = rgtClass1Pathway[match(list_class_1[[1]],rownames(rgtClass1Pathway)),]
rgtF1 = cor(rgtRadFeature2,t(rgtClass1Feature2),method= 'pearson')
rgtF1Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF1), length(rgtRadFeature2)))
rgtRadFeature2 = myNormalize(rgtRadFeature2)

list_class = data.frame(feature9 = rownames(rgtClass1RadFeature)[order(abs(rgtClass1RadFeature[,3]),decreasing=TRUE)[1:10]])
list_class_1[[2]] = list_class[c(1,2,3,4,5),]
rgtRadFeature3 = rgtRadFeature[,3]  
rgtClass1Feature3 = rgtClass1Pathway[match(list_class_1[[2]],rownames(rgtClass1Pathway)),]
rgtF3 = cor(rgtRadFeature3,t(rgtClass1Feature3),method= 'pearson')
rgtF3Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF3), length(rgtRadFeature3)))
rgtRadFeature3 = myNormalize(rgtRadFeature3)

list_class = data.frame(feature10 = rownames(rgtClass1RadFeature)[order(abs(rgtClass1RadFeature[,5]),decreasing=TRUE)[1:10]])
list_class_1[[3]] = list_class[c(1,2,3,4,5),]
rgtRadFeature5 = rgtRadFeature[,5] 
rgtClass1Feature5 = rgtClass1Pathway[match(list_class_1[[3]],rownames(rgtClass1Pathway)),]
rgtF5 = cor(rgtRadFeature5,t(rgtClass1Feature5),method= 'pearson')
rgtF5Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF5), length(rgtRadFeature5)))
rgtRadFeature5 = myNormalize(rgtRadFeature5)

list_class = data.frame(feature11 = rownames(rgtClass1RadFeature)[order(abs(rgtClass1RadFeature[,12]),decreasing=TRUE)[1:10]])
list_class_1[[4]] = list_class[c(1,2,3,4,5),]
rgtRadFeature12 = rgtRadFeature[,12] 
rgtClass1Feature12 = rgtClass1Pathway[match(list_class_1[[4]],rownames(rgtClass1Pathway)),]
rgtF12 = cor(rgtRadFeature12,t(rgtClass1Feature12),method= 'pearson')
rgtF12Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF12), length(rgtRadFeature12)))
rgtRadFeature12 = myNormalize(rgtRadFeature12)

list_class = data.frame(feature13 = rownames(rgtClass1RadFeature)[order(abs(rgtClass1RadFeature[,13]),decreasing=TRUE)[1:10]])
list_class_1[[5]] = list_class[c(1,2,3,4,5),]
rgtRadFeature13 = rgtRadFeature[,13]  
rgtClass1Feature13 = rgtClass1Pathway[match(list_class_1[[5]],rownames(rgtClass1Pathway)),]
rgtF13 = cor(rgtRadFeature13,t(rgtClass1Feature13),method= 'pearson')
rgtF13Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF13), length(rgtRadFeature13)))
rgtRadFeature13 = myNormalize(rgtRadFeature13)

list_class_2=list()
list_class = data.frame(feature4 = rownames(rgtClass2RadFeature)[order(abs(rgtClass2RadFeature[,1]),decreasing=TRUE)[1:10]])
list_class_2[[1]] = list_class[c(1,2,3,4,5),]
rgtRadFeature1 = rgtRadFeature[,1]
rgtClass2Feature1 = rgtClass2Pathway[match(list_class_2[[1]],rownames(rgtClass2Pathway)),]
rgtF1 = cor(rgtRadFeature1,t(rgtClass2Feature1),method= 'pearson')
rgtF1value = as.data.frame(corPvalueStudent(as.matrix(rgtF1), length(rgtRadFeature1)))
rgtRadFeature1 = myNormalize(rgtRadFeature1)

list_class = data.frame(feature6 = rownames(rgtClass2RadFeature)[order(abs(rgtClass2RadFeature[,9]),decreasing=TRUE)[1:10]])
list_class_2[[2]] = list_class[c(1,2,3,4,5),]
rgtRadFeature9 = rgtRadFeature[,9] 
rgtClass2Feature9 = rgtClass2Pathway[match(list_class_2[[2]],rownames(rgtClass2Pathway)),]
rgtF9 = cor(rgtRadFeature9,t(rgtClass2Feature9),method= 'pearson')
rgtF9Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF9), length(rgtRadFeature9)))
rgtRadFeature9 = myNormalize(rgtRadFeature9)

list_class_3=list()
list_class = data.frame(feature2 = rownames(rgtClass3RadFeature)[order(abs(rgtClass3RadFeature[,4]),decreasing=TRUE)[1:10]])
list_class_3[[1]] = list_class[c(1,2,3,4,5),]
rgtRadFeature4 = rgtRadFeature[,4]
rgtClass3Feature4 = rgtClass3Pathway[match(list_class_3[[1]],rownames(rgtClass3Pathway)),]
rgtF4 = cor(rgtRadFeature4,t(rgtClass3Feature4),method= 'pearson')
rgtF4Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF4), length(rgtRadFeature4)))
rgtRadFeature4 = myNormalize(rgtRadFeature4)

list_class = data.frame(feature12 = rownames(rgtClass3RadFeature)[order(abs(rgtClass3RadFeature[,6]),decreasing=TRUE)[1:10]])
list_class_3[[2]] = list_class[c(1,2,3,4,5),]
rgtRadFeature6 = rgtRadFeature[,6] 
rgtClass3Feature6 = rgtClass3Pathway[match(list_class_3[[2]],rownames(rgtClass3Pathway)),]
rgtF6 = cor(rgtRadFeature6,t(rgtClass3Feature6),method= 'pearson')
rgtF6Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF6), length(rgtRadFeature6)))
rgtRadFeature6 = myNormalize(rgtRadFeature6)

list_class_4=list()
list_class = data.frame(feature1 = rownames(rgtClass4RadFeature)[order(abs(rgtClass4RadFeature[,7]),decreasing=TRUE)[1:10]])
list_class_4[[1]] = list_class[c(1,2,3,4,5),]
rgtRadFeature7 = rgtRadFeature[,7]
rgtClass4Feature7 = rgtClass4Pathway[match(list_class_4[[1]],rownames(rgtClass4Pathway)),]
rgtF7 = cor(rgtRadFeature7,t(rgtClass4Feature7),method= 'pearson')
rgtF7Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF7), length(rgtRadFeature7)))
rgtRadFeature7 = myNormalize(rgtRadFeature7)

list_class = data.frame(feature3 = rownames(rgtClass4RadFeature)[order(abs(rgtClass4RadFeature[,8]),decreasing=TRUE)[1:10]])
list_class_4[[2]] = list_class[c(1,2,3,4,5),]
rgtRadFeature8 = rgtRadFeature[,8]
rgtClass4Feature8 = rgtClass4Pathway[match(list_class_4[[2]],rownames(rgtClass4Pathway)),]
rgtF8 = cor(rgtRadFeature8,t(rgtClass4Feature8),method= 'pearson')
rgtF8Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF8), length(rgtRadFeature8)))
rgtRadFeature8 = myNormalize(rgtRadFeature8)

list_class = data.frame(feature5 = rownames(rgtClass4RadFeature)[order(abs(rgtClass4RadFeature[,10]),decreasing=TRUE)[1:10]])
list_class_4[[3]] = list_class[c(1,2,3,4,5),]
rgtRadFeature10 = rgtRadFeature[,10] 
rgtClass4Feature10 = rgtClass4Pathway[match(list_class_4[[3]],rownames(rgtClass4Pathway)),]
f10 = cor(rgtRadFeature10,t(rgtClass4Feature10),method= 'pearson')
f10Pvalue = as.data.frame(corPvalueStudent(as.matrix(f10), length(rgtRadFeature10)))
rgtRadFeature10 = myNormalize(rgtRadFeature10)

list_class = data.frame(feature8 = rownames(rgtClass4RadFeature)[order(abs(rgtClass4RadFeature[,11]),decreasing=TRUE)[1:10]])
list_class_4[[4]] = list_class[c(1,2,3,4,5),]
rgtRadFeature11 = rgtRadFeature[,11] 
rgtClass4Feature11 = rgtClass4Pathway[match(list_class_4[[4]],rownames(rgtClass4Pathway)),]
rgtF11 = cor(rgtRadFeature11,t(rgtClass4Feature11),method= 'pearson')
rgtF11Pvalue = as.data.frame(corPvalueStudent(as.matrix(rgtF11), length(rgtRadFeature11)))
rgtRadFeature11 = myNormalize(rgtRadFeature11)


rgtClassFeature = rbind(rgtRadFeature1,rgtRadFeature2,rgtRadFeature3,rgtRadFeature4,rgtRadFeature5,
                        rgtRadFeature6,rgtRadFeature7,rgtRadFeature8,rgtRadFeature9,
                        rgtRadFeature10,rgtRadFeature11,rgtRadFeature12,rgtRadFeature13,
                      rgtClass2Feature1,rgtClass1Feature2,rgtClass1Feature3,rgtClass3Feature4,rgtClass1Feature5,
                      rgtClass3Feature6,rgtClass4Feature7,
                      rgtClass4Feature8,rgtClass2Feature9,
                      rgtClass4Feature10,rgtClass4Feature11,rgtClass1Feature12,rgtClass1Feature13)

rgtClinical = data.frame(rgtSig = rgtSig, rgtOs = rgtOs, rgtStatus = rgtStatus)
rgtPathway = data.frame(rgtClinical, t(rgtClassFeature))
rgtPathway = rgtPathway[order(rgtPathway[,1],decreasing=T),]


etClass1Pathway = etGsvaAllPathway[match(class_1[,1],rownames(etGsvaAllPathway)),]
etClass1Pathway = na.omit(etClass1Pathway)

etClass2Pathway = etGsvaAllPathway[match(class_2[,1],rownames(etGsvaAllPathway)),]
etClass2Pathway = na.omit(etClass2Pathway)

etClass3Pathway = etGsvaAllPathway[match(class_3[,1],rownames(etGsvaAllPathway)),]
etClass3Pathway = na.omit(etClass3Pathway)

etClass4Pathway = etGsvaAllPathway[match(class_4[,1],rownames(etGsvaAllPathway)),]
etClass4Pathway = na.omit(etClass4Pathway)


etRadFeature2 = etRadFeature[,2]
etClass1Feature2 = etClass1Pathway[match(list_class_1[[1]],rownames(etClass1Pathway)),]
etF2 = cor(etRadFeature2,t(etClass1Feature2),method= 'pearson')
etF2Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF2), length(etRadFeature2)))
etRadFeature2 = myNormalize(etRadFeature2)

etRadFeature3 = etRadFeature[,3]
etClass1Feature3 = etClass1Pathway[match(list_class_1[[2]],rownames(etClass1Pathway)),]
etF3 = cor(etRadFeature3,t(etClass1Feature3),method= 'pearson')
etF3Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF3), length(etRadFeature3)))
etRadFeature3 = myNormalize(etRadFeature3)

etRadFeature5 = etRadFeature[,5]
etClass1Feature5 = etClass1Pathway[match(list_class_1[[3]],rownames(etClass1Pathway)),]
etF5 = cor(etRadFeature5,t(etClass1Feature5),method= 'pearson')
etF5Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF5), length(etRadFeature5)))
etRadFeature5 = myNormalize(etRadFeature5)

etRadFeature12 = etRadFeature[,12]
etClass1Feature12 = etClass1Pathway[match(list_class_1[[4]],rownames(etClass1Pathway)),]
etF12 = cor(etRadFeature12,t(etClass1Feature12),method= 'pearson')
etF12Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF12), length(etRadFeature12)))
etRadFeature12 = myNormalize(etRadFeature12)

etRadFeature13 = etRadFeature[,13]
etClass1Feature13 = etClass1Pathway[match(list_class_1[[5]],rownames(etClass1Pathway)),]
etF13 = cor(etRadFeature13,t(etClass1Feature13),method= 'pearson')
etF13Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF13), length(etRadFeature13)))
etRadFeature13 = myNormalize(etRadFeature13)

etRadFeature1 = etRadFeature[,1]
etClass2Feature1 = etClass2Pathway[match(list_class_2[[1]],rownames(etClass2Pathway)),]
etF1 = cor(etRadFeature1,t(etClass2Feature1),method= 'pearson')
etF1Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF1), length(etRadFeature1)))
etRadFeature1 = myNormalize(etRadFeature1)

etRadFeature9 = etRadFeature[,9]
etClass2Feature9 = etClass2Pathway[match(list_class_2[[2]],rownames(etClass2Pathway)),]
etF9 = cor(etRadFeature9,t(etClass2Feature9),method= 'pearson')
etF9Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF9), length(etRadFeature9)))
etRadFeature9 = myNormalize(etRadFeature9)

etRadFeature4 = etRadFeature[,4]
etClass3Feature4 = etClass3Pathway[match(list_class_3[[1]],rownames(etClass3Pathway)),]
etF4 = cor(etRadFeature4,t(etClass3Feature4),method= 'pearson')
etF4Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF4), length(etRadFeature4)))
etRadFeature4 = myNormalize(etRadFeature4)

etRadFeature6 = etRadFeature[,6]
etClass3Feature6 = etClass3Pathway[match(list_class_3[[2]],rownames(etClass3Pathway)),]
etF6 = cor(etRadFeature6,t(etClass3Feature6),method= 'pearson')
etF6Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF6), length(etRadFeature6)))
etRadFeature6 = myNormalize(etRadFeature6)

etRadFeature7 = etRadFeature[,7]
etClass4Feature7 = etClass4Pathway[match(list_class_4[[1]],rownames(etClass4Pathway)),]
etF7 = cor(etRadFeature7,t(etClass4Feature7),method= 'pearson')
etF7Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF7), length(etRadFeature7)))
etRadFeature7 = myNormalize(etRadFeature7)

etRadFeature8 = etRadFeature[,8]
etClass4Feature8 = etClass4Pathway[match(list_class_4[[2]],rownames(etClass4Pathway)),]
etF8 = cor(etRadFeature8,t(etClass4Feature8),method= 'pearson')
etF8Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF8), length(etRadFeature8)))
etRadFeature8 = myNormalize(etRadFeature8)

etRadFeature10 = etRadFeature[,10]
etClass4Feature10 = etClass4Pathway[match(list_class_4[[3]],rownames(etClass4Pathway)),]
etF10 = cor(etRadFeature10,t(etClass4Feature10),method= 'pearson')
etF10Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF10), length(etRadFeature10)))
etRadFeature10 = myNormalize(etRadFeature10)

etRadFeature11 = etRadFeature[,11]
etClass4Feature11 = etClass4Pathway[match(list_class_4[[4]],rownames(etClass4Pathway)),]
etF11 = cor(etRadFeature11,t(etClass4Feature11),method= 'pearson')
etF11Pvalue = as.data.frame(corPvalueStudent(as.matrix(etF10), length(etRadFeature11)))
etRadFeature11 = myNormalize(etRadFeature11)


etClassFeature = rbind(etRadFeature1,etRadFeature2,etRadFeature3,etRadFeature4,etRadFeature5,
                        etRadFeature6,etRadFeature7,etRadFeature8,etRadFeature9,
                        etRadFeature10,etRadFeature11,etRadFeature12,etRadFeature13,
                      etClass2Feature1,etClass1Feature2,etClass1Feature3,etClass3Feature4,etClass1Feature5,
                      etClass3Feature6,etClass4Feature7,
                      etClass4Feature8,etClass2Feature9,
                      etClass4Feature10,etClass4Feature11,etClass1Feature12,etClass1Feature13)

etClinical = data.frame(etSig = etSig, etOs = etOs, etStatus = etStatus)
etPathway = data.frame(etClinical, t(etClassFeature))
etPathway = etPathway[order(etPathway[,1],decreasing=T),]
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
heatmap2PathwayGSVA(rgtPathway,etPathway,resultDir)




# ==================== heatmap keygene ====================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
resultDir =  paste(getwd(),"/output/radiogenomics/",sep = "")
setwd(resultDir)
rgtRedGene = keyGene(rgtData, rgtSig, 'red', rgtModuleColor, rgtMM)
rgtBlueGene = keyGene(rgtData, rgtSig, 'blue', rgtModuleColor, rgtMM)

rgtRedKeyGene = rgtRed[,na.omit(match(rgtRedGene[,1],colnames(rgtRed)))]
rgtBlueKeyGene = rgtBlue[,na.omit(match(rgtBlueGene[,1],colnames(rgtBlue)))] 
rgtKeyGene = cbind(rgtRedKeyGene, rgtBlueKeyGene)
rgtKeyGeneMean = as.numeric(apply(rgtKeyGene,1,mean))
drawKM.Rad_2(status = rgtStatus, time = rgtOs, sig = rgtKeyGeneMean, cutoff=-0.5, "rgt")

etRedKeyGene = etRed[,na.omit(match(rgtRedGene[,1],colnames(etRed)))]
etBlueKeyGene = etBlue[,na.omit(match(rgtBlueGene[,1],colnames(etBlue)))] 
etKeyGene = cbind(etRedKeyGene, etBlueKeyGene)
etKeyGeneMean = as.numeric(apply(etKeyGene,1,mean))
drawKM.Rad_2(status = etStatus, time = etOs, sig = etSig, cutoff=-0.5, "et")

tcgaRedKeyGene = tcgaRed[,na.omit(match(rgtRedGene[,1],colnames(tcgaRed)))]
tcgaBlueKeyGene = tcgaBlue[,na.omit(match(rgtBlueGene[,1],colnames(tcgaBlue)))] 
tcgaKeyGene = cbind(tcgaRedKeyGene, tcgaBlueKeyGene)
tcgaKeyGeneMean = as.numeric(apply(tcgaKeyGene,1,mean))
drawKM.Rad_2(status = tcgaStatus, time = tcgaOs, sig = tcgaSig, cutoff=-0.5, "tcga")

rgtClinical = data.frame(rgtSig = rgtSig, rgtOs = rgtOs, rgtStatus = rgtStatus)
rgtKey = cbind(rgtClinical,rgtKeyGeneMean,rgtKeyGene)
etClinical = data.frame(etSig = etSig, etOs = etOs, etStatus = etStatus)
etKey = cbind(etClinical,etKeyGeneMean,etKeyGene)
heatmapKeyGene(radioGenTrain_heatmap3_key_gene = rgtKey,externalTest_heatmap3_key_gene = etKey,resultDir = resultDir)

