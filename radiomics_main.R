# Radiomics analysis for survival prediction
# Date: 2021-04-19
# Authors: Qiuchang Sun and Zhi-Cheng Li


# ================================================import library=====================================================
# set path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./radiomics_function.R")
import_library()
# ================================================ load feature data and clincal data =====================================================
# read train data and test data
dataTrain = read.csv("./input/radiomics/train.csv")
dataTest = read.csv("./input/radiomics/test.csv")
dataAll = rbind(dataTrain,dataTest)
rowsAll = nrow(dataAll)
colsAll = ncol(dataAll)
rowTrain = nrow(dataTrain)
rowTest = nrow(dataTest)
# clinical data : Status and os
statusTrain = dataTrain[1:rowTrain,2]
osTrain = dataTrain[1:rowTrain,3]
statusTest = dataTest[1:rowTest,2]
osTest = dataTest[1:rowTest,3]

# path for figure
resultDir = paste(getwd(),"/output/radiomics",sep = "")
if(!(file.exists(resultDir))){
  dir.create(resultDir)}
# diffExam : calculated  P value for the difference in patient characteristics between RDC and RVC using Wilcoxon test
# output: P value (OS, sex, age, kps, IDH1, Radiation, Chemo, SurgResection)
diffExam(dataTrain=dataTrain, rowTrain=rowTrain, dataTest=dataTest, rowTest=rowTest)
# clinical 
clinicalAll = dataAll[1:rowsAll,6:12]
# imaging features
feeatureAll = dataAll[1:rowsAll,13:ncol(dataAll)]
numFeatures = ncol(feeatureAll)
Status = as.double(dataAll$Status)
os = as.double(dataAll$OS)

# =============================== feature selection ,Radiomics signature construction and validation =====================================================
# preprocessing : univariate feature selection , zscore feature Normalization
# input: feeatureAll, numFeatures, rowsAll, rowTrain, rowTest, Status, os
# output: featureSelectedTrain,featureSelectedTest,numFeatureSelected,result1:univariate prognostic result
dataList = preProcessing(feeatureAll, numFeatures, rowsAll, rowTrain, rowTest, Status, os)
featureSelectedTrain = dataList[[1]]
featureSelectedTest = dataList[[2]]
numFeatureSelected = dataList[[3]]
result1 = dataList[[4]]
# predictmodel.Rad.OS : a radiomics signature was built on the training subset by linearly combining remaining features using LASSO
# input: featureSelectedTrain, statusTrain, osTrain, featureSelectedTest, statusTest, osTest
# output: rad_CI_train,sigTrain,rad_CI_test,sigTest,active.Coefficients1,active.Index1
modelRadPredict = predictModel.Rad.OS(featureSelectedTrain, statusTrain, osTrain, featureSelectedTest, statusTest, osTest)
rad_CI_train = modelRadPredict[[1]]
rad_CI_test = modelRadPredict[[3]]
sigTrain = modelRadPredict[[2]]
sigTest = modelRadPredict[[4]]
active.Coefficients1 = modelRadPredict[[5]]
active.Index1 = modelRadPredict[[6]]


# =============================== km curve, calibration curve, decision curve =====================================================
# Rad.Xtilefile.make : make Xtilefile 
# input: statusTrain, osTrain, statusTest, osTest, sigTrain, sigTest,resultDir
# output: a txt 
rad.Xtilefile.Make(statusTrain, osTrain, statusTest, osTest, sigTrain, sigTest,resultDir)
# drawKM.Rad : draw km curve
# statusTrain, osTrain, statusTest, osTest, sigTrain, sigTest,cutoff:X-tile choose -4 here. The users should replace it with their cutoff (maybe calculated using X-tile or other tools)
# km curve
drawKM.Rad(statusTrain, osTrain, statusTest, osTest, sigTrain, sigTest,-4,resultDir) 
# calibration curve , decision curve
train = data.frame(cbind(osTrain, statusTrain, sigTrain))
name = '/cal_train_OS.tiff'
drawCal(train, resultDir,name, survival="OS", 44)
name = '/dc_train_OS.tiff'
drawDCA(train, resultDir,name, survival="OS")

test = data.frame(cbind(osTest, statusTest, sigTest))
name = '/cal_test_OS.tiff'
drawCal(test, resultDir, name, survival="OS", 23)
name = '/dc_test_OS.tiff'
drawDCA(test, resultDir, name, survival="OS")



# =================================================== Forest plot ===========================================================================
tiff(file = paste(resultDir,'/figure2_A.tiff',sep=''), res = 600, width = 8000, height = 4000, compression = "lzw")
feature_CI_Plot = read.csv("./input/radiomics/featureci.csv")
a = 1 - feature_CI_Plot[1:9,2] 
b = 1 - feature_CI_Plot[1:9,3]
c = 1 - feature_CI_Plot[1:9,4]
feature_CI_Plot[1:9,2] = a
feature_CI_Plot[1:9,3] = c
feature_CI_Plot[1:9,4] = b
tabletext<-cbind(
  c("", "f1", "f2", "f3", "f4", "f5", "f6", "f7", 
    "f8", "f9", "f10", "f11", 
    "f12", "f13"),
  c("Image Features", "GLCM.Inverse Difference Moment", "GLSZM.Large Zone Emphasis ","GLSZM.Large Zone Low Gray Level Emphasis","GLCM.Sum Average",
    "GLCM.Inverse Difference Moment","GLSZM.Zone Size Variance",
    "GLCM.Sum Entropy","GLSZM.Large Zone Emphasis ", "GLCM.Cluster Shade", "GLSZM.Small Zone Low Gray Level Emphasis", "GLSZM.Gray Level Variance","GLCM.Autocorrelation","GLCM.Sum Average"),
  c("Modality", "T1w", "T1c", "Flair","T2w", "Flair", "Flair", "T1c","T1w", "Flair", "T2w", "T1c",
    "T1c","T1c"),
  c("Subregion", "Core", "Core", "Whole", "Whole","Whole","Whole","Core", "Core", "Edema", "Core",
    "Core", "Core","Whole"))

forestplot(tabletext, hrzl_lines = list("2" = gpar(lwd=1)), 
           mean = c(NA,feature_CI_Plot[,2]),
           lower = c(NA,feature_CI_Plot[,3]),
           upper = c(NA,feature_CI_Plot[,4]),
           align = "c",
           boxsize = .3,
           xlab = "Concordance Index(95% CI)",
           zero = 1,
           is.summary=c(TRUE,rep(FALSE,13)),
           clip =c(0.5, 0.8),
           graphwidth = unit(8, "cm"),
           lineheight = unit(1, "cm"), 
           col=fpColors(box=c("blue")),
           txt_gp = fpTxtGp(label = list(gpar(fontface = "bold.italic", cex = 1),
                                         gpar(fontface = "bold", cex = 1),
                                         gpar(fontface = "bold", cex = 1),
                                         gpar(fontface = "bold", cex = 1)),
                            summary = gpar(fontface = "bold", cex=1.2),
                            ticks = gpar(fontface = "bold", cex=1.2),
                            xlab  = gpar(fontface = "bold", cex = 1.2)),
           new_page = TRUE,
           lwd.zero = FALSE)
dev.off()

