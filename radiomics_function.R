# Radiomics analysis for survival prediction
# Date: 2021-04-19
# Authors: Qiuchang Sun and Zhi-Cheng Li


# ====================function1 import library ====================
import_library = function() {
  wants <- c("survival","prodlim","glmnet","Matrix","foreach","gplots","MASS","stringi","simPH","survcomp","km.ci","readxl","dlstats","rms",
           "Hmisc","lattice","Formula","ggplot2","SparseM","GGally","rmda","tidyverse","tibble","survminer","survIDINRI","forestplot")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  has <- wants %in% rownames(installed.packages())
  if(any(!has)) BiocManager::install(wants[!has])
  sapply(wants, require, character.only = TRUE)
  }

# =============================================== function1 diff_Exam ===================================================
# calculated  P value for the difference in patient characteristics between RDC and RVC using Wilcoxon test
# ========input========
# dataTrain, rowTrain, dataTest, rowTest
# ========output========
# P value (OS, sex, age, kps, IDH1, Radiation, Chemo, SurgResection)

diffExam <- function(dataTrain, rowTrain, dataTest, rowTest) {
  trainStatus <- matrix(as.double(dataTrain[1:rowTrain,2]))
  testStatus <- matrix(as.double(dataTest[1:66,2]))

  trainOS <- matrix(as.double(dataTrain[1:rowTrain,3]))
  testOS <- matrix(as.double(dataTest[1:66,3]))
  
  surv1 = data.frame(status = c(trainStatus,testStatus),os = c(trainOS,testOS),set = c(rep("train",132),rep("test",66)))
  
  sexTrain <- matrix(as.double(dataTrain[1:rowTrain,6]))
  sexTest  <- matrix(as.double(dataTest[1:rowTest,6]))
  
  ageTrain <- matrix(as.double(dataTrain[1:rowTrain,7]))
  ageTest  <- matrix(as.double(dataTest[1:rowTest,7]))
  
  kpsTrain <- matrix(as.double(dataTrain[1:rowTrain,8]))
  kpsTest  <- matrix(as.double(dataTest[1:rowTest,8]))
  
  IDH1Train <- matrix(as.double(dataTrain[1:rowTrain,9]))
  IDH1Test  <- matrix(as.double(dataTest[1:rowTest,9]))
  
  radiationTrain <- matrix(as.double(dataTrain[1:rowTrain,10]))
  radiationTest  <- matrix(as.double(dataTest[1:rowTest,10]))
  
  chemoTrain <- matrix(as.double(dataTrain[1:rowTrain,11]))
  chemoTest  <- matrix(as.double(dataTest[1:rowTest,11]))
  
  surgresectionTrain <- matrix(as.double(dataTrain[1:rowTrain,12]))
  surgresectionTest  <- matrix(as.double(dataTest[1:rowTest,12]))
  
  htestOS = survdiff(Surv(os, status) ~ set, data = surv1)
  htestSex  = wilcox.test(sexTrain,sexTest)
  htestAge  = wilcox.test(ageTrain,ageTest)
  htestKps  = wilcox.test(kpsTrain,kpsTest)
  htestIDH1  = wilcox.test(IDH1Train,IDH1Test)
  htestRadiation = wilcox.test(radiationTrain,radiationTest)
  htestChemo = wilcox.test(chemoTrain,chemoTest)
  htestSurgResection = wilcox.test(surgresectionTrain,surgresectionTest)
  
  print(paste("os p value:",htestOS$chisq))
  print(paste("sex p value:",htestSex$p.value))
  print(paste("age p value:",htestAge$p.value))
  print(paste("kps p value:",htestKps$p.value))
  print(paste("IDH1 p value:",htestIDH1$p.value))
  print(paste("radiation p value:",htestRadiation$p.value))
  print(paste("chemo p value:",htestChemo$p.value))
  print(paste("surgResection p value:",htestSurgResection$p.value))
}


# =============================================== function2 feature selection ===================================================
# univariate feature selection , zscore feature Normalization
# ========input========
# featureAll,numFeatures,rowsAll,rowTrain, rowTest,status,OS
# ========output========
# featureSelectedTrain,featureSelectedTest,numFeatureSelected,result1

preProcessing <- function(featureAll,numFeatures,rowsAll,rowTrain, rowTest,status,OS) {
  dataFeature=data.matrix(featureAll)
  dataEvent=data.matrix(status)
  dataTime=data.matrix(OS)
  nFeature=ncol(dataFeature)
  
  result1 <- matrix(nrow=numFeatures,ncol=1)
  
  for(i in 1:nFeature)
  {ciList=concordance.index(x=dataFeature[1:rowsAll,i],surv.time=as.numeric(dataTime),surv.event=dataEvent,method="noether")
  result1[i,1]=ciList[['c.index']]
  }
  rm(ciList)
  # univariate prognostic analysis
  rownames(result1)=t(colnames(dataFeature))
  nFeaSelected = 1
  for(i in 1:nFeature)
  {
    if (((result1[i,1])-0.5>=0.05) | (0.5-(result1[i,1])>=0.05))
    {
      nFeaSelected = nFeaSelected+1
    }
    else
      result1[i,1] = 0
  }
  featureSelected = matrix(NA,rowsAll,(nFeaSelected-1))
  j =1
  for (i in 1:(nFeature))
  {
    if ((result1[i,1])!=0)
    {
      featureSelected[1:rowsAll,j] = dataFeature[1:rowsAll,i]
      j = j+1
    }
  }
  numFeatureSelected = ncol(featureSelected)
  featureSelectedTrain = featureSelected[1:rowTrain,1:numFeatureSelected]
  featureSelectedTest = featureSelected[(rowTrain+1):(rowTest+rowTrain),1:numFeatureSelected]
  
  # zscore feature Normalization
  featureSelectedTrainNorm = featureSelectedTrain
  for (i in 1:numFeatureSelected)
  {
    meanTrain=mean(featureSelectedTrain[1:rowTrain,i])
    sdTrain=sd(featureSelectedTrain[1:rowTrain,i])
    featureSelectedTrainNorm[i]=(featureSelectedTrain[1:rowTrain,i]-meanTrain)/sdTrain
  }
  featureSelectedTrain[1:rowTrain,1:numFeatureSelected] = featureSelectedTrainNorm[1:rowTrain,1:numFeatureSelected]
  
  featureSelectedTestNorm = featureSelectedTest
  for (i in 1:numFeatureSelected)
  {
    meanTest=mean(featureSelectedTest[1:rowTest,i])
    sdTest=sd(featureSelectedTest[1:rowTest,i])
    featureSelectedTestNorm[i]=(featureSelectedTest[1:rowTest,i]-meanTest)/sdTest
  }
  featureSelectedTest[1:rowTest,1:numFeatureSelected] = featureSelectedTestNorm[1:rowTest,1:numFeatureSelected]
  
  return(list(featureSelectedTrain,featureSelectedTest,numFeatureSelected,result1))
}


# =============================================== function3 construct predict model ===================================================
# a radiomics signature was built on the training subset by linearly combining remaining features using LASSO
# ========input========
# featureSelectedTrain, statusTrain, osTrain, featureSelectedTest, statusTest, osTest
# ========output========
# cindexTrain, sigTrain, cindexTest, sigTest, active.Coefficients, active.Index

predictModel.Rad.OS <- function(featureSelectedTrain, statusTrain, osTrain, featureSelectedTest, statusTest, osTest) {
  rad_feature_Train <- as.matrix(featureSelectedTrain)
  trainStatus <- as.double(statusTrain)
  trainTime <- as.double(osTrain)
  sTrain <- Surv(time = trainTime,event = trainStatus,type='right')
  
  rad_feature_Test <- as.matrix(featureSelectedTest)
  testStatus <- as.double(statusTest)
  testTime <- as.double(osTest)
  sTest <- Surv(time = testTime,event = testStatus,type='right')
#using 10-fold CV to select lambda:
  res1<-cv.glmnet(rad_feature_Train, sTrain, family="cox", alpha = 1, nfolds=10, nlambda=100)
  plot(res1,xlab=expression(log(lambda)))
  log(res1$lambda)
  
  res2<-glmnet(rad_feature_Train, sTrain, family="cox", alpha = 1,lambda=res1$lambda.min)
  NonZero=res2$df
  if (!NonZero)
    next
  Coefficients <- coef(res2, s = res2$lambda.min)
  active.Index <- which(Coefficients != 0)
  active.Coefficients <- Coefficients[active.Index]
  nPatientTrain = nrow(rad_feature_Train)
  nFeature = ncol(rad_feature_Train)
  sigTrain <- matrix(nrow=nPatientTrain,ncol=1)
  sigTrain1 <- matrix(nrow=res2$df,ncol=1)
  for (i in 1:nPatientTrain)
  {
    for(j in 1:res2$df)
    {
      sigTrain1[j,1] = active.Coefficients[j] * rad_feature_Train[i,active.Index[j]]
      sigTrain[i,1] =  sum(sigTrain1)
    }
  }
  nPatientTest = nrow(rad_feature_Test)
  nFeature = ncol(rad_feature_Test)
  sigTest <- matrix(nrow=nPatientTest,ncol=1)
  sigTest1 <- matrix(nrow=res2$df,ncol=1)
  
  for (i in 1:nPatientTest)
  {
    for(j in 1:res2$df)
    {
      sigTest1[j,1] = active.Coefficients[j] * rad_feature_Test[i,active.Index[j]]
      sigTest[i,1] =  sum(sigTest1)
    }
  }
  
  sample.Data1 <- data.frame(feature1=sigTrain,os1 = trainTime,death1 = trainStatus)
  sum.Surv1 <- summary(coxph(Surv(os1, death1) ~ feature1,data = sample.Data1))
  cindexTrain <-sum.Surv1$concordance[1]
  
  sample.Data2 <- data.frame(feature2=sigTest,os2 = testTime,death2 = testStatus)
  sum.Surv2 <- summary(coxph(Surv(os2, death2) ~ feature2,data = sample.Data2))
  cindexTest <-sum.Surv2$concordance[1]
  
  return(list(cindexTrain, sigTrain, cindexTest, sigTest, active.Coefficients, active.Index))
}


# =============================================== function4 make Xtilefile ===================================================
# make Xtilefile
# ========input========
# statusTrain, timeTrain, statusTest, timeTest, sigTrain, sigTest, resultDir
# ========output========
# Xtilefile

rad.Xtilefile.Make <- function(statusTrain, timeTrain, statusTest, timeTest, sigTrain, sigTest, resultDir) {
    trainStatus <- as.double(statusTrain)
    trainTime <- as.double(timeTrain)
    testStatus <- as.double(statusTest)
    testTime <- as.double(timeTest)

    time<-vector(mode="numeric",length = 198)
    time[1:132]<-trainTime
    time[133:198]<-testTime
    status<-vector(mode="character",length = 198)
    for(i in 1:132){
      if(trainStatus[i]==1){
        status[i]<-"uncensored"
      }
      else{
        status[i]<-"censored"
      }
    }
    for(i in 1:66){
      if(testStatus[i]==1){
        status[132+i]<-"uncensored"
      }
      else{
        status[132+i]<-"censored"
      }
    }
    signature<-vector(mode="numeric",length = 198)
    signature[1:132]<-sigTrain
    signature[133:198]<-sigTest

    valid<-vector(mode="character",length = 198)
    valid[1:132]<-"Training"
    valid[133:198]<-"Validation"

    x<-data.frame("Time"=time,"Status"=status,"Signature"=signature,"valid"=valid)
    write.table(x,file=paste(resultDir,"/Train.txt",sep = ""),row.names=F,quote=F,sep="\t")
  }


# =============================================== function5 draw km ===================================================
# draw km curve
# ========input========
# statusTrain, timeTrain, statusTest, timeTest, sigTrain, sigTest,cutoff,resultDir
# ========output========
# km curve
# ========figure========
# 1.km curve

# ========to beautify the km curve========
km.coxph.plot_1 <-
  function(formula.s, data.s, weight.s, x.label, y.label, main.title, sub.title, leg.text, leg.pos="bottomright", leg.bty="o", leg.inset=0.05, o.text, v.line, h.line, .col=1:4, .lty=1, .lwd=1, show.n.risk=FALSE, n.risk.step, n.risk.cex=0.85, verbose=TRUE, ...) {
    
    if (missing(sub.title)) { sub.title <- NULL }
    if (missing(leg.text)) { leg.text <- NULL }
    if (missing(weight.s)) { weight.s <- array(1, dim=nrow(data.s), dimnames=list(rownames(data.s))) }
    ## weights should be > 0
    data.s <- data.s[!is.na(weight.s) & weight.s > 0, , drop=FALSE]
    weight.s <- weight.s[!is.na(weight.s) & weight.s > 0]
    pos <- 1
    envir = as.environment(pos)
    assign("weight.s", weight.s, envir = envir)
    weighted <- length(sort(unique(weight.s))) > 1
    
    ng <- length(leg.text)
    old.mar <- par("mar")
    on.exit( par( mar = old.mar ) )
    .xaxt="s"
    .xlab=x.label
    if (show.n.risk) {
      par(mar = old.mar + c(ng,8,3,0))
      .xaxt="n"
      .xlab = ""
    }
    
    plot(survfit(formula.s, data=data.s, weights=weight.s), xaxt=.xaxt, col=.col, lty=.lty, lwd=.lwd, xlab=.xlab, ylab=y.label, ... )
    title(main.title)
    
    if (!missing(v.line) && !is.null(v.line)) { abline(v=v.line, lty=3, col="purple") }
    if (!missing(h.line) && !is.null(h.line)) { abline(h=h.line, lty=3, col="purple") }
    
    if (!is.null(leg.text)) { legend(x=leg.pos, xjust=0, yjust=1, legend=leg.text, text.font =2,col=.col, lty=.lty, lwd=.lwd, cex=0.8, bg="white", inset=leg.inset, bty=leg.bty) }
    if (!is.null(sub.title)) { mtext(sub.title, line=-4, outer=TRUE) }
    if (missing(o.text)) {
      sdf <- summary(survival::coxph(formula.s, data=data.s, weights=weight.s))
      if(verbose) { print(sdf) }
      p.val <- sdf$sctest["pvalue"]
      o.text <- sprintf("Logrank P = %.1E", p.val)
    }
    if (is.null(o.text)) { o.text <- FALSE }
    text(0,0, o.text, cex=0.85, pos=4,font = 2)
    
    if (show.n.risk) {
      usr.xy <- par( "usr" )
      nrisk <- no.at.risk(formula.s=formula.s, data.s=data.s, sub.s="all", t.step=n.risk.step, t.end=floor(usr.xy[2]) )
      at.loc <- seq(0, usr.xy[2], n.risk.step)
      axis(1, at=at.loc, font = 2)
      mtext(x.label, side=1, line=2, font=2)
      mtext("No. At Risk", side=1, line=3, at=-0.5*n.risk.step, adj=1, cex=n.risk.cex, font=2)
      #nrsk.lbs <- sapply( strsplit(levels(nrisk[,1]),"="), FUN=function(x) x[2] )
      #if( any(is.na(nrsk.lbs)) ) nrsk.lbs <- leg.text
      for( i in 1:nrow(nrisk) ) {
        mtext(leg.text[i], side=1, line=3+i, at=-0.5*n.risk.step, adj=1, cex=n.risk.cex, font=2)
        mtext(nrisk[i,-1], side=1, at=at.loc, line=3+i, adj=0.5, cex=n.risk.cex, font=2)
      }
    }
    
    if( exists("weight.s", envir=.GlobalEnv) ) remove("weight.s", envir=.GlobalEnv)
  }
# ========draw km ========
drawKM.Rad <- function(statusTrain, timeTrain, statusTest, timeTest, sigTrain, sigTest,cutoff,resultDir) {
  
  trainStatus <- as.double(statusTrain)
  trainTime <- as.double(timeTrain)
  testStatus <- as.double(statusTest)
  testTime <- as.double(timeTest)
  
  stratRadTrain = sigTrain
  stratRadTrain[which(sigTrain >= cutoff)]=1
  stratRadTrain[which(sigTrain < cutoff)]=0

  dd1 <- data.frame("surv.time"=trainTime, "surv.event"=trainStatus, "strat"=stratRadTrain)
  tiff(file = paste(resultDir,'/rad_train_km.tiff',sep=''), res = 600, width = 6000, height = 4000, compression = "lzw")
  par(mar=c(0,0,0,0))
  km.coxph.plot_1(formula.s=Surv(trainTime, trainStatus) ~ stratRadTrain, data.s=dd1,
                x.label="Time (days)", y.label="Probability of survival",mark.time = T,pch="+",cex=2,.lwd = 4,
                main.title="",.col=c("darkgreen", "red3"), n.risk.step = 200,
                .lty=c(1,1), show.n.risk=T, n.risk.cex=0.85, verbose=T,
                leg.text=paste(c("Low risk", "High risk"), "   ", sep=""),leg.bty = "n",
                leg.pos="topright",font.axis=2,cex.axis=1.2,font.lab=2,cex.lab=1.2,font.legend=2)
  dev.off()
  
  #HR1=hazard.ratio(x=SigTrain,surv.time = trainTime,surv.event = trainStatus,strat = Strat1,method.test = "logrank")
  coxphSigTrain=coxph(Surv(trainTime,trainStatus)~sigTrain)
  coxphStratRadTrain=coxph(Surv(trainTime,trainStatus)~stratRadTrain)
  print(summary(coxphSigTrain))
  print(summary(coxphStratRadTrain))
  
  
  stratRadTest = sigTest
  stratRadTest[which(sigTest >= cutoff)]=1
  stratRadTest[which(sigTest < cutoff)]=0
  
  dd2 <- data.frame("surv.time"=testTime, "surv.event"=testStatus, "strat"=stratRadTest)
  tiff(file = paste(resultDir,'/rad_test_km.tiff',sep=''), res = 600, width = 6000, height = 4000, compression = "lzw")
  par(mar=c(0,0,0,0))
  km.coxph.plot_1(formula.s=Surv(testTime, testStatus) ~ stratRadTest, data.s=dd2,
                x.label="Time (days)", y.label="Probability of survival",
                main.title="",.col=c("darkgreen", "red3"), n.risk.step = 200,mark.time = T,pch="+",cex=2,.lwd = 4,
                .lty=c(1,1), show.n.risk=T, n.risk.cex=0.85, verbose=T,
                leg.text=paste(c("Low risk", "High risk"), "   ", sep=""),leg.bty = "n",
                leg.pos="topright",font.axis=2,cex.axis=1.2,font.lab=2,cex.lab=1.2,font.legend=2)
  dev.off()
  #HR2=hazard.ratio(x=SigTest,surv.time = testTime,surv.event = testStatus,strat = Strat2,method.test = "logrank")
  coxphSigTest=coxph(Surv(testTime,testStatus)~sigTest)
  coxphStratRadTest=coxph(Surv(testTime,testStatus)~stratRadTest)
  print(summary(coxphSigTest))
  print(summary(coxphStratRadTest))
}


# =============================================== function6 draw calibration ===================================================
# draw calibration curve
# ========input========
# data, path, name,survival="OS", M
# ========output========
# calibration curve
# ========figure========
# 1.calibration curve

# ========beautify the calibration curve========
drawCal <- function(data, path, name,survival="OS", M)
{
  radSignature = data[,3]
  os = data[,1]
  status = data[,2]
  M = M
  path = path
  name = name
  radClinicAll = data.frame(cbind(os, status, radSignature))
  y <- function(x) x
 
  if (survival == "OS"){
    # clinic
    f2.1 <- psm(Surv(os, status) ~ radSignature,
                data = radClinicAll, x=T, y=T, dist='lognormal')
    cal1 <- calibrate(f2.1, cmethod='KM', method="boot", u=180, m=M, B=200)
    cal2 <- calibrate(f2.1, cmethod='KM', method="boot", u=360, m=M, B=200)
    cal3 <- calibrate(f2.1, cmethod='KM', method="boot", u=540, m=M, B=200)
    
    tiff(file = paste(path,name, sep=''), res = 600,width = 8000, height = 6000, compression = "lzw")
    par(oma = c(1, 2, 1, 1), mar = c(5, 5, 3, 0) + 0.1)
    myCalPlot(cal1,
              errbar.col=c(rgb(192,98,83,maxColorValue=255)),
              xlab="Predicted Probability of OS",
              ylab="Observed fraction OS proportion",
              font.lab = 2,
              cex.lab = 2,
              xlim=c(0.0,1.0),ylim=c(0.0,1.0),
              col=c(rgb(192,98,83,maxColorValue=255)))
    myCalPlot(cal2,
              errbar.col=c(rgb(98,192,83,maxColorValue=255)),
              col=c(rgb(98,192,83,maxColorValue=255)),
              add = TRUE)
    myCalPlot(cal3,
              errbar.col=c(rgb(98,83,192,maxColorValue=255)),
              col=c(rgb(98,83,192,maxColorValue=255)),
              add = TRUE)
    plot(y, lwd=4, add=TRUE, col="gray")
    legend(0,1,legend=c("ideal",
                        "6-month OS",
                        "12-month OS",
                        "18-month OS"),
           col=c("Gray",rgb(192,98,83,maxColorValue=255),
                 rgb(98,192,83,maxColorValue=255),
                 rgb(98,83,192,maxColorValue=255)),
           text.font=2,
           bty = "o", ncol = 1, cex = 1.5, lwd = 4)
    dev.off()
  }}


# =============================================== function7 draw decision curve ===================================================
# draw decision curve
# ========input========
# data, path, name,survival="OS"
# ========output========
# decision curve
# ========figure========
# 1.decision curve

# ========beautify the decision curve========
drawDCA <- function(data, path, name, survival="OS")
{
  radSignature = data[,3]
  OS = data[,1]
  status = data[,2]
  path = path
  name = name
  radClinicAll = data.frame(cbind(OS, status, radSignature))
  as_tibble(radClinicAll)
  
  if (survival=="OS"){
    clinic_model = decision_curve(status ~ radSignature,
                                   data = radClinicAll,
                                   thresholds = seq(0, 1, by = .01),
                                   bootstraps = 100,
                                   family = binomial(link ='logit'))
    par(c(5,4,4,2)+0.1)
    tiff(file = paste(path,name, sep=''), res = 600,width = 8000, height = 6000, compression = "lzw")
    # postscript(file = paste(path, '/dca_radtest_OS.eps',sep=''))
    plot_decision_curve_new(clinic_model,
                            curve.names = "Radiomics signature",
                            xlab="High Risk Threshold",
                            ylab="Net Benefit",
                            font.lab = 2,
                            cex.lab = 1.5,
                            col = "blue",
                            xlim = c(0.6, 1),
                            ylim = c(-0.1, 0.8),
                            confidence.intervals = F, #remove confidence intervals
                            cost.benefit.axis = FALSE, #remove cost benefit axis
                            legend.position = "topright", #remove the legend
                            legend(1,legend=c("Radionomogram",
                                              "All",
                                              "None"),
                                   col=c("blue","grey66","black"),
                                   text.font=4,
                                   ncol = 1, cex = 10, lwd = 4),
                            standardize = F)
    dev.off()
  }
}
myCalPlot <- function (x, xlab="", ylab="", subtitles = FALSE, conf.int = TRUE, cex.subtitles = 0.7,
                       riskdist = TRUE, add = FALSE, scat1d.opts = list(nhistSpike = 200), 
                       par.corrected = NULL, width = 4, ...) 
{
  at <- attributes(x)
  u <- at$u
  units <- at$units
  if (length(par.corrected) && !is.list(par.corrected)) 
    stop("par.corrected must be a list")
  z <- list(col = "blue", lty = 2, lwd = 10, pch = 10)
  if (!length(par.corrected)) 
    par.corrected <- z
  else for (n in setdiff(names(z), names(par.corrected))) par.corrected[[n]] <- z[[n]]
  predicted <- at$predicted
  if ("KM" %in% colnames(x)) {
    type <- "stratified"
    pred <- x[, "mean.predicted"]
    cal <- x[, "KM"]
    cal.corrected <- x[, "KM.corrected"]
    se <- x[, "std.err"]
  }
  else {
    type <- "smooth"
    pred <- x[, "pred"]
    cal <- x[, "calibrated"]
    cal.corrected <- x[, "calibrated.corrected"]
    se <- NULL
  }
  un <- if (u == 1) 
    paste(units, "s", sep = "")
  else units
  
  if (length(se) && conf.int) {
    ciupper <- function(surv, d) ifelse(surv == 0, 0, pmin(1, 
                                                           surv * exp(d)))
    cilower <- function(surv, d) ifelse(surv == 0, 0, surv * 
                                          exp(-d))
    # ÐÞ¸Älwd=width
    my_errbar(pred, cal, cilower(cal, 1.959964 * se), ciupper(cal, 1.959964 * se), lwd = width,
              xlab= xlab, ylab = ylab, font.lab = 1, cex.lab = 1.5, type = "b",
              bty = "n", xaxt = "n", yaxt ="n", 
              add = add, ...)
    axis(1, seq(0, 1, 0.1), seq(0, 1, 0.1), font.axis = 2, cex.axis = 1.8, lwd = 2)
    axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1), font.axis = 2, cex.axis = 1.8, lwd = 2)
    # axis(3)
    # axis(4)
    box(lwd = 2)
  }
  else if (add) 
    lines(pred, cal, type = if (type == "smooth") 
      "l"
      else "b")
  else plot(pred, cal, xlab = xlab, ylab = ylab, type = if (type == 
                                                            "smooth") 
    "l"
    else "b", ...)
  err <- NULL
  # ÐÞ¸Äpch
  if (type == "stratified")
    points(pred, cal.corrected, pch = 4, lwd = 4, cex = 2,
           col = "blue")
  else lines(pred, cal.corrected, col = par.corrected$col, 
             lty = par.corrected$lty, lwd = par.corrected$lwd)
  invisible()
}

my_errbar <- function (x, y, yplus, yminus, cap = 0.015, main = NULL, sub = NULL, 
                       xlab = as.character(substitute(x)), ylab = if (is.factor(x) || 
                                                                      is.character(x)) "" else as.character(substitute(y)), 
                       add = FALSE, lty = 2, type = "b", ylim = NULL, lwd = 2, pch = 16, width = 4,
                       errbar.col = par("fg"), Type = rep(1, length(y)), ...) 
{
  if (is.null(ylim)) 
    ylim <- range(y[Type == 1], yplus[Type == 1], yminus[Type == 
                                                           1], na.rm = TRUE)
  if (is.factor(x) || is.character(x)) {
    x <- as.character(x)
    n <- length(x)
    t1 <- Type == 1
    t2 <- Type == 2
    n1 <- sum(t1)
    n2 <- sum(t2)
    omai <- par("mai")
    mai <- omai
    mai[2] <- max(strwidth(x, "inches")) + 0.25
    par(mai = mai)
    on.exit(par(mai = omai))
    plot(NA, NA, xlab = ylab, ylab = "", xlim = ylim, ylim = c(1, 
                                                               n + 1), axes = FALSE, main = main, sub = sub, ...)
    axis(1)
    w <- if (any(t2)) 
      n1 + (1:n2) + 1
    else numeric(0)
    axis(2, at = c(seq.int(length.out = n1), w), labels = c(x[t1], 
                                                            x[t2]), las = 1, adj = 1)
    points(y[t1], seq.int(length.out = n1), pch = pch, type = type, 
           ...)
    segments(yplus[t1], seq.int(length.out = n1), yminus[t1], 
             seq.int(length.out = n1), lwd = lwd, lty = lty, col = errbar.col)
    if (any(Type == 2)) {
      abline(h = n1 + 1, lty = 2, ...)
      offset <- mean(y[t1]) - mean(y[t2])
      if (min(yminus[t2]) < 0 & max(yplus[t2]) > 0) 
        lines(c(0, 0) + offset, c(n1 + 1, par("usr")[4]), 
              lty = 2, ...)
      points(y[t2] + offset, w, pch = pch, type = type, 
             ...)
      segments(yminus[t2] + offset, w, yplus[t2] + offset, 
               w, lwd = lwd, lty = lty, col = errbar.col)
      at <- pretty(range(y[t2], yplus[t2], yminus[t2]))
      axis(side = 3, at = at + offset, labels = format(round(at, 
                                                             6)))
    }
    return(invisible())
  }
  
  # ÐÞ¸Älwd=width
  if (add) 
    points(x, y, pch = pch, type = type, lwd=width, ...)
  else plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, pch = pch, lwd=width,
            type = type, ...)
  
  xcoord <- par()$usr[1:2]
  smidge <- cap * (xcoord[2] - xcoord[1])/2
  segments(x, yminus, x, yplus, lty = lty, lwd = lwd, col = errbar.col)
  if (par()$xlog) {
    xstart <- x * 10^(-smidge)
    xend <- x * 10^(smidge)
  }
  else {
    xstart <- x - smidge
    xend <- x + smidge
  }
  segments(xstart, yminus, xend, yminus, lwd = lwd, lty = lty, 
           col = errbar.col)
  segments(xstart, yplus, xend, yplus, lwd = lwd, lty = lty, 
           col = errbar.col)
  return(invisible())
}

plot_decision_curve_new <- function(x, curve.names,
                                    cost.benefit.axis = FALSE,
                                    n.cost.benefits = 6,
                                    cost.benefits,
                                    standardize = FALSE,
                                    confidence.intervals,
                                    col,
                                    lty, lwd = 4,
                                    xlim, ylim,
                                    xlab, ylab,
                                    cost.benefit.xlab,
                                    legend.position = c("topright", "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "none"),
                                    ...){
  
  legend.position <- match.arg(legend.position)
  
  if(missing(curve.names)) curve.names  <- NA
  if(missing(confidence.intervals)) confidence.intervals <- NA
  
  prepData <- preparePlotData(x = x,
                              curve.names = curve.names,
                              confidence.intervals = confidence.intervals)
  
  predictors <- prepData$predictors
  dc.data <- prepData$dc.data
  confidence.intervals <- prepData$confidence.intervals
  rm(prepData)
  
  #set some defaults if needed
  if(missing(xlim)) xlim = range(dc.data$thresholds)
  
  if(missing(lty)) lty = rep(1, length(predictors) + 2)
  if(length(lty) ==1) lty = rep(lty, length(predictors) + 2)
  if(length(lty) == length(predictors)) lty = c(lty, 1, 1)
  
  if(missing(col)) col  = c(rainbow(length(predictors), v = .8), "grey66", "black")
  if(length(col) == length(predictors)) col <- c(col, "grey66", "black")
  
  if(missing(lwd)) lwd = 4
  if(length(lwd) ==4) lwd <- rep(lwd, length(predictors))
  if(length(lwd) == length(predictors)) lwd = c(lwd, 4, 4)
  
  if(missing(ylab)) ylab <- ifelse(standardize, "Standardized Net Benefit", "Net Benefit")
  policy = ifelse(class(x)=="decision_curve", x$policy, x[[1]]$policy)
  if(missing(xlab)) xlab <- ifelse(policy  == 'opt-in', "High Risk Threshold", "Low Risk Threshold")
  if(missing(cost.benefit.xlab)) cost.benefit.xlab <- "Cost:Benefit Ratio"
  if(missing(ylim)){
    
    if(standardize) ylim = c(-0.05, 1)
    else ylim = c(-0.05, 1.1*max(dc.data[["NB"]][is.finite(dc.data[["NB"]])]))
    
  }
  
  plot_generic(xx = dc.data,
               predictors = predictors,
               value = ifelse(standardize, "sNB", "NB"),
               plotNew = TRUE,
               standardize = standardize,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = lty, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = legend.position,
               policy = policy,
               ...)
  
}
preparePlotData   <- function(x, curve.names, confidence.intervals){
  #  if(missing(curve.names)) curve.nams <- NA
  #extract data from x,
  #there is only one DecisionCurve to plot
  if(class(x) == "decision_curve"){
    
    dc.data <- x$derived.data
    
    if(is.na(confidence.intervals[1])){
      confidence.intervals <- x$confidence.intervals
    }else{
      confidence.intervals <- ifelse(confidence.intervals, 1, "none")
    }
    
    if(is.na(curve.names[1])){
      predictors <- unique(dc.data$model)
      predictors <- predictors[!is.element(predictors, c("None", "All"))]
    }else{
      dc.data$model[!is.element(dc.data$model, c("None", "All"))] <- curve.names[1]
      predictors <- curve.names[1]
    }
    
  }else if(class(x)=="list"){
    xx <- NULL
    #check to make sure each element of the list is a decision_curve object,
    # or else we get funky results.
    if(!all(sapply(x, FUN = function(xx) class(xx) == "decision_curve")) ){
      stop("One or more elements of the list provided is not an object of class 'decision_curve' (output from the function 'decision_curve').")
    }
    policy.vec <- sapply(x, FUN = function(x) x$policy)
    if(length(unique(policy.vec))>1) stop("Comparing decision curves with different opt-in/opt-out policies is not valid.")
    if(is.na(confidence.intervals[1])){
      
      ci.list <- sapply(x, FUN = function(x) x$confidence.intervals)
      ci.list.log <- sapply(x, FUN = function(x) is.numeric(x$confidence.intervals))
      ci.list.num <- as.numeric(ci.list[ci.list.log])
      
      if(length(unique(ci.list.num))>1){warning("Confidence intervals of different sizes are being plotted on the same figure.")}
      
      confidence.intervals <- ifelse(any(ci.list.log), 1, "none")
    }else{
      confidence.intervals <- ifelse(confidence.intervals, 1, "none")
    }
    
    if(policy.vec[1] == "opt-in"){
      message("Note: When multiple decision curves are plotted, decision curves for 'All' are calculated using the prevalence from the first DecisionCurve object in the list provided.")
    }else{
      message("Note: When multiple decision curves are plotted, decision curves for 'None' are calculated using the prevalence from the first DecisionCurve object in the list provided.")
    }
    model <- NULL #appease check
    #multiple dc's to plot
    #pull the "all' and 'none' curves from the first element in x
    dc.data <- subset(x[[1]]$derived.data, is.element(model, c("All", "None")))
    
    #fill in ci variables for if confidence intervals weren't calculated using DecisionCurve
    if(ncol(dc.data) == 9 ) dc.data <- add.ci.columns(dc.data)
    
    predictors <- NULL
    
    #loop through the remaining curves
    i = 0
    for(ll in x){
      i = i + 1
      #extract data to add
      newdata <-  subset(ll$derived.data, !is.element(model, c("All", "None")))
      #predictor name
      
      if(is.na(curve.names[1])){
        #check to make sure the name is different
        newpred <- unique(newdata$model)
        if(is.element(newpred, predictors)) stop("After extracting the curve names from the decision_curve object, the names of the decision curves provided are the same for two or more decision_curve objects. Please set curve.names to avoid errors in plotting.")
      }else{
        newdata$model <- curve.names[i]
        newpred <- unique(newdata$model)
      }
      
      predictors <- c(predictors, newpred)
      
      #if confidence intervals weren't calculated
      if(ncol(newdata) == 9 ){
        if(confidence.intervals) warning(paste("confidence interval plotting were requested for curve '", newpred, "' but not calculated using decision_curve", sep = ''))
        #fill in ci variables for if confidence intervals weren't calculated using DecisionCurve
        newdata <- add.ci.columns(newdata)
      }
      
      dc.data <- rbind(dc.data, newdata)
    }
    
    
  }
  
  return(list(dc.data = dc.data,
              predictors = predictors,
              confidence.intervals = confidence.intervals))
  
}

plot_generic<- function(xx, predictors, value, plotNew,
                        standardize, confidence.intervals,
                        cost.benefit.axis = FALSE, cost.benefits, n.cost.benefits,
                        cost.benefit.xlab, xlab, ylab,
                        col, lty, lwd,
                        xlim, ylim, legend.position,
                        lty.fpr = 2, lty.tpr = 1,
                        tpr.fpr.legend = FALSE,
                        impact.legend = FALSE,
                        impact.legend.2 = FALSE,
                        population.size = 1000,
                        policy = policy, ...){
  ## xx is output from get_DecisionCurve,
  ## others are directly from the function call
  
  #save old par parameters and reset them once the function exits.
  old.par<- par("mar"); on.exit(par(mar = old.par))
  
  
  xx.wide <- reshape::cast(xx, thresholds~model, value =  value, add.missing = TRUE, fill = NA)
  xx.wide$thresholds <- as.numeric(as.character(xx.wide$thresholds))
  
  if(is.numeric(confidence.intervals)){
    
    val_lower <- paste(value, "lower", sep = "_")
    val_upper <- paste(value, "upper", sep = "_")
    
    xx.lower <- cast(xx, thresholds~model, value = val_lower, add.missing = TRUE, fill = NA)
    xx.upper <- cast(xx, thresholds~model, value = val_upper, add.missing = TRUE, fill = NA)
    xx.lower$thresholds <- as.numeric(as.character(xx.lower$thresholds))
    xx.upper$thresholds <- as.numeric(as.character(xx.upper$thresholds))
  }
  
  
  # adjust margins to add extra x-axis
  if(cost.benefit.axis) par(mar = c(7.5, 6, 3, 2) + 0.1)
  
  #set default ylim if not provided
  
  
  #initial call to plot and add gridlines
  if(plotNew){
    
    plot(xx.wide$thresholds, xx.wide$None, type = "n", ylim = ylim, axes=FALSE,
         col = "black", xlim = xlim,  xlab = xlab, ylab = ylab, frame.plot = FALSE, ...)
    
    grid(lty = 1, col = "grey92")
  }
  
  if(is.element(value, c("NB", "sNB"))){
    #plot none and all
    lines(xx.wide$thresholds, xx.wide$None, type = "l",
          col = col[length(predictors)+ 2],
          lty = lty[length(predictors)+ 2],
          lwd = lwd[length(predictors)+ 2])
    
    lines(xx.wide$threshold, xx.wide$All, type = "l",
          col = col[length(predictors)+ 1],
          lty = lty[length(predictors)+ 1],
          lwd = lwd[length(predictors)+ 1])
    
    if(is.numeric(confidence.intervals)){
      if(policy == "opt-in"){
        
        lines(xx.lower[,c("thresholds", "All")],
              col = col[length(predictors)+ 1],
              lty = lty[length(predictors)+ 1],
              lwd = lwd[length(predictors)+ 1]/2)
        
        lines(xx.upper[,c("thresholds", "All")],
              col = col[length(predictors)+ 1],
              lty = lty[length(predictors)+ 1],
              lwd = lwd[length(predictors)+ 1]/2)
        
        
      }else{
        lines(xx.lower[,c("thresholds", "None")],
              col = col[length(predictors)+ 2],
              lty = lty[length(predictors)+ 2],
              lwd = lwd[length(predictors)+ 2]/2)
        
        lines(xx.upper[,c("thresholds", "None")],
              col = col[length(predictors)+ 2],
              lty = lty[length(predictors)+ 2],
              lwd = lwd[length(predictors)+ 2]/2)
      }
      
    }
  }
  
  #the clinical impact plots are on a different scale
  if(is.element(value, c("DP" ,"nonDP", "prob.high.risk", "prob.low.risk"))){
    #population.size
    ps <- population.size
  }else{
    ps <- 1
  }
  
  #plot each predictor
  for(i in 1:length(predictors)){
    #plot ci's if asked for
    
    j <- ifelse(is.element(value, c("TPR", "prob.high.risk", "prob.low.risk")), 1, i)
    j <- ifelse(is.element(value, c("FPR", "DP", "nonDP")), 2, i)
    if(is.numeric(confidence.intervals)){
      #get rid of cases missing for that predictor; this sometimes
      #happens due to different thresholds for each predictor
      cc <- complete.cases(xx.lower[,c("thresholds", predictors[i])])
      
      lines(x = xx.lower[cc, c("thresholds")],
            y = xx.lower[cc, c(predictors[i])]*ps,
            type = "l",  col = col[j], lty = lty[i], lwd = lwd[i]/2)
      
      cc <- complete.cases(xx.upper[,c("thresholds", predictors[i])])
      
      lines(x = xx.upper[cc, c("thresholds")],
            y = xx.upper[cc, c(predictors[i])]*ps,
            type = "l",  col = col[j], lty = lty[i], lwd = lwd[i]/2)
      
    }
    cc <- complete.cases(xx.wide[,c("thresholds", predictors[i])])
    
    lines(x = xx.wide[cc, c("thresholds")],
          y = xx.wide[cc, c(predictors[i])]*ps,
          type = "l",  col = col[j], lty = lty[i], lwd = lwd[i])
    
    
  }
  
  #add legend
  if(is.element(legend.position, c("bottomright", "topright", "bottomleft", "topleft", "right", "left", "top", "bottom"))){
    
    if(value == "NB" | value == "sNB"){
      legend(legend.position, lty = lty, col = col, lwd = lwd, legend = c(predictors, "All", "None"),  text.font=2,
             ncol = 1, cex = 1.5, bg  = "white")
    }else if(tpr.fpr.legend){
      n.preds <- length(predictors)
      legend(legend.position,
             lty = c( lty.tpr, lty.fpr), bg  = "white",
             col = col,
             lwd = lwd, legend = c("True positive rate", "False positive rate"))
      
      
    } else if(impact.legend){
      legend(legend.position,
             lty = c( 1, 2),
             col = col, bg  = "white",
             lwd = lwd, legend = c("Number high risk", "Number high risk with event"))
      
    }else if(impact.legend.2){
      legend(legend.position,
             lty = c( 1, 2),
             col = col, bg  = "white",
             lwd = lwd, legend = c("Number low risk", "Number low risk without event"))
      
    }
    
  }
  axis(1, seq(0.6, 1, 0.1), seq(0.6, 1, 0.1), font.axis = 2, cex.axis = 1.5, lwd = 2)
  axis(2, seq(-0.1, 0.8, 0.1), seq(-0.1, 0.8, 0.1), font.axis = 2, cex.axis = 1.5, lwd = 2)
  # axis(3)
  # axis(4)
  box(lwd = 2)
  #add cost benefit axis if wanted
  # if(cost.benefit.axis){
  # 
  #   tmp <- Add_CostBenefit_Axis(xlim = xlim,
  #                               cost.benefits = cost.benefits,
  #                               n.cost.benefits = n.cost.benefits,
  #                               line = 4,
  #                               policy = policy)
  #   mtext(xlab, 1, 2.2)
  #   mtext(cost.benefit.xlab, side = 1, 6.1)
  # }else{
  #   mtext(xlab, side = 1, 3)
  # }
}


