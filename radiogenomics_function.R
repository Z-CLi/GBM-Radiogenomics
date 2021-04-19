# Radiogenomics analysis
# Date: 2021-04-19
# Authors: Qiuchang Sun and Zhi-Cheng Li


# ====================function1 import library ====================
import_library = function() {
  wants <- c("WGCNA", "impute",  "preprocessCore", "reshape2","stringr","dynamicTreeCut","fastcluster","plyr","clusterProfiler","org.Hs.eg.db","pathview","STRINGdb","ReactomePA","ggplot2",
           "DOSE","enrichplot","dplyr","GSVA","GSEABase","GSVAdata","Biobase","genefilter","limma","RColorBrewer","survival","glmnet","Matrix","foreach","gplots",
           "MASS","stringi","simPH","prodlim","survcomp","km.ci","pheatmap","circlize","ComplexHeatmap","unikn","colorspace")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  has <- wants %in% rownames(installed.packages())
  if(any(!has)) BiocManager::install(wants[!has])
  sapply(wants, require, character.only = TRUE)
  }

# ====================function2 wgcna_data ====================
# data cleaning and gene network construction
# ========input========
# rnaData: (row:gene; col:patient) log(FPKM+1)
# setPower: set 'TRUE', the code pick SoftThreshold(power),else use setted SoftThreshold(power)
# cutHeight :Choose a height cut to remove the offending sample
# resultDir: path for saving geneblock and figure 
# name: name for saving geneblock and figure 
# eg: setwd(resultDir)
# tiff(file = paste('./',name,'wgcna_clust.tiff',sep=''), res = 300, width = 9000, height = 5000, compression = "lzw")
# ========output========
# sampleTree: result of clustering the samples
# clust: Determine cluster under the line
# data: samples in the biggest clust
# power: SoftThreshold
# net: gene network
# MEs:The module eigengene E is defined as the first principal component of a given module. It can be considered a representative of the gene expression profiles in a module
# MM: module membership
# moduleColors,colors for every gene
# gsg$goodGenes,  result of removing the genes with too many missing value
# ========figure========
# 1.pick SoftThreshold
# 2.dendrograms
wgcnaData <- function(rnaData,setPower,cutHeight,resultDir,name) {
  gsg = goodSamplesGenes(rnaData, verbose = 3);
  gsg$allOK
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(rnaData)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(rnaData)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    rnaData = rnaData[gsg$goodSamples, gsg$goodGenes]
  }
  sampleTree = hclust(dist(rnaData), method = "average")
  tiff(file = file.path(resultDir,name,'wgcna_clust.tiff',sep=''), res = 300, width = 9000, height = 5000, compression = "lzw")
  par(cex=0.6);
  par(mar=c(0,4,2,0))
  plot(sampleTree,main="Sample clustering to detect outliers", sub="",xlab="",
       cex.lab=2.5, cex.axis=2.5, cex.main=3,font.axis=3,cex.axis=2.5,font.lab=3,cex.lab=2.5)
  abline(h=cutHeight,col="red");
  dev.off()
  clust=cutreeStatic(sampleTree, cutHeight= cutHeight ,minSize = 10)
  keepSamples = (clust==1)
  data = rnaData[keepSamples, ]

  if (setPower == 'TRUE'){
    print("wait a second, pick power")
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(data, powerVector=powers,verbose=5)
    tiff(file = file.path(resultDir,name,'wgcna_power.tiff',sep=''), res = 300, width = 8000, height = 4000, compression = "lzw")
    par(mfrow = c(1,2))
    cex1 = 0.9
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")
    abline(h=0.9,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
         cex=cex1, col="red")
    dev.off()
    power = sft$powerEstimate
  }
  else {
    print(paste("power is",setPower))
    power = setPower}
  cor = WGCNA::cor
  net = blockwiseModules(data, maxBlockSize = 8000, power = power, TOMType = "unsigned", 
                         minModuleSize = 50, 
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.15, numericLabels = TRUE, 
                         pamRespectsDendro = FALSE, saveTOMs = TRUE, 
                         saveTOMFileBase = file.path(resultDir,name,"GBMTOM",sep=""), verbose = 3)
  
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  MEs = moduleEigengenes(data, moduleColors)$eigengenes
  MEs = orderMEs(MEs)
  MEs = MEs[, -which(colnames(MEs) %in% c("MEgrey"))]
  MM = as.data.frame(cor(data, MEs, use ="p"))
  tiff(file = file.path(resultDir,name,'dendrograms.tiff'), res = 300, width = 15000, height = 5000, compression = "lzw")
  color = moduleColors[net$blockGenes[[1]]]
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Cluster Dendrogram",
                      font.main=3,cex.main=2.5,
                      cex.colorLabels=2.2, font.axis=3,cex.axis=2.5,font.lab=3,cex.lab=2.5)
  dev.off()
  
  
  return(list(sampleTree,clust,data,power,net,MEs,MM,moduleColors,gsg$goodGenes))
}
# ====================function3 zsummary ====================
#  module preservation between two independent data sets
# ========input========
# dataTrain: discovery dataset (output[3] from function1 wgcna_data)
# dataTest: valid datasets (output[3] from function1 wgcna_data)
# colorTrain: color for every gene in discovery dataset (output[8] from function1 wgcna_data)
# colorTest: color for every gene in valid dataset (output[8] from function1 wgcna_data)
# resultDir,name: same with function1
# ========output========
# zsummaryResult: the preservation medianRank and Zsummary statistics for both datasets
# mp: the function of module preservation
# ========figure========
# 1.zsummary.tiff

zsummaryData <- function(dataTrain, dataTest, colorTrain, colorTest,resultDir,name) {
  setLabels = c("dataTrain", "dataTest");
  multiExpr = list(dataTrain = list(data = dataTrain), dataTest = list(data = dataTest));
  multiColor = list(dataTrain = colorTrain, dataTest = colorTest);
  nSets = 2
  system.time( {
    mp = modulePreservation(multiExpr, multiColor,
                            referenceNetworks = c(1:2),
                            nPermutations = 200,
                            randomSeed = 1,
                            verbose = 3)
  } );
  
  # save(mp, file = "BxHLiverFemaleOnly-modulePreservation.RData");
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
  print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
  zsummaryResult = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                          signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
  # Module labels and module sizes are also contained in the results
  modColors = rownames(mp$preservation$observed[[ref]][[test]])
  moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
  # leave grey and gold modules out
  plotMods = !(modColors %in% c( "qqq"));
  # !(modColors %in% c("grey", "gold"));
  # Text labels for points
  text = modColors[plotMods];
  # Auxiliary convenience variable
  plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
  # Main titles for the plot
  mains = c("Preservation Median rank", "Preservation Zsummary");
  # Start the plot
  tiff(file = file.path(resultDir,name,'zsummary.tiff',sep=''), res = 300, width = 10000, height = 5000, compression = "lzw")

  #pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
  par(mfrow = c(1,2))
  par(mar = c(4.5,4.5,2.5,1))
  for (p in 1:2)
  {
    min = min(plotData[, p], na.rm = TRUE);
    max = max(plotData[, p], na.rm = TRUE);
    # Adjust ploting ranges appropriately
    if ((p==1)|(p==2))
    {
      if (min > -max/10) min = -max/10
      ylim = c(min - 0.01 * (max-min), max + 0.01 * (max-min))
    } else
      ylim = c(max + 0.01 * (max-min), min - 0.01 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
         main = mains[p],
         cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = ylim,
         xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    # labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
    # labelPoints(moduleSizes[plotMods][2], plotData[plotMods, p][2], text[2],offs = 0.08, cex = 1);
    text(moduleSizes[plotMods][2], plotData[plotMods, p][2], text[2],pos=1,cex = 1.5);
    # For Zsummary, add threshold lines
    if (p==2)
    {
      abline(h=0)
      abline(h=2, col = "blue", lty = 2)
      abline(h=10, col = "darkgreen", lty = 2)
    }
  }
  dev.off()
  return(list(zsummaryResult,mp))
}
# ====================function3 makeGsvaSet ====================
# make GSVA set using gene moudle
# ========input========
# color: color for modules
# moduleColors: color for genes
# geneName: genes names
# ========output========
# module gene set
makeGsvaSet <-function(color,moduleColors,geneName){
  list1 = list()
  for(w in 1:length(color)){
    a =  geneName[which(moduleColors == color[w])]
    list1[[w]] = data.frame(t(a))}
  return(list(list1))
}
# ====================function4 genesetPathwayModule ====================
# A pathway enrichment analysis within each radiomics-associated gene module identified
# ========input========
# color: identified gene module color
# moduleColors: colors for every gene
# geneName: genes names
# MEs: same with function1 wgcna_data output
# MM: same with function1 wgcna_data output
# name: same with function1 wgcna_data 
# resultDir: same with function1 wgcna_data
# ========output========
# pathway enrichment results (GO KEGG REACTOME HALLMARK CPG CP PID) in csv

genesetPathwayModule <-function(color,moduleColors,geneName,MEs,MM,name,resultDir) {
  for(w in 1:length(color)){
    if(!(file.exists(file.path(resultDir,color[w])))){
      dir.create(file.path(resultDir,color[w]))}
    setwd(file.path(resultDir,color[w]))
    a =  geneName[which(moduleColors == color[w])]
    modNames = substring(names(MEs), 3)
    column = match(color[w], modNames);
    print(color[w])
    print(column)
    ss = data.frame(a,MM=MM[a, column])
    pos = a[which(ss$MM>0)]
    neg = a[which(ss$MM<0)]
    print(length(pos))
    print(length(neg))
    # set the number 20 to make sure the module genes can enrich pathway
    if (length(pos)>20){
      genelistPos <- bitr(pos,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db) # groupGOĬ??ʹ??"ENTREZID"??Ϊ???????ݣ????????Ƚ???SYMBOL??ת??Ϊ??ENTREZID?????Է��???
      
      goPos <-enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                        qvalueCutoff = 0.2, TERM2GENE = geneset_GO)
      keggPos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                           qvalueCutoff = 0.2, TERM2GENE = geneset_KEGG)
      reactPos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                            qvalueCutoff = 0.2, TERM2GENE = geneset_REACTOME)
      
      hallmarkPos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                               qvalueCutoff = 0.2, TERM2GENE = geneset_hallmark)
      # biocarta_pos <- enricher(genelist_pos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 10, maxGSSize = 500,
      # qvalueCutoff = 0.2, TERM2GENE = geneset_biocarta)
      CPG_pos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                          qvalueCutoff = 0.2, TERM2GENE = geneset_CPG)
      CP_pos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                         qvalueCutoff = 0.2, TERM2GENE = geneset_CP)
      PID_pos <- enricher(genelistPos$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                          qvalueCutoff = 0.2, TERM2GENE = geneset_PID)

      fileGoPos = paste(name,color[w],"goPos",".csv",sep = "")
      write.csv(as.data.frame(goPos),fileGoPos,row.names =FALSE)
      fileKeggPos = paste(name,color[w],"keggPos",".csv",sep = "")
      write.csv(as.data.frame(keggPos),fileKeggPos,row.names =FALSE)
      fileReactPos = paste(name,color[w],"reactomePos",".csv",sep = "")
      write.csv(as.data.frame(reactPos),fileReactPos,row.names =FALSE)
      fileHallPos = paste(name,color[w],"hallmarkPos",".csv",sep = "")
      write.csv(as.data.frame(hallmarkPos),fileHallPos,row.names =FALSE)
      # filebiopos = paste(name,color[w],"biocarta_pos",".csv",sep = "")
      # write.csv(summary(biocarta_pos),filebiopos,row.names =FALSE)
      fileCpgPos = paste(name,color[w],"CPG_pos",".csv",sep = "")
      write.csv(as.data.frame(CPG_pos),fileCpgPos,row.names =FALSE)
      fileCpPos = paste(name,color[w],"CP_pos",".csv",sep = "")
      write.csv(as.data.frame(CP_pos),fileCpPos,row.names =FALSE)
      filePidPos = paste(name,color[w],"PID_pos",".csv",sep = "")
      write.csv(as.data.frame(PID_pos),filePidPos,row.names =FALSE)
    }
    
    if (length(neg)>20){
      genelistNeg <- bitr(neg,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db) 
      
      goNeg <-enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                        qvalueCutoff = 0.2, TERM2GENE = geneset_GO)
      keggNeg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                           qvalueCutoff = 0.2, TERM2GENE = geneset_KEGG)
      reactNeg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                            qvalueCutoff = 0.2, TERM2GENE = geneset_REACTOME)
      
      hallmarkNeg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                               qvalueCutoff = 0.2, TERM2GENE = geneset_hallmark)
      # biocarta_neg <- enricher(genelist_neg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 10, maxGSSize = 500,
      #                          qvalueCutoff = 0.2, TERM2GENE = geneset_biocarta)
      CPG_neg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                          qvalueCutoff = 0.2, TERM2GENE = geneset_CPG)
      CP_neg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                         qvalueCutoff = 0.2, TERM2GENE = geneset_CP)
      PID_neg <- enricher(genelistNeg$SYMBOL,pvalueCutoff = 0.1, pAdjustMethod = "BH",  minGSSize = 5, maxGSSize = 1000,
                          qvalueCutoff = 0.2, TERM2GENE = geneset_PID)
      
      fileGoNeg = paste(name,color[w],"goNeg",".csv",sep = "")
      write.csv(as.data.frame(goNeg),fileGoNeg,row.names =FALSE)
      fileKeggNeg = paste(name,color[w],"keggNeg",".csv",sep = "")
      write.csv(as.data.frame(keggNeg),fileKeggNeg,row.names =FALSE)
      fileReactNeg = paste(name,color[w],"reactomeNeg",".csv",sep = "")
      write.csv(as.data.frame(reactNeg),fileReactNeg,row.names =FALSE)
      fileHallNeg = paste(name,color[w],"hallmarkNeg",".csv",sep = "")
      write.csv(as.data.frame(hallmarkNeg),fileHallNeg,row.names =FALSE)
      # filebioneg = paste(name,color[w],"biocarta_neg",".csv",sep = "")
      # write.csv(summary(biocarta_neg),filebioneg,row.names =FALSE)
      fileCpgNeg = paste(name,color[w],"CPG_neg",".csv",sep = "")
      write.csv(as.data.frame(CPG_neg),fileCpgNeg,row.names =FALSE)
      fileCpNeg = paste(name,color[w],"CP_neg",".csv",sep = "")
      write.csv(as.data.frame(CP_neg),fileCpNeg,row.names =FALSE)
      filePidNeg = paste(name,color[w],"PID_neg",".csv",sep = "")
      write.csv(as.data.frame(PID_neg),filePidNeg,row.names =FALSE)}
  }
}

# ====================function5 makePathFile ====================
# Merge the CSV within the folder
# ========input========
# file_path: the forder 
# ========output========
# result of merging 

makePathFile = function(filePath){
  setwd(filePath)
  file = list.files(filePath)
  n = length(file) 
  file1 = read.csv(file[1], header = TRUE)
  for (i in 2:n){
    file2 = read.csv(file[i], header = TRUE)
    file1 = rbind(file1,file2)
  }
  return(file1)
}


# ====================function6 keyGene ====================
# genes with |MM| > 0.80 and GS > 0.20 were selected as key genes
# ========input========
# data: data function1 wgcna_data output
# sig: radiomics siganture
# module: pick a module (eg "blue")
# moduleColor: colors for every gene
# MM: same with function1 wgcna_data output
# ========output========
# a list of keygene

keyGene <- function(data, sig, module, moduleColor,MM){
  geneSig = as.data.frame(cor(data,sig, method= 'pearson'))
  geneSigPvalue = as.data.frame(corPvalueStudent(as.matrix(geneSig), length(sig)))
  modNames = substring(names(MM), 3)
  column = match(module, modNames);
  moduleGenes = moduleColor==module;
  geneKey = data.frame(name = colnames(data)[moduleGenes], MM = abs(MM[moduleGenes, column]), geneSig = abs(geneSig[moduleGenes, 1]), geneSigPvalue = geneSigPvalue[moduleGenes, 1])
  geneKey = geneKey[order(-geneKey$geneSig),]
  geneKey = geneKey[geneKey$MM>0.8,]
  geneKey = geneKey[geneKey$geneSigPvalue<0.2,]
  return(geneKey)
}


# ====================function7 draw km ====================
# draw km curve
# ========input========
# status, time(os), signature, cutoff, name
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
drawKM.Rad_2 <- function(status, time, sig, cutoff, name) {
  
  Status <- as.double(status)
  Time <- as.double(time)
  
  stratRad = sig
  stratRad[which(sig >= cutoff)]=1
  stratRad[which(sig < cutoff)]=0
  
  dd1 <- data.frame("surv.time"=Time, "surv.event"=Status, "strat"=stratRad)
  tiff(file = paste('./',name,'_figure_km.tiff',sep=''), res = 600, width = 6000, height = 4000, compression = "lzw")
  par(mar=c(0,0,0,0))
  km.coxph.plot_1(formula.s=Surv(Time, Status) ~ stratRad, data.s=dd1,
                  x.label="Time (days)", y.label="Probability of survival",mark.time = T,pch="+",cex=2,.lwd = 4,
                  main.title="",.col=c("darkgreen", "red3"), n.risk.step = 200,
                  .lty=c(1,1), show.n.risk=T, n.risk.cex=0.85, verbose=T,
                  leg.text=paste(c("Low risk", "High risk"), "   ", sep=""),leg.bty = "n",
                  leg.pos="topright",font.axis=2,cex.axis=1.2,font.lab=2,cex.lab=1.2,font.legend=2)
  dev.off()
  
  coxphSig=coxph(Surv(Time,Status)~sig)
  coxphStratRad=coxph(Surv(Time,Status)~stratRad)
  print(summary(coxphSig))
  print(summary(coxphStratRad))
}
 
# ==================== function8 Heatmap of z-score-normalized gene expressions ====================
# Heatmap of z-score-normalized gene expressions of two highly reproducible radiomics-associated modules (red and blue)
# ========input========
# radioGenTrain_heatmap1_gene_expressions
# externalTest_heatmap1_gene_expressions
# resultDir
# ========output========
# heatmap_gene.tiff

heatmapGene <- function(radioGenTrain_heatmap1_gene_expressions,externalTest_heatmap1_gene_expressions,resultDir){
  data = radioGenTrain_heatmap1_gene_expressions
  radSig <- as.numeric(data[,1])
  survTime <- as.numeric(data[,2])
  survStatus <- as.numeric(data[,3])
  gapsCol = 53
  group = c(rep("1",gapsCol),rep("0",nrow(data)-gapsCol))
  dataRGT = as.data.frame(t(data[,4:ncol(data)]))
  
  rowSpiltModule = c(rep("Red Module",446),rep("Blue Module",369))
  colorGene <- colorRampPalette(c("navy", "white", "firebrick3"))(256)
  colorSig <- colorRamp2(c(max(radSig),median(radSig),min(radSig)),rev(brewer.pal(9, "Blues"))[c(1,3,6)])
  colorOs <- colorRamp2(c(max(survTime),median(survTime),min(survTime)),usecol(pal = rev(pal_petrol),3))
  
  ha_0 = HeatmapAnnotation(group = group,
                           sig = radSig,
                           os = survTime,
                           status = survStatus,
                           border = TRUE,
                           # gp = gpar(col = "grey"),
                           col = list(group = c("1" = "red", "0" = "forestgreen"),
                                      sig = colorSig,
                                      os = colorOs,
                                      status =  c("1" = "white", "0" = "dimgrey")),
                           gap = unit(c(rep(0.5,3),2), "mm"),
                           show_annotation_name = FALSE,
                           annotation_height =   unit(c(1, 1, 1,0.1), c("mm", "mm","mm","mm", "mm")),
                           annotation_legend_param = list(group = list(at = c(0, 1),
                                                                       labels = c("low", "high"),
                                                                       title = "risk group",
                                                                       gp = gpar(fontsize = 16,fontface = "bold"),
                                                                       direction = "horizontal"),
                                                          sig = list(title = "RadRisk score",
                                                                     gp = gpar(fontsize = 16,fontface = "bold"),
                                                                     legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          os = list(title = "overall survival (days)",
                                                                    gp = gpar(fontsize = 16,fontface = "bold"),
                                                                    legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          status = list(at = c(0, 1),
                                                                        labels = c("censored", "uncensored"),
                                                                        title = "survival status",gp = gpar(fontsize = 16,fontface = "bold"),
                                                                        border = TRUE,
                                                                        direction = "horizontal")))
  redHighRGT = apply(dataRGT[1:446,1:53],2,mean)
  redLowRGT = apply(dataRGT[1:446,54:93],2,mean)
  blueHighRGT = apply(dataRGT[447:815,1:53],2,mean)
  blueLowRGT = apply(dataRGT[447:815,54:93],2,mean)
  m = c("Red Module","Red Module","Blue Module","Blue Module")
  meanValueRGT = c(mean(redHighRGT),mean(redLowRGT),mean(blueHighRGT),mean(blueLowRGT))
  meanRGT = data.frame(m = m,value = meanValueRGT)
  panel_fun_g2 = function(index, nm) {
    pushViewport(viewport(xscale = c(0, 2), yscale = c(0, 2)))
    ht = anno_barplot(meanRGT[which(meanRGT[,1]==nm),2],baseline = 0,height = unit(4, "cm"),
                      bar_width = 0.5,gp = gpar(fill = c("red", "forestgreen"),fontface = "bold"),
                      axis_param = list(
                        side = "right"),
                      show_legend = TRUE)
    g = grid.grabExpr(draw(ht))
    grid.rect()
    grid.draw(g)
    popViewport()
  }
  anno_g2 = anno_zoom(align_to = rowSpiltModule, which = "row", panel_fun = panel_fun_g2, 
                      size = unit(4, "cm"), gap = unit(1, "cm"), width = unit(6, "cm"))
  # Heatmap(m, name = "mat", left_annotation = rowAnnotation(foo = anno), row_split = subgroup)
  colnames(dataRGT) = seq(1:93)
  heat1_RGT = Heatmap(as.matrix(dataRGT), name = "gene2 ",show_heatmap_legend = FALSE,
                    column_title = "RDC",
                    column_title_gp = gpar(fill = "grey", col = "black", border = "grey", 
                                           fontsize = 12,fontface = "bold", height=unit(0.01, "mm")),
                    col = colorGene,
                    show_row_names = FALSE,
                    # row_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(direction = "horizontal"),
                    width = unit(16, "cm"), height = unit(12, "cm"),
                    cluster_rows = FALSE, show_column_names = FALSE,
                    column_split = factor(group, levels = c("1", "0")),
                    row_split = factor(rowSpiltModule, levels = c("Red Module", "Blue Module")),
                    column_gap = unit(5, "mm"),
                    row_gap = unit(5, "mm"),
                    cluster_columns = FALSE,top_annotation = ha_0, right_annotation = rowAnnotation(foo2 = anno_g2))
  
  data = externalTest_heatmap1_gene_expressions
  radSig <- as.numeric(data[,1])
  survTime <- as.numeric(data[,2])
  survStatus <- as.numeric(data[,3])
  gapsCol = 39
  group = c(rep("1",gapsCol),rep("0",nrow(data)-gapsCol))
  dataET = as.data.frame(t(data[,4:ncol(data)]))
  
  colorSig <- colorRamp2(c(max(radSig),median(radSig),min(radSig)),rev(brewer.pal(9, "Blues"))[c(1,3,6)])
  colorOs <- colorRamp2(c(max(survTime),median(survTime),min(survTime)),usecol(pal = rev(pal_petrol),3))
  
  ha_0 = HeatmapAnnotation(group = group,
                           sig = radSig,
                           os = survTime,
                           status = survStatus,
                           border = TRUE,
                           # gp = gpar(col = "grey"),
                           col = list(group = c("1" = "red", "0" = "forestgreen"),
                                      sig = colorSig,
                                      os = colorOs,
                                      status =  c("1" = "white", "0" = "dimgrey")),
                           gap = unit(c(rep(0.5,3),2), "mm"),
                           show_annotation_name = FALSE,
                           annotation_height =   unit(c(1, 1, 1, 1, 0.1), c("mm", "mm","mm","mm", "mm")),
                           annotation_legend_param = list(group = list(at = c(0, 1),
                                                                       labels = c("low", "high"),
                                                                       title = "risk group",
                                                                       gp = gpar(fontsize = 16,fontface = "bold"),
                                                                       direction = "horizontal"),
                                                          sig = list(title = "RadRisk score",
                                                                     gp = gpar(fontsize = 16,fontface = "bold"),
                                                                     legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          os = list(title = "overall survival(days)",
                                                                    gp = gpar(fontsize = 16,fontface = "bold"),
                                                                    legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          status = list(at = c(0, 1),
                                                                        labels = c("censored", "uncensored"),
                                                                        title = "survival status",gp = gpar(fontsize = 16,fontface = "bold"),
                                                                        border = TRUE,
                                                                        direction = "horizontal")))
  
  redHighET = apply(dataET[1:446,1:39],2,mean)
  redLowET = apply(dataET[1:446,40:76],2,mean)
  blueHighET = apply(dataET[447:815,1:39],2,mean)
  blueLowET = apply(dataET[447:815,40:76],2,mean)
  m = c("Red Module","Red Module","Blue Module","Blue Module")
  meanValueET = c(mean(redHighET),mean(redLowET),mean(blueHighET),mean(blueLowET))
  meanET = data.frame(m = m,value = meanValueET)
  panel_fun_g3 = function(index, nm) {
    pushViewport(viewport(xscale = c(0, 2), yscale = c(0, 2)))
    ht = anno_barplot(meanET[which(meanET[,1]==nm),2],baseline = 0,height = unit(4, "cm"),
                      bar_width = 0.5,gp = gpar(fill = c("red", "forestgreen"),fontface = "bold"),
                      axis_param = list(
                        side = "right"))
    g = grid.grabExpr(draw(ht))
    grid.rect()
    grid.draw(g)
    popViewport()
  }
  anno_g3 = anno_zoom(align_to = rowSpiltModule, which = "row", panel_fun = panel_fun_g3, 
                      size = unit(4, "cm"), gap = unit(1, "cm"), width = unit(6, "cm"))
  colnames(dataET) = seq(1:76)
  
  heat1_ET = Heatmap(as.matrix(dataET), name = "gene expression",
                    column_title = "RVC",
                    column_title_gp = gpar(fill = "grey", col = "black", border = "grey", 
                                           fontsize = 12,fontface = "bold", height=unit(0.01, "mm")),
                    col = colorGene,
                    show_row_names = FALSE,
                    # row_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(direction = "horizontal"),
                    width = unit(13, "cm"), height = unit(12, "cm"),
                    cluster_rows = FALSE, show_column_names = FALSE,
                    column_split = factor(group, levels = c("1", "0")),
                    row_split = factor(rowSpiltModule, levels = c("Red Module", "Blue Module")),
                    row_title = c("red module", "blue module"),
                    row_title_rot = 0,
                    column_gap = unit(5, "mm"),
                    row_gap = unit(5, "mm"),
                    cluster_columns = FALSE,top_annotation = ha_0, right_annotation = rowAnnotation(gene_expression = anno_g3))
  
  heat1_list = heat1_RGT + heat1_ET 
  
  setwd(resultDir)
  tiff(file = paste('./','heatmap_gene.tiff',sep=''), res = 600, width = 13000, height = 5000, compression = "lzw")
  par(mar=c(0,0,0,0))

  draw(heat1_list,merge_legend = TRUE,ht_gap = unit(25, "mm"),
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  annotation_titles = c(group = "risk group",
                        sig = "RadRisk score",
                        os = "overall survival",
                        status = "survival status")

  for(an in names(annotation_titles)) {
    decorate_annotation(an, {
      grid.text(annotation_titles[an], unit(-1.38,'npc'), just = "left",gp = gpar(fontsize = 12,fontface = "bold"))
    })
  }
  dev.off()
}

# ====================function10 pathway GSVA score across RDC and RVC patients ====================
# 13 imaging features and their significantly associated pathways
# ========input========
# radioGenTrain_heatmap2_GSVA_features,externalTest_heatmap2_GSVA_features,resultDir
# ========output========
# figureHeatmapGSVA.tiff

heatmap2PathwayGSVA <- function(radioGenTrain_heatmap2_GSVA_features,externalTest_heatmap2_GSVA_features,resultDir){
  dataRGT = radioGenTrain_heatmap2_GSVA_features
  dataET = externalTest_heatmap2_GSVA_features
  radSig <- as.numeric(c(dataRGT[,1],dataET[,1]))
  survTime <- as.numeric(c(dataRGT[,2],dataET[,2]))
  survStatus <- c(as.numeric(dataRGT[,3]),as.numeric(dataET[,3]))
  gapsCol = 53
  group = c(rep("A",53),rep("B",40),rep("C",39),rep("D",37))
  feature1 <- as.numeric(c(myNormalize(dataRGT[,4]),myNormalize(dataET[,4])))
  feature2 <- as.numeric(c(myNormalize(dataRGT[,5]),myNormalize(dataET[,5])))
  feature3 <- as.numeric(c(myNormalize(dataRGT[,6]),myNormalize(dataET[,6])))
  feature4 <- as.numeric(c(myNormalize(dataRGT[,7]),myNormalize(dataET[,7])))
  feature5 <- as.numeric(c(myNormalize(dataRGT[,8]),myNormalize(dataET[,8])))
  feature6 <- as.numeric(c(myNormalize(dataRGT[,9]),myNormalize(dataET[,9])))
  feature7 <- as.numeric(c(myNormalize(dataRGT[,10]),myNormalize(dataET[,10])))
  feature8 <- as.numeric(c(myNormalize(dataRGT[,11]),myNormalize(dataET[,11])))
  feature9 <- as.numeric(c(myNormalize(dataRGT[,12]),myNormalize(dataET[,12])))
  feature10 <- as.numeric(c(myNormalize(dataRGT[,13]),myNormalize(dataET[,13])))
  feature11 <- as.numeric(c(myNormalize(dataRGT[,14]),myNormalize(dataET[,14])))
  feature12 <- as.numeric(c(myNormalize(dataRGT[,15]),myNormalize(dataET[,15])))
  feature13 <- as.numeric(c(myNormalize(dataRGT[,16]),myNormalize(dataET[,16])))
  dataRGT_GSVA = dataRGT[,17:ncol(dataRGT)]
  dataRGT_GSVA = as.data.frame(t(dataRGT_GSVA))
  dataET_GSVA = dataET[,17:ncol(dataET)]
  dataET_GSVA = as.data.frame(t(dataET_GSVA))
  data = cbind(dataRGT_GSVA,dataET_GSVA)
  feature1Pathway = data[1:5,]
  feature2Pathway = data[6:10,]
  feature3Pathway = data[11:15,]
  feature4Pathway = data[16:20,]
  feature5Pathway = data[21:25,]
  feature6Pathway = data[26:30,]
  feature7Pathway = data[31:35,]
  feature8Pathway = data[36:40,]
  feature9Pathway = data[41:45,]
  feature10Pathway = data[46:50,]
  feature11Pathway = data[51:55,]
  feature12Pathway = data[56:60,]
  feature13Pathway = data[61:65,]
  
  colorGsva <- colorRampPalette(rev(divergingx_hcl(100,"RdBu")))(256)
  colorSig <- colorRamp2(c(max(radSig),median(radSig),min(radSig)),rev(brewer.pal(9, "Blues"))[c(1,3,6)])
  colorOs <- colorRamp2(c(max(survTime),median(survTime),min(survTime)),usecol(pal = rev(pal_petrol),3))
  
  colorFeature2 = colorRamp2(c(max(feature2),median(feature2),min(feature2)),sequential_hcl(3,"Burg"))
  colorFeature3 = colorRamp2(c(max(feature3),median(feature3),min(feature3)),sequential_hcl(3,"Burg"))
  colorFeature5 = colorRamp2(c(max(feature5),median(feature5),min(feature5)),sequential_hcl(3,"Burg"))
  colorFeature12 = colorRamp2(c(max(feature12),median(feature12),min(feature12)),sequential_hcl(3,"Burg"))
  colorFeature13 = colorRamp2(c(max(feature13),median(feature13),min(feature13)),sequential_hcl(3,"Burg"))
  
  colorFeature1 = colorRamp2(c(max(feature1),median(feature1),min(feature1)),sequential_hcl(7,"Purples")[c(1,3,5)])
  colorFeature9 = colorRamp2(c(max(feature9),median(feature9),min(feature9)),sequential_hcl(7,"Purples")[c(1,3,5)])
  
  colorFeature4 = colorRamp2(c(max(feature4),median(feature4),min(feature4)),sequential_hcl(7,"YlOrRd")[4:6])
  colorFeature6 = colorRamp2(c(max(feature6),median(feature6),min(feature6)),sequential_hcl(7,"YlOrRd")[4:6])
  
  colorFeature7 = colorRamp2(c(max(feature7),median(feature7),min(feature7)),sequential_hcl(3,"Green-Yellow"))
  colorFeature8 = colorRamp2(c(max(feature8),median(feature8),min(feature8)),sequential_hcl(3,"Green-Yellow"))
  colorFeature10 = colorRamp2(c(max(feature10),median(feature10),min(feature10)),sequential_hcl(3,"Green-Yellow"))
  colorFeature11 = colorRamp2(c(max(feature11),median(feature11),min(feature11)),sequential_hcl(3,"Green-Yellow"))
  
  
  ha_1 = HeatmapAnnotation(feature1 = feature1,col = list(feature1 = colorFeature1),show_annotation_name = FALSE,border = FALSE,
                           annotation_legend_param = list(feature1 = list(title = "proliferative",column_gap = unit(5, "cm"),legend_width = unit(5, "cm"),gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 20,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 120,fontface = "bold"),annotation_height = unit(0.1, "mm")) 
  ha_3 = HeatmapAnnotation(feature3 = feature3,col = list(feature3 = colorFeature3),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature3 = list(title = expression(italic(f[3])),column_gap = unit(5, "cm"),legend_width = unit(5, "cm"),gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(0.1, "mm")) 
  ha_4 = HeatmapAnnotation(feature4 = feature4,col = list(feature4 = colorFeature4),show_annotation_name = FALSE,border = FALSE, 
                           annotation_legend_param = list(feature4 = list(title = "treatment responsive",column_gap = unit(5, "cm"),legend_width = unit(5, "cm"),gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 20,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(0.1, "mm")) 
  ha_5 = HeatmapAnnotation(feature5 = feature5,col = list(feature5 = colorFeature5),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature5 = list(title = expression(italic(f[5])),legend_width = unit(5, "cm"),gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(0.1, "mm")) 
  ha_6 = HeatmapAnnotation(feature6 = feature6,col = list(feature6 = colorFeature6),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature6 = list(title = expression(italic(f[6])),legend_width = unit(5, "cm"),gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),height = unit(0.1, "mm")) 
  ha_7 = HeatmapAnnotation(feature7 = feature7,col = list(feature7 = colorFeature7),show_annotation_name = FALSE,border = FALSE,
                           annotation_legend_param = list(feature7 = list(title = "celluar function",legend_width = unit(5, "cm"),row_gap = unit(40, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 20,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(0.1, "mm")) 
  ha_8 = HeatmapAnnotation(feature8 = feature8,col = list(feature8 = colorFeature8),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature8 = list(title = expression(italic(f[8])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),
                           annotation_height = unit(0.1, "mm")) 
  ha_9 = HeatmapAnnotation(feature9 = feature9,col = list(feature9 = colorFeature9),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature9 = list(title = expression(italic(f[9])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                    title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(1, "mm")) 
  ha_10 = HeatmapAnnotation(feature10 = feature10,col = list(feature10 = colorFeature10),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                            annotation_legend_param = list(feature10 = list(title = expression(italic(f[10])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                      title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                            annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(1, "mm")) 
  ha_11 = HeatmapAnnotation(feature11 = feature11,col = list(feature11 = colorFeature11),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                            annotation_legend_param = list(feature11 = list(title = expression(italic(f[11])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                      title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                            annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(1, "mm")) 
  
  ha_12= HeatmapAnnotation(feature12 = feature12,col = list(feature12 = colorFeature12),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                           annotation_legend_param = list(feature12 = list(title = expression(italic(f[12])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                     title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                           annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(1, "mm")) 
  ha_13 = HeatmapAnnotation(feature13 = feature13,col = list(feature13 = col_feature13),show_annotation_name = FALSE,border = FALSE,show_legend = FALSE,
                            annotation_legend_param = list(feature13 = list(title = expression(italic(f[13])),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize =18,fontface = "bold"),
                                                                      title_gp = gpar(fontsize = 25,fontface = "bold"),direction = "horizontal")),
                            annotation_name_gp = gpar(fontsize = 20,fontface = "bold"),annotation_height = unit(1, "mm")) 
  
  ha_2 = columnAnnotation(risk_group = group,
                          RadRisk_score = radSig,
                          overall_survival = survTime,
                          survival_status = survStatus,
                          feature2 = feature2,
                          border = c(risk_group=TRUE,RadRisk_score=TRUE,overall_survival=TRUE,survival_status=TRUE),
                          # gp = gpar(col = "grey"),
                          col = list(risk_group = c("A" = "red", "B" = "forestgreen","C" = "red", "D" = "forestgreen"),
                                     RadRisk_score = colorSig,
                                     overall_survival = colorOs,
                                     survival_status =  c("1" = "white", "0" = "dimgrey"),
                                     f2 = colorFeature2),
                          show_annotation_name = FALSE,
                          gap = unit(c(rep(0.5,3),2), "mm"),
                          annotation_height =   unit(c(4, 4, 4, 4, 0.1), c("mm", "mm","mm","mm", "mm")),
                          annotation_legend_param = list(risk_group = list(at = c("A", "B"),
                                                                           labels = c("high", "low"),
                                                                           title = "risk group",
                                                                           title_gp = gpar(fontsize = 20,fontface = "bold"),
                                                                           labels_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                           legend_height = unit(2, "cm"),
                                                                           direction = "horizontal"),
                                                         RadRisk_score = list(title = "RadRisk score",
                                                                              title_gp = gpar(fontsize = 20,fontface = "bold"),labels_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                              legend_width = unit(5, "cm"),direction = "horizontal"),
                                                         overall_survival = list(title = "overall survival(days)",at = c(0,1000,2000),
                                                                                 title_gp = gpar(fontsize = 20,fontface = "bold"),labels_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                                 legend_width = unit(5, "cm"),direction = "horizontal"),
                                                         survival_status = list(at = c(0, 1),column_gap = unit(5, "mm"),
                                                                                labels = c("censored", "uncensored"),
                                                                                title = "survival status",title_gp = gpar(fontsize = 20,fontface = "bold"),
                                                                                labels_gp = gpar(fontsize = 15,fontface = "bold"),
                                                                                border = TRUE,column_gap = unit(5, "mm"),
                                                                                direction = "horizontal"),
                                                         feature2 = list(title = "immune",
                                                                   labels_gp = gpar(fontsize = 18,fontface = "bold"),
                                                                   title_gp = gpar(fontsize = 20,fontface = "bold"),legend_width = unit(5, "cm"),direction = "horizontal")))
  
  ht1 = Heatmap(as.matrix(feature1Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_1)
  
  row_labels_ht3 = gsub("|\\d+$", "", rownames(feature3Pathway)) 
  ht3 = Heatmap(as.matrix(feature3Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_3)
  row_labels_ht4 = gsub("|\\d+$", "", rownames(feature4Pathway)) 
  ht4 = Heatmap(as.matrix(feature4Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_4)
  row_labels_ht5 = gsub("|\\d+$", "", rownames(feature5Pathway)) 
  ht5 = Heatmap(as.matrix(feature5Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_5)
  row_labels_ht6 = gsub("|\\d+$", "", rownames(feature6Pathway)) 
  ht6 = Heatmap(as.matrix(feature6Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_6)
  row_labels_ht7 = gsub("|\\d+$", "", rownames(feature7Pathway)) 
  ht7 = Heatmap(as.matrix(feature7Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_7)
  row_labels_ht8 = gsub("|\\d+$", "", rownames(feature8Pathway)) 
  ht8 = Heatmap(as.matrix(feature8Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_8)
  row_labels_ht9 = gsub("|\\d+$", "", rownames(feature9Pathway)) 
  ht9 = Heatmap(as.matrix(feature9Pathway),
                col = colorGsva,
                show_row_names = FALSE,
                width = unit(29, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE,show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                show_column_dend = FALSE,top_annotation = ha_9)
  row_labels_ht10 = gsub("|\\d+$", "", rownames(feature10Pathway)) 
  ht10 = Heatmap(as.matrix(feature10Pathway),
                 col = colorGsva,
                 show_row_names = FALSE,
                 width = unit(29, "cm"), height = unit(2, "cm"),
                 cluster_rows = FALSE,show_column_names = FALSE,
                 show_heatmap_legend = FALSE,
                 show_column_dend = FALSE,top_annotation = ha_10)
  row_labels_ht11 = gsub("|\\d+$", "", rownames(feature11Pathway)) 
  ht11 = Heatmap(as.matrix(feature11Pathway),
                 name = " Pathway GSVA",
                 col = colorGsva,
                 show_row_names = FALSE,
                 heatmap_legend_param = list(legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 18,fontface = "bold"),title_gp = gpar(fontsize = 20,fontface = "bold"),title_gap = unit(25, "cm"),direction = "horizontal"),
                 width = unit(29, "cm"), height = unit(2, "cm"),
                 cluster_rows = FALSE,show_column_names = FALSE,
                 # show_heatmap_legend = FALSE,
                 show_column_dend = FALSE,top_annotation = ha_11)
  row_labels_ht12 = gsub("|\\d+$", "", rownames(feature12Pathway)) 
  ht12 = Heatmap(as.matrix(feature12Pathway),
                 col = colorGsva,
                 show_row_names = FALSE,
                 width = unit(29, "cm"), height = unit(2, "cm"),
                 cluster_rows = FALSE,show_column_names = FALSE,
                 show_heatmap_legend = FALSE,
                 show_column_dend = FALSE,top_annotation = ha_12)
  row_labels_ht13 = gsub("|\\d+$", "", rownames(feature13Pathway)) 
  ht13 = Heatmap(as.matrix(feature13Pathway),
                 col = colorGsva,
                 show_row_names = FALSE,
                 width = unit(29, "cm"), height = unit(2, "cm"),
                 cluster_rows = FALSE,show_column_names = FALSE,
                 show_heatmap_legend = FALSE,
                 show_column_dend = FALSE,top_annotation = ha_13)
  
  ht2 = Heatmap(as.matrix(feature2Pathway), name = "GSVA",
                column_title = "Radiogenomics Discovery Cohort                                                  Radiogenomics Validation Cohort",
                column_title_gp = gpar(
                  fontsize = 24,fontface = "bold", height=unit(0.01, "mm")),
                col = colorGsva,
                show_row_names = FALSE,
                heatmap_legend_param = list(direction = "horizontal"),
                width = unit(50, "cm"), height = unit(2, "cm"),
                cluster_rows = FALSE, show_column_names = FALSE,
                column_split = factor(group, levels = c("A", "B","C", "D")),
                column_gap =  unit(c(2, 10, 2), "mm"),
                show_heatmap_legend = FALSE,
                cluster_columns = FALSE,top_annotation = ha_2)
  
  ht_list_g2 = ht2 %v% ht3 %v% ht5 %v% ht12 %v% ht13 %v% ht1 %v% ht9 %v% ht4 %v% ht6 %v% ht7 %v% ht8 %v% ht10 %v% ht11 
  
  lgd_list = list(Legend(labels = c("inflammatory","proliferative","DNA-damage responsive"),legend_gp =  gpar(fill = c( "palevioletred3", "goldenrod3","palegreen4")), title = "type",
                         direction = "horizontal"))
  setwd(resultDir)
  tiff(file = paste('./','figureHeatmapGSVA.tiff',sep=''), res = 600, width = 15500, height = 12000, compression = "lzw")
  par(mar=c(0,0,0,0))
  draw(ht_list_g2,ht_gap =  unit(c(rep(2,4),8,2,8,2,8,rep(2,3)), "mm"),
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  decorate_title("GSVA", {
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    x1 = c(first_index(group == "A")-1, first_index(group == "C")+0.5) 
    x2 = c(last_index(group == "B")-2, last_index(group == "D"))
    grid.rect(x = x1/length(group), width = (x2 - x1)/length(group), just = "left",
              default.units = "npc", gp = gpar(fill = c("Grey", "Grey"), col = NA))
    grid.text("Radiogenomics Discovery Cohort                                                  Radiogenomics Validation Cohort",gp = gpar(fontsize = 24,fontface = "bold"))
  })
  
  annotation_titles = c(risk_group = "risk group",
                        RadRisk_score = "RadRisk score",
                        overall_survival = "overall survival(days)",
                        survival_status = "survival status")
  for(an in names(annotation_titles)) {
    decorate_annotation(an, {
      grid.text(annotation_titles[an], unit(-0.02,'npc'), just = "right",gp = gpar(fontsize = 18,fontface = "bold"))},
      envir = as.environment(1L)
    )
  }
  
  annotation_featues = c(feature1 = expression(italic(f[1])),
                         feature2 = expression(italic(f[2])),
                         feature3 = expression(italic(f[3])),
                         feature4 = expression(italic(f[4])),
                         feature5 = expression(italic(f[5])),
                         feature6 = expression(italic(f[6])),
                         feature7 = expression(italic(f[7])),
                         feature8 = expression(italic(f[8])),
                         feature9 = expression(italic(f[9])),
                         feature10 = expression(italic(f[10])),
                         feature11 = expression(italic(f[11])),
                         feature12 = expression(italic(f[12])),
                         feature13 = expression(italic(f[13])))
  for(an in names(annotation_featues)) {
    decorate_annotation(an, {
      grid.text(annotation_featues[an], y = unit(0, "npc") ,unit(3.29,'npc'), just = "left",gp = gpar(fontsize = 32,fontface = "bold"))},
      envir = as.environment(1L)
    )
  }
  dev.off()
  }


# ====================function11 Heatmap of the z-score normalized expression of 30 key genes ====================
#  Heatmap of the z-score normalized expression of 30 key genes across RDC and RVC patients 
# ========input========
# radioGenTrain_heatmap3_key_gene,externalTest_heatmap3_key_gene,resultDir
# ========output========
# figureKeyGene.tiff

heatmapKeyGene <- function(radioGenTrain_heatmap3_key_gene,externalTest_heatmap3_key_gene,resultDir){
  data = radioGenTrain_heatmap3_key_gene
  radSig <- as.numeric(data[,1])
  survTime <- as.numeric(data[,2])
  survStatus <- as.numeric(data[,3])
  geneMean = as.numeric(data[,4])
  gapsCol = 53
  group = c(rep("1",gapsCol),rep("0",nrow(data)-gapsCol))
  data = data[,5:ncol(data)]
  dataRGT = as.data.frame(t(data))
  colnames(dataRGT) = seq(1:93)
  
  rowSpiltGene = c(rep("red module",20),rep("blue module",10))
  colorGene <- colorRampPalette(c("navy", "white", "firebrick3"))(256)
  colorSig <- colorRamp2(c(max(radSig),median(radSig),min(radSig)),rev(brewer.pal(9, "Blues"))[c(1,3,6)])
  colorOs <- colorRamp2(c(max(survTime),median(survTime),min(survTime)),usecol(pal = rev(pal_petrol),3))
  colorGeneMean = colorRamp2(c(max(geneMean),median(geneMean),min(geneMean)),sequential_hcl(3,"BluYl"))
  
  ha_0 = HeatmapAnnotation(group_1 = group,
                           sig_1 = radSig,
                           os_1 = survTime,
                           status_1 = survStatus,
                           gene = geneMean,
                           border =  TRUE,
                           # gp = gpar(col = "grey"),
                           col = list(group_1 = c("1" = "red", "0" = "forestgreen"),
                                      sig_1 = colorSig,
                                      os_1 = colorOs,
                                      status_1 =  c("1" = "white", "0" = "dimgrey"),
                                      gene = colorGeneMean),
                           gap = unit(c(rep(0.5,3),2), "mm"),
                           show_annotation_name = FALSE,
                           annotation_height =   unit(c(1, 1, 1, 1, 0.1), c("mm", "mm","mm","mm", "mm")),
                           annotation_legend_param = list(group_1 = list(at = c(0, 1),
                                                                         labels = c("low", "high"),
                                                                         title = "risk group",
                                                                         gp = gpar(fontsize = 16,fontface = "bold"),
                                                                         direction = "horizontal"),
                                                          sig_1 = list(title = "RadRisk score",
                                                                       gp = gpar(fontsize = 16,fontface = "bold"),
                                                                       legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          os_1 = list(title = "overall survival",
                                                                      gp = gpar(fontsize = 16,fontface = "bold"),
                                                                      legend_width = unit(4, "cm"),direction = "horizontal"),
                                                          status_1 = list(at = c(0, 1),
                                                                          labels = c("censored", "uncensored"),
                                                                          title = "survival status",gp = gpar(fontsize = 16,fontface = "bold"),
                                                                          border = TRUE,
                                                                          direction = "horizontal"),
                                                          gene = list(title = "RadGene score",
                                                                      gp = gpar(fontsize = 16,fontface = "bold"),
                                                                      legend_width = unit(4, "cm"),direction = "horizontal")))
  
  heat3_RGT = Heatmap(as.matrix(dataRGT), name = "gene expression",
                    column_title = "Radiogenomics Discovery Cohort",
                    column_title_gp = gpar(fill = "grey", col = "black", border = "grey", 
                                           fontsize = 12,fontface = "bold", height=unit(0.01, "mm")),
                    col = colorGene,
                    show_row_names = FALSE,
                    # row_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm")),
                    width = unit(16, "cm"), height = unit(14, "cm"),
                    cluster_rows = FALSE,
                    show_column_names = FALSE,
                    column_split = factor(group, levels = c("1", "0")),
                    row_split = factor(rowSpiltGene, levels = c("red module", "blue module")),
                    show_row_dend = FALSE,
                    column_gap = unit(4, "mm"),
                    row_gap = unit(2, "mm"),
                    cluster_columns = FALSE,top_annotation = ha_0,left_annotation = rowAnnotation(module = anno_block(gp = gpar(fill = c("red", "blue")))))
  
  
  data = externalTest_heatmap3_key_gene
  radSig <- as.numeric(data[,1])
  survTime <- as.numeric(data[,2])
  survStatus <- as.numeric(data[,3])
  geneMean = as.numeric(data[,4])
  gapsCol = 39
  group = c(rep("1",gapsCol),rep("0",nrow(data)-gapsCol))
  data = data[,5:ncol(data)]
  dataET = as.data.frame(t(data))
  colnames(dataET) = seq(1:76)
  
  colorGeneMean = colorRamp2(c(max(geneMean),median(geneMean),min(geneMean)),sequential_hcl(3,"BluYl"))
  
  colorSig <- colorRamp2(c(max(radSig),median(radSig),min(radSig)),rev(brewer.pal(9, "Blues"))[c(1,3,6)])
  colorOs <- colorRamp2(c(max(survTime),median(survTime),min(survTime)),usecol(pal = rev(pal_petrol),3))
  
  ha_1 = HeatmapAnnotation(group = group,
                           sig = radSig,
                           os = survTime,
                           status = survStatus,
                           gene = geneMean,
                           border =  TRUE,
                           # gp = gpar(col = "grey"),
                           col = list(group = c("1" = "red", "0" = "forestgreen"),
                                      sig = col_sig,
                                      os = col_os,
                                      status =  c("1" = "white", "0" = "dimgrey"),
                                      gene = colorGeneMean),
                           gap = unit(c(rep(0.5,3),2), "mm"),
                           show_annotation_name = FALSE,
                           show_legend = FALSE,
                           annotation_height =   unit(c(1, 1, 1, 1, 0.1), c("mm", "mm","mm","mm", "mm")),
                           annotation_legend_param = list(group = list(at = c(0, 1),
                                                                       labels = c("low", "high"),
                                                                       title = "risk group",
                                                                       direction = "horizontal"),
                                                          sig = list(title = "RadRisk score",direction = "horizontal"),
                                                          os = list(title = "OS",direction = "horizontal"),
                                                          status = list(at = c(0, 1),
                                                                        labels = c("censored", "uncensored"),
                                                                        title = "survival status",
                                                                        border = TRUE,
                                                                        direction = "horizontal"),
                                                          gene = list(title = "RadGene score",
                                                                      gp = gpar(fontsize = 16,fontface = "bold"),
                                                                      legend_width = unit(4, "cm"),direction = "horizontal")))
  
  heat3_ET = Heatmap(as.matrix(dataET), name = "gene ",
                    column_title = "Radiogenomics Validation Cohort",
                    column_title_gp = gpar(fill = "grey", col = "black", border = "grey", 
                                           fontsize = 12,fontface = "bold", height=unit(0.01, "mm")),
                    col = colorGene,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 12,fontface = "bold"),
                    show_heatmap_legend = FALSE,
                    # row_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(direction = "horizontal"),
                    width = unit(13, "cm"), height = unit(14, "cm"),
                    cluster_rows = FALSE, show_column_names = FALSE,
                    column_split = factor(group, levels = c("1", "0")),
                    row_split = factor(rowSpiltGene, levels = c("red module", "blue module")),
                    column_gap = unit(4, "mm"),
                    row_gap = unit(2, "mm"),
                    cluster_columns = FALSE,top_annotation = ha_1)
  
  heat3_list = heat3_RGT + heat3_ET
  
  
  
  
  setwd(resultDir)
  tiff(file = paste('./','figureKeyGene.tiff',sep=''), res = 600, width = 12000, height = 6000, compression = "lzw")
  # png(file = paste('./','Heatmap3.png.',sep=''), width = 1000, height = 2000)
  par(mar=c(0,0,0,0))
  # draw(ht1)
  draw(heat3_list,merge_legend = TRUE,ht_gap = unit(3, "mm"),
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  annotation_titles = c(group = "risk group",
                        sig = "RadRisk score",
                        os = "overall survival",
                        status = "survival status",
                        gene = "RadGene score")
  
  for(an in names(annotation_titles)) {
    decorate_annotation(an, {
      grid.text(annotation_titles[an], unit(2.03,'npc'), just = "left",gp = gpar(fontsize = 12,fontface = "bold"))
    })
  }
  dev.off()
}

myNormalize <- function (target) {
  # 2*(target - min(target))/(max(target) - min(target))-1
  scale(target)}


