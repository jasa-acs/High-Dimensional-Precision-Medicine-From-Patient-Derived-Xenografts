# pdx.bdc1.R

# prognostic gene screening

library(stringi)
library(energy)
library(CCP)

# function to screen prognostic biomarkers
# takes as imput a tumor type and a treatment type
# input_dir is the directory where data files are (main directory)
# output_dir is the directory where output files should be sent (scratch directory)
pdx.bdc1 = function(tumor_type, trt_type, input_dir, output_dir) {
  
  setwd(input_dir)
  load("split.cm.data.rda")
  load("trts.by.cancer.rda")
  load("genes.by.cancer.rda")
  
  if (tumor_type == "overall") {
    source("create.full.data.R")
    dat = create.full.data(split.cm.data, trts.by.cancer, genes.by.cancer)
    cm.data = dat[ , 1:5]
    gene.data = dat[ , 7:ncol(dat)]
  }
  
  # extract gene data for given tumor type
  if (tumor_type == "BRCA") {
    gene.data = split.cm.data$BRCA[ , 18:ncol(split.cm.data$BRCA)]
    cm.data = split.cm.data$BRCA[ , 1:17]
  }
  if (tumor_type == "PDAC") { 
    gene.data = split.cm.data$PDAC[ , 18:ncol(split.cm.data$PDAC)]
    cm.data = split.cm.data$PDAC[ , 1:17]
  }
  if (tumor_type == "NSCLC") {
    gene.data = split.cm.data$NSCLC[ , 18:ncol(split.cm.data$NSCLC)]
    cm.data = split.cm.data$NSCLC[ , 1:17]
  }
  if (tumor_type == "CRC") {
    gene.data = split.cm.data$CRC[ , 18:ncol(split.cm.data$CRC)]
    cm.data = split.cm.data$CRC[ , 1:17]
  }
  if (tumor_type == "CM") {
    gene.data =split.cm.data$CM[ , 18:ncol(split.cm.data$CM)]
    cm.data = split.cm.data$CM[ , 1:17]
  }
  
  gene.data = unique(gene.data)
  rownames(gene.data) = unique(cm.data$Model)
  
  gene.data = as.data.frame(gene.data)
  
  # remove variables that aren't needed here
  gene.data = gene.data[ , !(names(gene.data) %in% c("ArmLevelCNScore.cn", "FocalCNScore.cn"))]

  # list of gene names
  split.fun = function(x) return(strsplit(x, "\\.")[[1]][1])
  temp.names = matrix(colnames(gene.data))
  gene.names = apply(temp.names, 1, split.fun)

  # clinical data for selected trt only
  keep.trt = c(trt_type)
  cm.data = cm.data[which(cm.data$Treatment %in% keep.trt), ]

  # create two dimensional response variable
  Y = cbind(cm.data$BestAvgResponse, cm.data$TimeToDouble)
  if (tumor_type == "overall") Y = cbind(cm.data$RespScaled, cm.data$logSurvScaled)
  n = dim(Y)[1]
  tab = table(cm.data$Model)

  # screen all genes in the same file
  start_gene = 1; end_gene = length(gene.names)

  # to store the necessary results
  genes = matrix(NA, nrow = end_gene - start_gene + 1, ncol = 1)  # gene names
  cors = matrix(NA, nrow = end_gene - start_gene + 1, ncol = 1)  # BDC
  pvals = matrix(NA, nrow = end_gene - start_gene + 1, ncol = 1)  # BDC p-values
  ccps = matrix(NA, nrow = end_gene - start_gene + 1, ncol = 1)  # canonical correlation p-values
  cca.stat = matrix(NA, nrow = end_gene - start_gene + 1, ncol = 1)  # canonical correlation statistic

  # loop through genes
  for (i in start_gene:end_gene) {
  
    # get variable names for each platform for current gene
    cur.gene = gene.names[i]
    rna.var.name = paste(gene.names[i], ".rna", sep = "")
    #cn.var.name = paste(gene.names[i], ".cn", sep = "")
    #mut.var.name = paste(gene.names[i], ".mut", sep = "")
    #var.names = c(rna.var.name, cn.var.name, mut.var.name)
    var.names = c(rna.var.name)
  
    # gene data for current gene
    X.temp = as.matrix(gene.data[ , which(colnames(gene.data) %in% var.names)])
    rownames(X.temp) = rownames(gene.data)
    X.temp = as.matrix(X.temp[which(rownames(X.temp) %in% unique(cm.data$Model)), ])

    # when there is only one mouse with a given tumor/trt combination, we apply this to correct for the way as.matrix works
    if (dim(cm.data)[1] == 1) X.temp = t(X.temp)
    X = NULL
    for (j in 1:dim(X.temp)[1]) {
    
      temp = rep(X.temp[j, ], tab[j])
      temp = matrix(temp, nrow = tab[j], byrow = T)
    
      X = rbind(X, temp)
    
    }

    # calculate distance correlation for genes not entirely screened out
    # only runs if there is at least one platform left for current gene after unsupervised screening
    if (dim(X)[2] > 0) {
  
      # Brownian distance correlation
      gene.cor = dcor.ttest(X, Y)
      genes[i - start_gene + 1] = cur.gene
      cors[i - start_gene + 1] = gene.cor$estimate
      pvals[i - start_gene + 1] = gene.cor$p.value
    
      # canonical correlation analysis
      try(expr = (gene.ccp = p.perm(X, Y)), silent = T)
      try(expr = (ccps[i - start_gene + 1] = gene.ccp$p.value), silent = T)
      try(expr = (cca.stat[i - start_gene + 1] = max(cancor(X, Y)$cor)), silent = T)
      try(expr = (rm(gene.ccp)), silent = T)
    
    }
    
  }

  # make a data frame of results
  dcors = data.frame("gene.name" = as.character(genes), "correlation" = cors, "bdc.p.value" = pvals, "cca.p.value" = ccps, "cca.stat" = cca.stat)

  # remove gene names that are NA
  dcors = dcors[which(!is.na(dcors$gene.name)), ]

  # sort by distance correlation p-value
  dcors = dcors[order(dcors$bdc.p.value), ]

  # file name to write results to
  outfile = paste(tumor_type, "_", trt_type, ".screen.csv", sep = "") 

  # output results
  setwd(output_dir)
  write.csv(x = dcors, file = outfile)

}
