# import functions from main directory
source("~/PDX/pdx.ql1.smooth.eval.R")
source("~/PDX/pdx.ql2.smooth.eval.R")
source("~/PDX/pdx.ql.rf1.smooth.eval.R")
source("~/PDX/pdx.ql.rf2.smooth.eval.R")

library(caret)
library(glmnet)
library(stringi)

pdx.super.learning.4 = function(cancer.type, outcome, gene.data.file, numgenes, input_dir = '/mnt/home/yzhan234/PDX/PDX.CodeAndData/super-learning', output_dir = '/mnt/home/yzhan234/PDX/PDX.CodeAndData/super-learning', k = 5, seed = 1, strip = T, outstring = "_sl4.csv"){
  
  O1 = pdx.ql1.smooth.eval(cancer.type, outcome, gene.data.file, numgenes, input_dir, output_dir)
  O2 = pdx.ql2.smooth.eval(cancer.type, outcome, gene.data.file, numgenes, input_dir, output_dir)
  O3 = pdx.ql.rf1.smooth.eval(cancer.type, outcome, gene.data.file, numgenes, input_dir, output_dir)
  O4 = pdx.ql.rf2.smooth.eval(cancer.type, outcome, gene.data.file, numgenes, input_dir, output_dir)
  
  O.list = list(O1, O2, O3, O4)
  J = length(O.list)
  
  setwd(input_dir)
  load("split.cm.data.rda")
  
  load("trts.by.cancer.rda")
  
  # random forest predicted values -- use these as outcomes when estimating decision rule
  load("pred_vals_rf.rda")
  
  # extract clinical data for given cancer type
  if (cancer.type == "BRCA") {dat = split.cm.data$BRCA; clinical = dat[ , 1:17]}
  if (cancer.type == "CM") {dat = split.cm.data$CM;clinical = dat[ , 1:17]}
  if (cancer.type == "CRC") {dat = split.cm.data$CRC; clinical = dat[ , 1:17]}
  if (cancer.type == "NSCLC") {dat = split.cm.data$NSCLC; clinical = dat[ , 1:17]}
  if (cancer.type == "PDAC") {dat = split.cm.data$PDAC; clinical = dat[ , 1:17]}
  if (cancer.type == "overall") {
    load("full.data.rda")
    clinical = dat[ , 1:5]
  }
  clinical$RespScaled = -clinical$RespScaled  # reverse sign of best average response -- this way larger values are better
  
  biomarkers = read.csv(gene.data.file)
  if (strip) biomarkers = biomarkers[ , -1]
  
  # remove duplicated columns from biomarkers
  if (sum(duplicated(as.matrix(biomarkers), MARGIN = 2)) > 0) biomarkers = biomarkers[ , -which(duplicated(as.matrix(biomarkers), MARGIN = 2))]
  
  # center/scale biomarkers
  center.scale = function(x) return((x - mean(x)) / sd(x))
  biomarker.temp = apply(biomarkers, 2, center.scale)
  biomarkers = biomarker.temp
  
  # format data
  biomarkers = biomarkers[!duplicated(biomarkers), ]  # here biomarkers contains one observation per line
  
  new.resp.mat = matrix(NA, nrow = length(unique(clinical$Model)), ncol = length(unique(clinical$Treatment)))
  new.surv.mat = matrix(NA, nrow = length(unique(clinical$Model)), ncol = length(unique(clinical$Treatment)))
  rownames(new.resp.mat) = unique(clinical$Model); colnames(new.resp.mat) = unique(clinical$Treatment)
  rownames(new.surv.mat) = unique(clinical$Model); colnames(new.surv.mat) = unique(clinical$Treatment)
  for (dim1 in 1:nrow(new.resp.mat)) {
    for (dim2 in 1:ncol(new.resp.mat)) {
      row = which(clinical$Model == rownames(new.resp.mat)[dim1] & clinical$Treatment == colnames(new.resp.mat)[dim2])
      if (length(row) != 0) {
        new.resp.mat[dim1, dim2] = clinical$RespScaled[row]
        new.surv.mat[dim1, dim2] = clinical$logSurvScaled[row]
      }
    }
  }  # end format of clinical data
  
  # clinical.other contains outcomes not specified for estimation
  if (outcome == "BAR") {
    clinical = new.resp.mat 
    clinical.other = new.surv.mat
  }
  if (outcome == "Surv") {
    clinical = new.surv.mat 
    clinical.other = new.resp.mat
  } 
  
  # create folds for super-learning cross-validation
  set.seed(seed)
  folds.sl = createFolds(1:dim(biomarkers)[1], k = k, list = TRUE, returnTrain = FALSE)
  
  test.clin.list = list()
  test.clin.other.list = list()
  
  for (f.sl in 1:length(folds.sl)) {
    # select testing sets
    test.clin = clinical[folds.sl[[f.sl]], ]
    # replace na by the mean untreated value
    for(i in 1:nrow(test.clin)){
      if (is.na(test.clin[i,"untreated"]))
        test.clin[i,"untreated"] = mean(test.clin[ ,"untreated"], na.rm = T)
    }
    test.clin = test.clin - test.clin[,"untreated"]
    test.clin.list[[f.sl]] = test.clin
    
    # select testing sets
    test.clin.other = clinical.other[folds.sl[[f.sl]], ]
    # replace na by the mean untreated value
    for(i in 1:nrow(test.clin.other)){
      if (is.na(test.clin.other[i,"untreated"]))
        test.clin.other[i,"untreated"] = mean(test.clin.other[ ,"untreated"], na.rm = T)
    }
    test.clin.other = test.clin.other - test.clin.other[,"untreated"]
    test.clin.other.list[[f.sl]] = test.clin.other    
  }  
  
  
  # function to generate simplex of length n
  get.simplex = function(n){
    alpha = rep(0, n)
    sum = 0
    for(i in 1:(n-1)){
      alpha[i] = runif(1, 0, 1 - sum)
      sum = sum + alpha[i]
    }
    alpha[n] = 1 - sum
    
    return (alpha)
  }
  
  obj.function = function(alpha, O.list, test.clin.list, k){
    J = length(O.list)
    
    value = 0
    
    # for each super-learning fold
    for(v in 1:k){
      # calculated the weighted average of latent functions for all methods
      O = O.list[[1]][[v]] * alpha[1]
      for(j in 2:J){
        O = O + O.list[[j]][[v]] * alpha[j]
      }
      
      # for each test model
      for(i in 1:nrow(test.clin.list[[v]])){
        # make decision according to O
        # [1] is for tie case
        value.decision = test.clin.list[[v]][i,which(O[i,] == max(O[i,]))][1]
        # if our decision is missing in the data, then we choose among non-NA values
        if(is.na(value.decision)){
          index.not.NA = which(is.na(test.clin.list[[v]][i,]) == 0)
          # intersect is for the case when many O[i,] == max, and the first one is NA
          value.decision = test.clin.list[[v]][i,][intersect(which(O[i,] == max(O[i,index.not.NA])), index.not.NA)][1]
        }
        
        value = value - value.decision
      }
    }
    
    return(value)
  }
  
  # randomly generate n.try alpha
  # and pick the best n.pick alpha
  n.try = 1000
  n.pick = 10
  
  alpha.random = matrix(0, n.try, J)
  f.values = rep(0, n.try)
  
  for(i in 1:n.try){
    alpha.random[i,] = get.simplex(J)
    f.values[i] = obj.function(alpha.random[i,], O.list, test.clin.list, k)
  }
  
  # pick the best n.pick alpha and randomly choose from the tie case
  alpha.start = alpha.random[which(f.values %in% sort(f.values)[1:n.pick]),][1:n.pick,]
  
  n.stage = 60
  temp = rep(5, n.stage)
  for(i in 2:n.stage){
    temp[i] = 0.9 * temp[i - 1]
  }
  stage.length = rep(c(60,120,220), each = 20)
  e = 0.5
  
  
  
  # record the final alpha and function value
  alpha.final = alpha.start
  f.final = rep(0, n.pick)
  
  #f.record = rep(0, sum(stage.length))
  
  for(start in 1:n.pick){
    alpha.old = alpha.start[start,]
    f.old = obj.function(alpha.old, O.list, test.clin.list, k)
    
    #step = 1
    
    for(stage in 1:n.stage){
      for(i in 1:stage.length[stage]){
        alpha.tmp = get.simplex(J)
        alpha.new = alpha.old * (1-e) + alpha.tmp * e
        f.new = obj.function(alpha.new, O.list, test.clin.list, k)
        
        if(runif(1) < exp((f.old - f.new)/temp[stage])){
          f.old = f.new
          alpha.old = alpha.new
        }
        
        #f.record[step] = f.old
        #step = step + 1
      }
    }
    alpha.final[start,] = alpha.old
    f.final[start] = f.old
  }
  
  f.best = min(f.final)
  if(length(which(f.final == f.best)) == 1){
    alpha.best = alpha.final[which(f.final == f.best),]
  }else{
    alpha.best = alpha.final[which(f.final == f.best),][1,]
  }  
  
  # reevaluate and get the covariance
  main.outs = NULL
  other.outs = NULL
  
  # for each super-learning fold
  for(v in 1:k){
    main.value = 0
    other.value = 0
    
    # calculated the weighted average of latent functions for all methods
    O = O.list[[1]][[v]] * alpha.best[1]
    for(j in 2:J){
      O = O + O.list[[j]][[v]] * alpha.best[j]
    }
    
    # for each test model
    for(i in 1:nrow(test.clin.list[[v]])){
      # make decision according to 
      # [1] is for tie case
      main.value.decision = test.clin.list[[v]][i,which(O[i,] == max(O[i,]))][1]
      # if our decision is missing in the data, then we choose among non-NA values
      if(is.na(main.value.decision)){
        index.not.NA = which(is.na(test.clin.list[[v]][i,]) == 0)
        main.value.decision = test.clin.list[[v]][i,][intersect(which(O[i,] == max(O[i,index.not.NA])), index.not.NA)][1]
      }
      
      other.value.decision = test.clin.other.list[[v]][i,which(O[i,] == max(O[i,]))][1]
      # if our decision is missing in the data, then we choose among non-NA values
      if(is.na(other.value.decision)){
        index.not.NA = which(is.na(test.clin.other.list[[v]][i,]) == 0)
        other.value.decision = test.clin.other.list[[v]][i,][intersect(which(O[i,] == max(O[i,index.not.NA])), index.not.NA)][1]
      }      
      
      main.value = main.value + main.value.decision
      other.value = other.value + other.value.decision
    }
    
    main.value = main.value / nrow(test.clin.list[[v]])
    other.value = other.value / nrow(test.clin.list[[v]])
    
    main.outs = c(main.outs, main.value)
    other.outs = c(other.outs, other.value)
  }
  
  if (outcome == "BAR") {
    final.resp = mean(main.outs)
    final.surv = mean(other.outs)
  }
  if (outcome == "Surv") {
    final.resp = mean(other.outs)
    final.surv = mean(main.outs)
  }
  
  cov.mat = cov(matrix(c(main.outs, other.outs), ncol = 2), use = "complete.obs")
  
  # observed mean for primary outcome
  for(i in 1:nrow(clinical)){
    if (is.na(clinical[i,"untreated"]))
      clinical[i,"untreated"] = mean(clinical[ ,"untreated"], na.rm = T)
  }
  clinical = clinical - clinical[,"untreated"]  # subtract off mean of untreated group in each line
  # observed mean for secondary outcome
  for(i in 1:nrow(clinical.other)){
    if (is.na(clinical.other[i,"untreated"]))
      clinical.other[i,"untreated"] = mean(clinical.other[ ,"untreated"], na.rm = T)
  }
  clinical.other = clinical.other - clinical.other[,"untreated"]  # subtract off mean of untreated group in each line
  
  if (outcome == "BAR") {
    obs.resp = mean(clinical, na.rm = T)
    obs.surv = mean(clinical.other, na.rm = T)
    opt.resp = mean(apply(clinical, 1, max, na.rm = T)[!is.infinite(apply(clinical, 1, max, na.rm = T))], na.rm = T)
    opt.surv = mean(apply(clinical.other, 1, max, na.rm = T)[!is.infinite(apply(clinical.other, 1, max, na.rm = T))], na.rm = T)
    var.resp = cov.mat[1, 1]
    var.surv = cov.mat[2, 2]
    cov = cov.mat[1, 2]
  }
  if (outcome == "Surv") {
    obs.resp = mean(clinical.other, na.rm = T)
    obs.surv = mean(clinical, na.rm = T)
    opt.resp = mean(apply(clinical.other, 1, max, na.rm = T)[!is.infinite(apply(clinical.other, 1, max, na.rm = T))], na.rm = T)
    opt.surv = mean(apply(clinical, 1, max, na.rm = T)[!is.infinite(apply(clinical, 1, max, na.rm = T))], na.rm = T)
    var.resp = cov.mat[2, 2]
    var.surv = cov.mat[1, 1]
    cov = cov.mat[1, 2]     
  }
  
  res = data.frame(c1 = "/", c2 = "/", mean.response = final.resp, mean.survival = final.surv, 
                   var.response = var.resp, var.survival = var.surv, covariance = cov, observed.resp = obs.resp, 
                   observed.surv = obs.surv, optimal.resp = opt.resp, optimal.surv = opt.surv)
  
  rownames(res) = NULL
  setwd(output_dir)
  output.name = paste(cancer.type, "_", outcome, "_", stri_sub(gene.data.file, 1, -5), outstring, sep = "")
  write.csv(res, output.name)
  
}

