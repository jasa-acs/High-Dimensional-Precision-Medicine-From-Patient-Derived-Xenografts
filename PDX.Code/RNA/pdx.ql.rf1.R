# pdx.ql.rf1.R

# Q-learning on PDX data without pseudo-values

library(caret)
library(randomForest)
library(stringi)

# this function performs Q-learning on PDX data
# need to supply to the function a cancer type, one of 
# "BRCA", "CM", "CRC", "NSCLC", or "PDAC"
# and an outcome, one of "BAR" for best average response 
# or "Surv" for survival
# gene.data.file is the name of a csv file with gene data
# input_dir and output_dir are directories (both should be main directory)
# k is the number of folds for cross-validation
# also need to supply a seed to make training/testing sets for cross-validation 
# if strip is TRUE, we assume the first column of gene data contains row names and is dropped
pdx.ql.rf1 = function(cancer.type, outcome, gene.data.file, input_dir, output_dir, k = 5, seed = 1, strip = T, outstring = "_qlrf1_rna.csv") {
  
  setwd(input_dir)
  load("split.cm.data.rda")
  
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

  # for across-cancer analysis, add variable for cancer type
  if (cancer.type == "overall") {
    biomarkers = as.data.frame(biomarkers)
    biomarkers$brca = as.numeric(clinical$TumorType == "BRCA")
    biomarkers$cm = as.numeric(clinical$TumorType == "CM")
    biomarkers$crc = as.numeric(clinical$TumorType == "CRC")
    biomarkers$nsclc = as.numeric(clinical$TumorType == "NSCLC")
    biomarkers$pdac = as.numeric(clinical$TumorType == "PDAC")
    biomarkers = as.matrix(biomarkers)
  }
  
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
  
  # create folds for cross-validation
  set.seed(seed)
  folds = createFolds(1:dim(biomarkers)[1], k = k, list = TRUE, returnTrain = FALSE)
  
  # c1 is the number of trt's to group with untreated, c2 is the number of steps to take down the tree
  c1s = 0:(ncol(clinical) - 3); c1s = 0:4
  if (cancer.type == "overall") c1s = c(0)
  avg.main.outs = NULL  # store primary value functions across c1 and c2
  avg.other.outs = NULL  # store secondary value functions across c1 and c2
  parameters = matrix(NA, nrow = length(c1s) * ncol(clinical), ncol = 2)  # store pairs of c1 and c2
  colnames(parameters) = c("c1", "c2")
  cov.list = list()  # to store covariances of value functions
  ctr = 1  # count number of times through inner loop
  
  # loop through parameters
  for (c1 in c1s) {
    
    c2s = seq(1:(ncol(clinical) - c1 - 2))
    
    for (c2 in c2s) {  
      
      # store primary and secondary value functions across folds
      main.folds = NULL
      other.folds = NULL
      
      # loop through folds
      for (f in 1:length(folds)) {
        
        # select training and testing sets
        train.bio = biomarkers[-folds[[f]], ]
        train.clin = clinical[-folds[[f]], ]
        test.bio = biomarkers[folds[[f]], ]
        test.clin = clinical[folds[[f]], ]
        test.clin.other = clinical.other[folds[[f]], ]
        
        # find the c1 closest treatments to "untreated"
        dist_mat = as.matrix(dist(t(train.clin)))
        col_ind = which(colnames(dist_mat) == "untreated")
        ordered_dist_mat = dist_mat[order(dist_mat[ , col_ind]) , col_ind]
        untrt = names(ordered_dist_mat[1:(1 + c1)])
        
        # average outcomes aross "No treatment group"
        untrt.ind = which(colnames(train.clin) %in% untrt)
        means = apply(as.matrix(train.clin[ , untrt.ind]), 1, mean, na.rm = T)
        train.clin = train.clin - means  # subtract off mean of untreated group in each line
        train.clin = train.clin[ , -untrt.ind]
        
        # same for the test set
        # replace na by the mean untreated value
        for(i in 1:nrow(test.clin)){
          if (is.na(test.clin[i,"untreated"]))
            test.clin[i,"untreated"] = mean(test.clin[ ,"untreated"], na.rm = T)
        }
        test.clin = test.clin - test.clin[,"untreated"]
        untrt.ind = which(colnames(test.clin) %in% untrt)
        test.clin = test.clin[ , -untrt.ind]
        
        # same for the test set of the secondary outcome
        # replace na by the mean untreated value
        for(i in 1:nrow(test.clin.other)){
          if (is.na(test.clin.other[i,"untreated"]))
            test.clin.other[i,"untreated"] = mean(test.clin.other[ ,"untreated"], na.rm = T)
        }
        test.clin.other = test.clin.other - test.clin.other[,"untreated"]
        untrt.ind = which(colnames(test.clin.other) %in% untrt)
        test.clin.other = test.clin.other[ , -untrt.ind]
        
        # create treatment tree by clustering
        clusters = hclust(dist(t(train.clin))) # repeat distance matrix after removing untreated columns
        
        # full mouse level data 
        rownames(train.bio) = rownames(train.clin)
        full.bio = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = ncol(train.clin)), ncol = 1), ]
        full.clin = matrix(as.numeric(t(train.clin)), ncol = 1)
        rownames(full.clin) = rownames(full.bio)  # we don't need full.bio after this...
        full.trt = matrix(rep(colnames(train.clin), nrow(train.clin)))
        full.clin = cbind(full.clin, full.trt)
        
        # create treatment variables for each step of the tree
        num.steps = dim(clusters$merge)[1]
        merge.steps = clusters$merge
        trt.list = clusters$labels
        new.trt.vars = NULL  # new treatment variables
        for (j in 1:num.steps) {
          temp = rep(NA, dim(full.clin)[1])
          merge1 = merge.steps[j, 1]; merge2 = merge.steps[j, 2]
          if (merge1 < 0) {
            temp[which(full.clin[ , 2] == trt.list[-merge1])] = 1
          }
          if (merge2 < 0) {
            temp[which(full.clin[ , 2] == trt.list[-merge2])] = -1
          }
          if (merge1 > 0) {
            temp[which(!is.na(new.trt.vars[ , merge1]))] = 1
          }
          if (merge2 > 0) {
            temp[which(!is.na(new.trt.vars[ , merge2]))] = -1
          }
          new.trt.vars = cbind(new.trt.vars, temp)
        }  # end creation of trt variables
        # change coding of treatments from 1/-1 to 1/0
        new.trt.vars[which(new.trt.vars == -1)] = 0
        # assign column names to new trt variables
        #new.trt.names = NULL
        #for (j in 1:dim(new.trt.vars)[2]) {
        #  name = paste("step_", j, sep = "")
        #  new.trt.names = c(new.trt.names, name)
        #}
        #colnames(new.trt.vars) = new.trt.names
        rownames(new.trt.vars) = rownames(full.clin)
        
        # select trt variables for c2 steps down tree
        trt.vars = new.trt.vars[ , rev(rev(1:dim(new.trt.vars)[2])[1:c2])]
        trt.vars = as.matrix(trt.vars)
        
        # Q-learning up the tree
        list.models = list()  # store models at each step
        list.pred = list()  # store predicted values at each step
        for (d in 1:dim(trt.vars)[2]) {
          
          X = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = 2), ncol = 1), ]
          
          Y = matrix(NA, nrow = nrow(X), ncol = 2)
          Y[ , 2] = rep(c(0, 1), nrow(train.bio))
          rownames(Y) = rownames(X)
          
          # for first step, outcomes are mean outcomes in root nodes
          if (d == 1) {
            for (g in 1:dim(Y)[1]) {
              Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(trt.vars) == rownames(Y)[g]), 1]), na.rm = T)
            }
          }
          
          # for other steps, outcomes are the maximum over outcomes on lower steps of the tree (if a previous decision has been made)
          if (d > 1) {
            for (g in 1:dim(Y)[1]) {
              
              # most recent previous step of tree where decision was made
              index = max(which(!is.na(as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2]), 1:(d - 1)])[1,])))
              ##print(index)
              if (is.infinite(index)) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
              if (!is.infinite(index)){
                while(!is.infinite(index)){
                  d.temp = index
                  choice = list.pred[[d.temp]][which(list.pred[[d.temp]][ , 1] == rownames(Y)[g]), 3]
                  if(d.temp == 1){
                    index = -Inf
                  }else{
                    index = max(which(!is.na(as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d.temp] == choice), 1:(d.temp - 1)])[1,])))
                  }
                }
                Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d.temp] == choice & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
              }
            }
          }
          
          A = Y[ , 2]
          Y = Y[ , 1]
          
          # fit model for Q-learning
          tempx = cbind(X, A)
          tempy = Y
          if (sum(is.na(Y) > 0)) tempy = tempy[-which(is.na(Y))]
          if (sum(is.na(Y) > 0)) tempx = tempx[-which(is.na(Y)), ]
          fit = randomForest(tempx, tempy)
          list.models[[d]] = fit
          
          # create pseudo-values for next step
          preds = matrix(NA, nrow = nrow(train.bio), ncol = 3)
          for (t in 1:dim(train.bio)[1]) {
            predx0 = c(train.bio[t, ], 0) 
            predx1 = c(train.bio[t, ], 1)
            preds[t, 1] = rownames(train.bio)[t]
            preds[t, 2] = max(predict(fit, matrix(predx0, nrow = 1)), predict(fit, matrix(predx1, nrow = 1)))
            preds[t, 3] = as.numeric(predict(fit, matrix(predx0, nrow = 1)) < predict(fit, matrix(predx1, nrow = 1)))
          }
          list.pred[[d]] = preds
        }  # end loop through steps in tree
        
        # test on validation set 
        main.outs = NULL  # save mean outcomes for each line among mice treated consistent with treatment rule
        other.outs = NULL
        rownames(test.bio) = rownames(test.clin)
        
        # loop through lines in test set
        for (t in 1:dim(test.bio)[1]) {
          
          # determine sequence of decisions for each line
          moves = NULL
          for (i in rev(1:dim(trt.vars)[2])) {
            predx0 = c(train.bio[t, ], 0) 
            predx1 = c(train.bio[t, ], 1)
            temp.fit = list.models[[i]]
            moves = c(moves, as.numeric(predict(temp.fit, matrix(predx0, nrow = 1)) < predict(temp.fit, matrix(predx1, nrow = 1))))
          }
          cur.step = merge.steps[nrow(merge.steps), ]
          trt.ind = NULL  # save indices (in trt.list) of those trts that are consistent with decision rule
          
          temp = rep(NA, nrow(merge.steps) - length(moves))
          tmoves = c(temp, rev(moves))
          
          # determine root node for each line in testing set
          keeplooping = T
          cur.move = tmoves[length(tmoves)]
          while (keeplooping) {
            if (cur.move == 1) next.step = cur.step[1]
            if (cur.move == 0) next.step = cur.step[2]  
            if (next.step < 0) {
              trt.ind = c(trt.ind, -next.step)
              keeplooping = F
            }
            if (next.step > 0) { 
              cur.step = merge.steps[next.step, ]
              cur.move = tmoves[next.step]
              if (is.na(cur.move)) keeplooping = F
            }
          }
          
          # recursive function to construct list of indices for trts in root node
          get.trt.list = function(cur.step) {
            trt.ind = NULL
            if (cur.step[1] < 0) trt.ind = c(trt.ind, -cur.step[1])
            if (cur.step[1] > 0) trt.ind = c(trt.ind, get.trt.list(merge.steps[cur.step[1], ]))
            if (cur.step[2] < 0) trt.ind = c(trt.ind, -cur.step[2])
            if (cur.step[2] > 0) trt.ind = c(trt.ind, get.trt.list(merge.steps[cur.step[2], ]))
            return(trt.ind)
          }  # end recursive function
          
          # if root node is not a single trt, get list of trts in root node
          if (length(trt.ind) == 0) trt.ind = get.trt.list(cur.step)
          
          # list of trts in root node
          treatments = trt.list[trt.ind]
          
          # determine if line should be treated at all
          predx0 = c(test.bio[t, ], 0) 
          predx1 = c(test.bio[t, ], 1)
          temp.fit = list.models[[dim(trt.vars)[2]]]
          to.treat = max(predict(temp.fit, matrix(predx0, nrow = 1)), predict(temp.fit, matrix(predx1, nrow = 1))) > 0
          
          # if predicted outcome from decision rule is larger, take mean outcome in root node
          if (to.treat) {
            main.outs = c(main.outs, mean(test.clin[which(rownames(test.clin) == rownames(test.bio)[t]), which(colnames(test.clin) %in% treatments)]))
            other.outs = c(other.outs, mean(test.clin.other[which(rownames(test.clin.other) == rownames(test.bio)[t]), which(colnames(test.clin.other) %in% treatments)]))
          }
          
          # if predicted outcome in untreated group is larger, take mean outcome in untreated group
          if (!to.treat) {
            main.outs = c(main.outs, 0)
            other.outs = c(other.outs, 0)
          }
        }  # end loop through lines in test set
        
        # save mean outcomes from this fold
        main.folds = c(main.folds, mean(main.outs, na.rm = T))
        other.folds = c(other.folds, mean(other.outs, na.rm = T))
        
      }  # end loop through five folds
      
      # save means, variances, and covariances of outcomes across folds for one choice of c1/c2
      avg.main.outs = c(avg.main.outs, mean(main.folds, na.rm = T))
      avg.other.outs = c(avg.other.outs, mean(other.folds, na.rm = T))
      cov.mat = cov(matrix(c(main.folds, other.folds), ncol = 2), use = "complete.obs")
      cov.list[[ctr]] = cov.mat
      parameters[ctr, ] = c(c1, c2)  
      ctr = ctr + 1
      
    }  # end loop through c2s
    
  }   # end loop  through c1s
  
  # determine which c1/c2 maximize main outcome
  opt.ind = which(avg.main.outs == max(avg.main.outs, na.rm = T))
  if (length(opt.ind) > 1) opt.ind = opt.ind[which(avg.other.outs[opt.ind] == max(avg.other.outs[opt.ind]))]
  if (length(opt.ind) > 1) opt.ind = sample(opt.ind, 1)
  
  # select value functions at optimal c1/c2
  if (outcome == "BAR") {
    final.resp = avg.main.outs[opt.ind]
    final.surv = avg.other.outs[opt.ind]
  }
  if (outcome == "Surv") {
    final.resp = avg.other.outs[opt.ind]
    final.surv = avg.main.outs[opt.ind]
  }
  
  # select optimal c1/c2
  final.param = parameters[opt.ind, ]
  
  # select covariance matrix at optimal c1/c2
  final.covariance = cov.list[[opt.ind]]  # note that covariance matrix always has main outcome in upper left
  
  # calculate observed means (over all trts) after subtracting mean in untreated group
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
    var.resp = final.covariance[1, 1]
    var.surv = final.covariance[2, 2]
    cov = final.covariance[1, 2]
  }
  if (outcome == "Surv") {
    obs.resp = mean(clinical.other, na.rm = T)
    obs.surv = mean(clinical, na.rm = T)
    opt.resp = mean(apply(clinical.other, 1, max, na.rm = T)[!is.infinite(apply(clinical.other, 1, max, na.rm = T))], na.rm = T)
    opt.surv = mean(apply(clinical, 1, max, na.rm = T)[!is.infinite(apply(clinical, 1, max, na.rm = T))], na.rm = T)
    var.resp = final.covariance[2, 2]
    var.surv = final.covariance[1, 1]
    cov = final.covariance[1, 2]
  }
  
  res = data.frame(c1 = final.param[1], c2 = final.param[2], mean.response = final.resp, mean.survival = final.surv, 
                   var.response = var.resp, var.survival = var.surv, covariance = cov, observed.resp = obs.resp, 
                   observed.surv = obs.surv, optimal.resp = opt.resp, optimal.surv = opt.surv)
  
  rownames(res) = NULL
  setwd(output_dir)
  output.name = paste(cancer.type, "_", outcome, "_", stri_sub(gene.data.file, 1, -5), outstring, sep = "")
  write.csv(res, output.name)
  
  # return list of results
  return(list(parameters = final.param, mean.response = final.resp, mean.survival = final.surv, 
              covariance = final.covariance, observed.resp = obs.resp, observed.surv = obs.surv, 
              optimal.resp = opt.resp, optimal.surv = opt.surv))
  
}  # end pdx.ql.rf1 function


