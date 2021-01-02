# pdx.ql2.eval.R

# Q-learning with pseudo-values
# evaluate the latent function for super-learning algorithm

library(caret)
library(glmnet)
library(stringi)

pdx.ql2.eval = function(cancer.type, outcome, gene.data.file, input_dir, output_dir, k = 5, seed = 1, strip = T, outstring = "_ql2.csv") {
  
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
  
  
  # create folds for super-learning cross-validation
  set.seed(seed)
  folds.sl = createFolds(1:dim(biomarkers)[1], k = k, list = TRUE, returnTrain = FALSE)
  
  # record matrices for values of latent functions
  output = list()
  
  # loop through folds
  for (f.sl in 1:length(folds.sl)) {
    # clinical.other contains outcomes not specified for estimation
    if (outcome == "BAR") {
      clinical = new.resp.mat 
      clinical.other = new.surv.mat
    }
    if (outcome == "Surv") {
      clinical = new.surv.mat 
      clinical.other = new.resp.mat
    }  
    
    # select training and testing sets
    train.bio = biomarkers[-folds.sl[[f.sl]], ]
    train.clin = clinical[-folds.sl[[f.sl]], ]
    test.bio = biomarkers[folds.sl[[f.sl]], ]
    test.clin = clinical[folds.sl[[f.sl]], ]
    test.clin.other = clinical.other[folds.sl[[f.sl]], ]
    
    # training of the model using optimal c1 and c2 from the previous training step
    
    res.name = paste(cancer.type, "_", outcome, "_", stri_sub(gene.data.file, 1, -5), "_", f.sl, outstring, sep = "")
    res = read.csv(res.name)
    
    c1 = res[1,"c1"]
    c2 = res[1,"c2"]
    
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
    untrt.ind = which(colnames(test.clin) %in% untrt)
    means = apply(as.matrix(test.clin[ , untrt.ind]), 1, mean, na.rm = T)
    test.clin = test.clin - means  # subtract off mean of untreated group in each line
    test.clin = test.clin[ , -untrt.ind]
    
    # same for the test set of the secondary outcome
    untrt.ind = which(colnames(test.clin.other) %in% untrt)
    means = apply(as.matrix(test.clin.other[ , untrt.ind]), 1, mean, na.rm = T)
    test.clin.other = test.clin.other - means  # subtract off mean of untreated group in each line
    test.clin.other = test.clin.other[ , -untrt.ind]
    
    # create treatment tree by clustering
    clusters = hclust(dist(t(train.clin))) # repeat distance matrix after removing untreated columns
    
    # full mouse level data 
    rownames(train.bio) = rownames(train.clin)
    full.bio = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = ncol(train.clin)), ncol = 1), ]
    full.clin = matrix(as.numeric(t(train.clin)), ncol = 1)
    rownames(full.clin) = rownames(full.bio)
    full.trt = matrix(rep(colnames(train.clin), nrow(train.clin)))
    full.clin = cbind(full.clin, full.trt)
    
    # create treatment variables for each step of the tree
    num.steps = dim(clusters$merge)[1]
    merge.steps = clusters$merge
    trt.list = clusters$labels
    new.trt.vars = NULL  # new treatment variables
    # ys: merge < 0 <---- the number of treatment to be merge
    # ys: merge > 0 <---- the number of step (group of treatments) to be merge
    # ys: 1 means the left subgroup, -1 means the right subgroup, NA means not in the group
    for (j in 1:num.steps) {
      temp = rep(NA, dim(full.clin)[1])
      merge1 = merge.steps[j, 1]; merge2 = merge.steps[j, 2]
      if (merge1 < 0) {
        temp[which(full.clin[ , 2] == trt.list[-merge1])] = 1
      }
      if (merge2 < 0) {
        temp[which(full.clin[ , 2] == trt.list[-merge2])] = -1
      }
      # ys: label the treatments in the merging cluster
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
    rownames(new.trt.vars) = rownames(full.clin)
    
    # select trt variables for c2 steps down tree
    # ys: which are the last c2 columns
    trt.vars = new.trt.vars[ , rev(rev(1:dim(new.trt.vars)[2])[1:c2])]
    trt.vars = as.matrix(trt.vars)
    
    # Q-learning up the tree
    list.models = list()  # store models at each step
    list.pred = list()  # store predicted values at each step
    for (d in 1:dim(trt.vars)[2]) {
      
      X = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = 2), ncol = 1), ]
      
      # ys: each model has two replicates in X and Y, which means the left group and right group
      
      Y = matrix(NA, nrow = nrow(X), ncol = 2)
      Y[ , 2] = rep(c(0, 1), nrow(train.bio))
      rownames(Y) = rownames(X)
      
      # for first step, outcomes are mean outcomes in root nodes
      # ys: here root nodes means: the two lowest subgroups
      if (d == 1) {
        for (g in 1:dim(Y)[1]) {
          # ys: mean of responses: 
          # 1. with treatments belong to the (left or right) subgroup
          # 2. on the current model (mouse)
          Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(trt.vars) == rownames(Y)[g]), 1]), na.rm = T)
        }
      }
      
      # for other steps, outcomes are the maximum over outcomes on lower steps of the tree (if a previous decision has been made)
      if (d > 1) {
        for (g in 1:dim(Y)[1]) {
          
          # most reason previous step of tree where decision was made
          # ys: index <- number of step when the left/right subgroup was merged
          index = max(which(!is.na(as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2]), 1:(d - 1)])[1,])))
          # ys: if infinite, then both subgroups are leaf nodes
          # calculate the mean as same as if d == 1
          if (is.infinite(index)) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
          # ys: mean observed response of all the treatments in the larger sub-subgroup of left/right subgroup
          # here list.pred is only used to tell which sub-subgroup is larger
          # !!!!!!!!!!!!question: shouldn't be the maximum leaf node???
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
      
      # ys: A is the action and Y is the response 
      A = Y[ , 2]
      Y = Y[ , 1]
      
      # fit model for Q-learning
      tempx = cbind(X, A, X * A) # ys: all the predictors
      tempy = Y
      # ys: If either Rj0(t) or Rj1(t) are missing for a given line, 
      # meaning that PDX line j did not receive that particular treatment 
      # and did not observe any corresponding response for this treatment, 
      # these observations are removed before fitting the regression model
      if (sum(is.na(Y) > 0)) tempy = tempy[-which(is.na(Y))]
      if (sum(is.na(Y) > 0)) tempx = tempx[-which(is.na(Y)), ]
      # ys: fit the model, here "gaussian" means the type of data is quantitative
      fit = cv.glmnet(tempx, tempy, family = "gaussian", lambda = 2^seq(-10, -4, 0.25), nfolds = 5)
      lam = fit$lambda.min
      #print(lam)
      list.models[[d]] = fit
      
      # create pseudo-values for next step
      preds = matrix(NA, nrow = nrow(train.bio), ncol = 3)
      for (t in 1:dim(train.bio)[1]) {
        # ys: predictors of left and right subgroups (same X, different A)
        predx0 = c(train.bio[t, ], 0, rep(0, length(train.bio[t, ]))) 
        predx1 = c(train.bio[t, ], 1,  train.bio[t, ])
        # ys: preds[t, 1] <- model name
        # preds[t, 2] <- larger value of predicted response from left and right subgroups
        # preds[t, 3] <- index (left or right)
        preds[t, 1] = rownames(train.bio)[t]
        preds[t, 2] = max(predict(fit$glmnet.fit, matrix(predx0, nrow = 1), s = lam), predict(fit$glmnet.fit, matrix(predx1, nrow = 1), s = lam))
        preds[t, 3] = as.numeric(predict(fit$glmnet.fit, matrix(predx0, nrow = 1), s = lam) < predict(fit$glmnet.fit, matrix(predx1, nrow = 1), s = lam))
      }
      list.pred[[d]] = preds
    }  # end loop through steps in tree  
    
    # record values of latent functions
    O = matrix(0, nrow(test.bio), ncol(test.clin))
    
    # for each test model i
    for(i in 1:nrow(test.bio)){
      # find the treatment group for treatment t
      for(t in 1:ncol(test.clin)){
        pos = t
        isleft = TRUE
        ismerged = -1
        for(l in 1:nrow(merge.steps)){
          if((ismerged*pos) %in% merge.steps[l,]){
            if((ismerged*pos) == merge.steps[l,1])
              isleft = TRUE
            else
              isleft = FALSE
            ismerged = 1
            pos = l
            if(pos > nrow(merge.steps) - c2)
              break
          }
        }
        group = pos - (nrow(merge.steps) - c2)
        temp.fit = list.models[[group]]
        lam = temp.fit$lambda.min
        
        # predictors accroding to left or right
        if(isleft){
          predx = c(test.bio[i, ], 1,  test.bio[i, ])
        }else{
          predx = c(test.bio[i, ], 0, rep(0, length(test.bio[i, ])))
        }
        
        # latent function as the prediction of the treatment group
        O[i,t] = predict(temp.fit$glmnet.fit, matrix(predx, nrow = 1), s = lam)
      }
    }
    
    # add the untreated group back
    for(ut in untrt.ind){
      O = cbind(O[,1:(ut - 1)], rep(0, nrow(O)), O[,-(1:(ut - 1))])
    }
    
    colnames(O) = colnames(clinical)
    rownames(O) = rownames(test.clin)
    
    output[[f.sl]] = O
    
  }  
  
  return(output)
}











