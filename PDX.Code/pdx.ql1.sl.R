# pdx.ql1.R

# Q-learning with no pseudo-values

library(caret)
library(glmnet)
library(stringi)

# need to supply to the function a cancer type, one of 
# "BRCA", "CM", "CRC", "NSCLC", or "PDAC"
# and an outcome, one of "BAR" for best average response 
# or "Surv" for time to doubling
# gene.data.file is the name of a csv file with gene data
# input_dir and output_dir are directories (both should be main directory)
# k is the number of folds for cross-validation
# also need to supply a seed to make training/testing sets for cross-validation 
# if strip is TRUE, we assume the first column of gene data contains row names and is dropped
# outstring is an identifier for the output csv file
pdx.ql1 = function(cancer.type, outcome, gene.data.file, input_dir, output_dir, k = 5, seed = 1, strip = T, outstring = "_ql1.csv") {
  
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
    
    biomarkers.train = biomarkers[-folds.sl[[f.sl]], ]
    clinical.train = clinical[-folds.sl[[f.sl]], ]
    clinical.other.train = clinical.other[-folds.sl[[f.sl]], ]
    
    # create folds for cross-validation
    #set.seed(seed)
    #folds = createFolds(1:dim(biomarkers.train)[1], k = k, list = TRUE, returnTrain = FALSE)
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
      
      # ys: if c2s is not specified, the max number for each c1 are tried
      c2s = seq(1:(ncol(clinical) - c1 - 2))
      
      for (c2 in c2s) {  
        
        # store primary and secondary value functions across folds
        main.folds = NULL
        other.folds = NULL
        
        # loop through folds
        for (f in 1:length(folds)) {
          
          # select training and testing sets
          # train.bio = biomarkers.train[-folds[[f]], ]
          # train.clin = clinical.train[-folds[[f]], ]
          # test.bio = biomarkers.train[folds[[f]], ]
          # test.clin = clinical.train[folds[[f]], ]
          # test.clin.other = clinical.other.train[folds[[f]], ]
          
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
                # "which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2])" 
                # gives the lines of treatments in left/right group for the current model
                # "as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2]), 1:(d - 1)])" 
                # gives the matrix form of the history of the above treatments
                # "max(which(!is.na(as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2]), 1:(d - 1)])[1,])))" 
                # is the most recent merging time that the first treatment(in the subgroup) appears, which must be the merging time of whole left/right subgroup
                index = max(which(!is.na(as.matrix(trt.vars[which(rownames(Y)[g] == rownames(trt.vars) & trt.vars[ , d] == Y[g, 2]), 1:(d - 1)])[1,])))
                # ys: if infinite, then subgroup is leaf node
                # calculate the mean as same as if d == 1
                if (is.infinite(index)) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
                # ys: mean observed response of all the treatments in the larger sub-subgroup of left/right subgroup
                # here list.pred is only used to tell which sub-subgroup is larger
                # "list.pred[[index]][which(list.pred[[index]][ , 1] == rownames(Y)[g]), 3]" tells the choice of left/right based on predictions
                # "rownames(full.clin) == rownames(Y)[g])" is the current model
                # "which(trt.vars[ , index] == list.pred[[index]][which(list.pred[[index]][ , 1] == rownames(Y)[g]), 3] & rownames(full.clin) == rownames(Y)[g])" is the lines of the chosed treatments in current model
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
            print(lam)
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
          
          # test on validation set 
          main.outs = NULL  # save mean outcomes for each line among mice treated consistent with treatment rule
          other.outs = NULL
          rownames(test.bio) = rownames(test.clin)
          
          # loop through lines in test set
          for (t in 1:dim(test.bio)[1]) {
            
            # determine sequence of decisions for each line
            moves = NULL
            # ys: start from the root node
            for (i in rev(1:dim(trt.vars)[2])) {
              # ys: predictors of left and right subgroups (same X, different A)
              predx0 = c(test.bio[t, ], 0, rep(0, length(test.bio[t, ]))) 
              predx1 = c(test.bio[t, ], 1,  test.bio[t, ])
              temp.fit = list.models[[i]]
              lam = temp.fit$lambda.min
              # ys: move to the side with larger prediction
              moves = c(moves, as.numeric(predict(temp.fit$glmnet.fit, matrix(predx0, nrow = 1), s = lam) < predict(temp.fit$glmnet.fit, matrix(predx1, nrow = 1), s = lam)))
            }
            # ys: start at root node
            cur.step = merge.steps[nrow(merge.steps), ]
            trt.ind = NULL  # save indices (in trt.list) of those trts that are consistent with decision rule
            
            temp = rep(NA, nrow(merge.steps) - length(moves))
            tmoves = c(temp, rev(moves)) 
            # ys: tell the move at each merging point
            # if the merging point is within the subgroup of treatment
            # then the value is NA (don't need to move any more)
            
            # determine root node for each line in testing set
            keeplooping = T
            cur.move = tmoves[length(tmoves)]
            while (keeplooping) {
              if (cur.move == 1) next.step = cur.step[1]
              if (cur.move == 0) next.step = cur.step[2] 
              # ys: leaf node
              if (next.step < 0) {
                trt.ind = c(trt.ind, -next.step)
                keeplooping = F
              }
              if (next.step > 0) { 
                cur.step = merge.steps[next.step, ]
                cur.move = tmoves[next.step]
                # ys: f the merging point is within the subgroup of treatment, then stop
                # the treatment list is given by the following function get.trt.list
                if (is.na(cur.move)) keeplooping = F
              }
            }
            
            # recursive function to construct list of indices for trts in root node
            get.trt.list = function(cur.step) {
              trt.ind = NULL
              # ys: if leaf node, then add this into the list
              # ys: if a group, then add all the treatments into the list
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
            predx0 = c(test.bio[t, ], 0, rep(0, length(test.bio[t, ]))) 
            predx1 = c(test.bio[t, ], 1,  test.bio[t, ])
            temp.fit = list.models[[dim(trt.vars)[2]]]
            lam = temp.fit$lambda.min
            # ys: the mean value of untreated group is 0
            to.treat = max(predict(temp.fit$glmnet.fit, matrix(predx0, nrow = 1), s = lam), predict(temp.fit$glmnet.fit, matrix(predx1, nrow = 1), s = lam)) > 0
            
            # if predicted outcome from decision rule is larger, take mean outcome in root node
            if (to.treat) {
              # ys: average for all the response for the current model and treatments from decision rule
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
    # if multiple groups reaches same max main outcome, then choose max other outcome
    if (length(opt.ind) > 1) opt.ind = opt.ind[which(avg.other.outs[opt.ind] == max(avg.other.outs[opt.ind]))]
    # if there are still multiple groups, then randomly choose
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
    # adjust output function for now
    output.name = paste(cancer.type, "_", outcome, "_", stri_sub(gene.data.file, 1, -5), "_", f.sl, outstring, sep = "")
    write.csv(res, output.name)
    
  }
  
  
  
  # return list of results
  return(list(parameters = final.param, mean.response = final.resp, mean.survival = final.surv, 
              covariance = final.covariance, observed.resp = obs.resp, observed.surv = obs.surv, 
              optimal.resp = opt.resp, optimal.surv = opt.surv))
  
}  # end pdx.ql1 function

