# pdx.dl.owl.linear.smooth.R

# OWL with linear decision function, random forest smoothed outcomes, and deep learning reduced data

library(caret)
library(stringi)

# need to supply to the function a cancer type, one of 
# "BRCA", "CM", "CRC", "NSCLC", or "PDAC"
# and an outcome, one of "BAR" for best average response 
# or "Surv" for time to doubling
# gene.data.file is the name of a csv file with gene data
# numgenes is the number of genes to use for smoothing
# k is the number of folds for cross-validation
# also need to supply a seed to make training/testing sets for cross-validation 
# c1s and c2s are the tuning parameters to try; if c2s is not specified, the max number for each c1 are tried
# if strip = T the first column of gene.data.file is assumed to be rownames and is stripped off
# outstring is an identifier for the output csv file
pdx.dl.owl.linear.smooth = function(cancer.type, outcome, gene.data.file, numgenes, input_dir, output_dir, c1s = c(0), c2s = NA, k = 5, seed = 1, strip = T, outstring = "_owllinearsmoothdl.csv") {
  
  setwd(input_dir)
  load("split.cm.data.rda")
  
  load("trts.by.cancer.rda")
  
  # random forest predicted values -- use these as outcomes when estimating decision rule
  load("pred_vals_rf.rda")
  
  # extract clinical data for given cancer type
  if (cancer.type == "BRCA") dat = split.cm.data$BRCA
  if (cancer.type == "CM") dat = split.cm.data$CM
  if (cancer.type == "CRC") dat = split.cm.data$CRC
  if (cancer.type == "NSCLC") dat = split.cm.data$NSCLC
  if (cancer.type == "PDAC") dat = split.cm.data$PDAC
  clinical = dat[ , 1:17]
  clinical$RespScaled = -clinical$RespScaled  # reverse sign of best average response -- this way larger values are better

  ind1 = which(cancer.type == names(trts.by.cancer))
  ind2 = which(numgenes == c(50, 100, 500, 1000))
  
  clinical$new.resp = pred_vals_rf[[ind1]][[ind2]][[1]]
  clinical$new.surv = pred_vals_rf[[ind1]][[ind2]][[2]]

  biomarkers = read.csv(gene.data.file)
  if (strip) biomarkers = biomarkers[ , -1]
  ngene = dim(biomarkers)[2]
  
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
  
  # RF predicted values
  smooth.resp.mat = matrix(NA, nrow = length(unique(clinical$Model)), ncol = length(unique(clinical$Treatment)))
  smooth.surv.mat = matrix(NA, nrow = length(unique(clinical$Model)), ncol = length(unique(clinical$Treatment)))
  rownames(smooth.resp.mat) = unique(clinical$Model); colnames(smooth.resp.mat) = unique(clinical$Treatment)
  rownames(smooth.surv.mat) = unique(clinical$Model); colnames(smooth.surv.mat) = unique(clinical$Treatment)
  for (dim1 in 1:nrow(smooth.resp.mat)) {
    for (dim2 in 1:ncol(smooth.resp.mat)) {
      row = which(clinical$Model == rownames(smooth.resp.mat)[dim1] & clinical$Treatment == colnames(smooth.resp.mat)[dim2])
      if (length(row) != 0) {
        smooth.resp.mat[dim1, dim2] = clinical$new.resp[row]
        smooth.surv.mat[dim1, dim2] = clinical$new.surv[row]
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
    
    # smooth.clinical.other contains outcomes not specified for estimation
    if (outcome == "BAR") {
      smooth.clinical = smooth.resp.mat 
      smooth.clinical.other = smooth.surv.mat
    }
    if (outcome == "Surv") {
      smooth.clinical = smooth.surv.mat 
      smooth.clinical.other = smooth.resp.mat
    }  
    
    biomarkers.train = biomarkers[-folds.sl[[f.sl]], ]
    clinical.train = clinical[-folds.sl[[f.sl]], ]
    clinical.other.train = clinical.other[-folds.sl[[f.sl]], ]
    smooth.clinical.train = smooth.clinical[-folds.sl[[f.sl]], ]
    smooth.clinical.other.train = smooth.clinical.other[-folds.sl[[f.sl]], ]
    
    # create folds for cross-validation
    #set.seed(seed)
    folds = createFolds(1:dim(biomarkers.train)[1], k = k, list = TRUE, returnTrain = FALSE)
    
    # c1 is the number of trt's to group with untreated, c2 is the number of steps to take down the tree
    avg.main.outs = NULL  # store primary value functions across c1 and c2
    avg.other.outs = NULL  # store secondary value functions across c1 and c2
    parameters = matrix(NA, nrow = length(c1s) * ncol(clinical), ncol = 2)  # store pairs of c1 and c2
    colnames(parameters) = c("c1", "c2")
    cov.list = list()  # to store covariances of value functions
    ctr = 1  # count number of times through inner loop
    print ("Begin the loops of ntrt and nodes")
    # loop through parameters
    for (c1 in c1s) {
      
      if (is.na(c2s)) c2s = seq(1, (ncol(clinical) - c1 - 2))
      
      print(sprintf("  c1 = %d", c1))
      
      for (c2 in c2s) {  
        
        print(sprintf("    c2 = %d", c2))
        # store primary and secondary value functions across folds
        main.folds = NULL
        other.folds = NULL
        
        # loop through folds
        for (f in 1:length(folds)) {
          
          # select training and testing sets
          train.bio = biomarkers.train[-folds[[f]], ]
          train.clin = smooth.clinical.train[-folds[[f]], ]
          test.bio = biomarkers.train[folds[[f]], ]
          test.clin = clinical.train[folds[[f]], ]
          test.clin.other = clinical.other.train[folds[[f]], ]
          
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
          
          # select trt variables for c2 steps down tree
          row.names(new.trt.vars) = row.names(full.clin)
          trt.vars = new.trt.vars[ , rev(rev(1:dim(new.trt.vars)[2])[1:c2])]
          trt.vars = as.matrix(trt.vars)
          
          # OWL up the tree
          moves.mat = NULL  # this will be nrow(test.bio) by ncol(trt.vars) -- each row contains the moves up the tree for one line in test set
          list.pred = list()
          for (d in 1:dim(trt.vars)[2]) {
            
            X = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = 2), ncol = 1), ]
            
            Y = matrix(NA, nrow = nrow(X), ncol = 2)
            Y[ , 2] = rep(c(-1, 1), nrow(train.bio))
            rownames(Y) = rownames(X)
            
            # for first step, outcomes are mean outcomes in leaf nodes
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
                if (is.infinite(index)) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , d] == Y[g, 2] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
                if (!is.infinite(index)) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , index] == list.pred[[index]][which(list.pred[[index]][ , 1] == rownames(Y)[g]), 3] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
              }
            }
            
            A = Y[ , 2]
            Y = Y[ , 1]
            
            # fit model for OWL
            temp.mat = cbind(Y, A, X)
            trainname = sprintf("dlsmoothlinear_train_%s_%s_%i.csv", cancer.type, outcome, ngene)
            testname = sprintf("dlsmoothlinear_test_%s_%s_%i.csv", cancer.type, outcome, ngene)
            write.table(temp.mat, file = trainname, col.names = FALSE, row.names = FALSE, sep = ",")
            write.table(test.bio, testname, col.names = F, row.names = F, sep = ",")
            cmdstring = sprintf("python3 owl.temp.linear.py %s %s", trainname, testname)
            owl.res = system(cmdstring, intern = T)
            temp.res = read.csv(sprintf("temp_owl_%s", testname), header = F)  # this is the file that python.owl.temp returns
            temp.res = as.matrix(temp.res * 2 - 3)
            print(dim(temp.res))
            print(dim(moves.mat))
            moves.mat = cbind(moves.mat, temp.res)
            
            # predicted best treatment for training set
            # we're bowworing code from the QL files -- predicted treatments are in column 3
            traindatmoves = read.csv(sprintf("temp_owl_%s", trainname), header=F)
            traindatmoves = as.matrix(traindatmoves * 2 - 3)
            traindatmoves = traindatmoves[which(1:nrow(traindatmoves) %% 2 == 1), ]
            traindatmoves = as.matrix(traindatmoves)
            preds = matrix(NA, nrow = nrow(train.bio), ncol = 3)
            preds[ , 1] = rownames(train.bio)
            preds[ , 3] = traindatmoves
            list.pred[[d]] = preds
          }  # end loop through steps in tree
          
          # test on validation set 
          main.outs = NULL  # save mean outcomes for each line among mice treated consistent with treatment rule
          other.outs = NULL
          rownames(test.bio) = rownames(test.clin)
          
          # moves for each row in test set
          moves.all = moves.mat
          if (ncol(moves.mat) > 1) moves.all = t(apply(moves.mat, 1, rev))  # when we fill moves.mat we are doing it from the bottom of the tree up, not top down
          moves.all[which(moves.all == -1)] = 0
          
          # loop through lines in test set 
          for (t in 1:dim(test.bio)[1]) {
            
            # moves for line t in test set (these are the moves that would be taken under decision rule)
            moves = moves.all[t, ]
            
            cur.step = merge.steps[nrow(merge.steps), ]
            trt.ind = NULL  # save indices (in trt.list) of those trts that are consistent with decision rule
  
            temp = rep(NA, nrow(merge.steps) - length(moves))
            tmoves = c(temp, rev(moves))
            
            # determine leaf node for each line in testing set 
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
            
            # recursive function to construct list of indices for trts in leaf node
            get.trt.list = function(cur.step) {
              trt.ind = NULL
              if (cur.step[1] < 0) trt.ind = c(trt.ind, -cur.step[1])
              if (cur.step[1] > 0) trt.ind = c(trt.ind, get.trt.list(merge.steps[cur.step[1], ])) 
              if (cur.step[2] < 0) trt.ind = c(trt.ind, -cur.step[2])
              if (cur.step[2] > 0) trt.ind = c(trt.ind, get.trt.list(merge.steps[cur.step[2], ])) 
              return(trt.ind)
            }  # end recursive function
            
            # if leaf node is not a single trt, get list of trts in leaf node
            if (length(trt.ind) == 0) trt.ind = get.trt.list(cur.step)
            
            # list of trts in leaf node
            treatments = trt.list[trt.ind]
            
            # take mean outcome in leaf node
            main.outs = c(main.outs, mean(test.clin[which(rownames(test.clin) == rownames(test.bio)[t]), which(colnames(test.clin) %in% treatments)]))
            other.outs = c(other.outs, mean(test.clin.other[which(rownames(test.clin.other) == rownames(test.bio)[t]), which(colnames(test.clin.other) %in% treatments)]))
          }  # end loop through lines in test set
          
          # determine if any lines should have been left untreated by fitting OWL one more time
          X = train.bio[matrix(apply(matrix(1:nrow(train.bio), ncol = 1), 1, rep, times = 2), ncol = 1), ]
          Y = matrix(NA, nrow = nrow(X), ncol = 2)
          Y[ , 2] = rep(c(0, 1), nrow(train.bio))
          rownames(Y) = rownames(X)
          for (g in 1:dim(Y)[1]) {
            if (Y[g, 2] == 0) Y[g, 1] = 0  # 0 is untreated
            if (Y[g, 2] == 1) Y[g, 1] = mean(as.numeric(full.clin[which(trt.vars[ , ncol(trt.vars)] == list.pred[[ncol(trt.vars)]][which(list.pred[[ncol(trt.vars)]][ , 1] == rownames(Y)[g]), 3] & rownames(full.clin) == rownames(Y)[g]), 1]), na.rm = T)
          }
          
          A = Y[ , 2]
          Y = Y[ , 1]
          
          A[which(A == 0)] = -1  # trt coded as 1/-1 for OWL
          
          # fit model for OWL
          temp.mat = cbind(Y, A, X)
          trainname = sprintf("dlsmoothlinear_train_%s_%s_%i.csv", cancer.type, outcome, ngene)
          testname = sprintf("dlsmoothlinear_test_%s_%s_%i.csv", cancer.type, outcome, ngene)
          write.table(temp.mat, file=trainname, col.names=FALSE, row.names=FALSE, sep=",")
          write.table(test.bio, testname, col.names=F, row.names=F, sep=",")
          cmdstring = sprintf("python3 owl.temp.linear.py %s %s", trainname, testname)
          owl.res = system(cmdstring, intern=T)
          temp.res = read.csv(sprintf("temp_owl_%s", testname), header=F)  # this is the file that python.owl.temp returns
          temp.res = as.matrix(temp.res*2 - 3)
          
          # for any mice who should have been left untreated, replace outcome with mean in untreated group
          main.outs[which(temp.res == -1)] = 0
          
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
        
      }  # end loop through c1s
      
    }   # end loop  through c2s
    
    # determine which c1/c2 maximize main outcome
    opt.ind = which(avg.main.outs == max(avg.main.outs, na.rm = T))
    if (length(opt.ind) > 1) opt.ind = which(avg.other.outs[opt.ind] == max(avg.other.outs[opt.ind]))
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
    # find the c1 closest treatments to "untreated"
    dist_mat = as.matrix(dist(t(clinical)))
    col_ind = which(colnames(dist_mat) == "untreated")
    ordered_dist_mat = dist_mat[order(dist_mat[ , col_ind]) , col_ind]
    untrt = names(ordered_dist_mat[1:(1 + final.param[1])])
    # observed mean for primary outcome
    untrt.ind = which(colnames(clinical) %in% untrt)
    means = apply(as.matrix(clinical[ , untrt.ind]), 1, mean, na.rm = T)
    clinical = clinical - means  # subtract off mean of untreated group in each line
    clinical = clinical[ , -untrt.ind]
    # observed mean for secondary outcome
    untrt.ind = which(colnames(clinical.other) %in% untrt)
    means = apply(as.matrix(clinical.other[ , untrt.ind]), 1, mean, na.rm = T)
    clinical.other = clinical.other - means  # subtract off mean of untreated group in each line
    clinical.other = clinical.other[ , -untrt.ind]
    
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
    output.name = paste(cancer.type, "_", outcome, "_", stri_sub(gene.data.file, 1, -5), "_", f.sl, outstring, sep = "")
    write.csv(res, output.name)
  }
  # return list of results
  return(list(parameters = final.param, mean.response = final.resp, mean.survival = final.surv,
              covariance = final.covariance, observed.resp = obs.resp, observed.surv = obs.surv,
              optimal.resp = opt.resp, optimal.surv = opt.surv))
  
}  # end pdx.dl.owl.linear.smooth function
