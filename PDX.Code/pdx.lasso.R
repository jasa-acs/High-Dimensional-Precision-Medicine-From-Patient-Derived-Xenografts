# Function takes a specified biomarker datafile and passed clinical variable dataset, 
#   then carries out the LASSO analysis using glmnet. Refer to Michael's RF code.  Namely:
#   1) Process data and divide into 5 folds. For each outcome variable:
#       1a) For each fold: fit lasso to the rest of the data, then use the fitted model to predict optimal treatments 
#               in the held-out fold. Evaluate the decision rule via value function. Compute mean predicted value, 
#               SD predicted value, mean observed value, optimal value, and adjusted optimal value. Store results.
#       1b) Fit LASSO to entire data. Store variable importance.

#infile,clin,
#verbose=T

lasso.value <- function(infile,clin,verbose=F) {
  
  dat <- read.csv(infile,header=T)
  keepcols <- grepl(".cn|.rna|.mut",colnames(dat))
  dat <- dat[,keepcols]
  # Want to minimize BAR, so we plug -BAR into the maximization analysis
  BAR <- -clin[,"RespScaled"]
  TTD <- clin[,"logSurvScaled"]
  trt <- clin[,"Treatment"]
  untrt.level <- which(levels(trt)=="untreated")
  trt <- relevel(trt,ref="untreated")
  lines <- unique(clin[,"Model"])
  ntrt <- length(levels(trt))
  trt.mat <- model.matrix(~trt-1)[,-1]
  colnames(trt.mat) <- levels(trt)[which(levels(trt)!="untreated")]
  X <- cbind(trt.mat,dat)
  set.seed(1)
  folds = createFolds(1:length(lines), k = 5, list = TRUE, returnTrain = FALSE)
  for (i in 1:length(folds)) {
    folds[[i]] <- lines[folds[[i]]]
  }
  value <- sapply(c("trained_BAR","trained_TTD"),function(x) NULL)
  for (i in 1:2) {
    value[[i]] <- sapply(c("eval_BAR","eval_TTD"),function(x) NULL)
    for (j in 1:2) {
      value[[i]][[j]] <- rep(NA,5)
      names(value[[i]][[j]]) <- c("mean_value","sd_value","obs_value","opt_value","opt_value_adj")
    }
  }
  importance <- sapply(c("trained_BAR","trained_TTD"),function(x) NULL)
  pred.cent.BAR <- rep(NA,dim(clin)[1])
  pred.cent.TTD <- rep(NA,dim(clin)[1])
  for (t in 1:2) {                                    # begin loop over outcome variables
    outcomes <- c("RespScaled","logSurvScaled")
    Y <- clin[,outcomes[t]]
    Y2 <- clin[,outcomes[((t==1)+1)]]
    for (line in lines) {
      ind <- which(clin[,"Model"]==line)
      if (length(ind)>1) {
        subY <- Y[ind]
        subY2 <- Y2[ind]
        untrtY <- Y[which(clin[,"Model"]==line & clin[,"Treatment"]=="untreated")]
        untrtY2 <- Y2[which(clin[,"Model"]==line & clin[,"Treatment"]=="untreated")]
        if(length(untrtY)>0) {
          subY <- subY-untrtY
          subY2 <- subY2-untrtY2
        }
        if (length(untrtY)==0) {
          subY <- subY-min(subY)
          subY2 <- subY2 - min(subY2)
        }
        Y[ind] <- subY
        Y2[ind] <- subY2
      }
    }
    # Want to minimize BAR
    if (t==1) {
      Y <- -Y
    }
    if (t==2) {
      Y2 <- -Y2
    }
    if (verbose) {print(outcomes[t])}
    value_by_fold_BAR <- rep(NA,length(folds))
    value_by_fold_TTD <- rep(NA,length(folds))
    opt_val_by_fold_BAR <- rep(NA,length(folds))
    opt_val_by_fold_TTD <- rep(NA,length(folds))
    for (i in 1:length(folds)) {                     # begin loop over CV folds
      # Subset to -ith fold and train random forest
      cvind <- !(clin[,"Model"]%in%folds[[i]])
      X.temp <- X[cvind,]
      Y.temp <- Y[cvind]
      X.trt.temp = X.temp[, 1:(ntrt-1)]
      X.bio.temp = X.temp[, -(1:(ntrt-1))]
      Xt.temp = matrix(NA, ncol = ncol(X.trt.temp) * ncol(X.bio.temp), nrow = nrow(X.temp))
      colindex = 1
      for(it in 1:ncol(X.trt.temp)){
        for(ix in 1:ncol(X.bio.temp)){
          Xt.temp[,colindex] = X.trt.temp[,it] * X.bio.temp[,ix]
          colindex = colindex + 1
        }
      }
      X.temp.int <- cbind(X.temp, Xt.temp)
      #lasso_md <- glmnet(as.matrix(X.temp), Y.temp, alpha=1)
      cvlasso  = cv.glmnet(as.matrix(X.temp.int), Y.temp, alpha=1, type.measure = "mse", nfolds = 5)
      lamopt = cvlasso$lambda.min
      #rf <- randomForest(X.temp,Y.temp)
      
      # Use lasso to predict outcome for each treatment in the ith fold and store predicted best treatments
      pred.trt.outcome <- matrix(NA,nrow=length(folds[[i]]),ncol=ntrt)
      rownames(pred.trt.outcome) <- folds[[i]]
      colnames(pred.trt.outcome) <- levels(trt)
      for (j in 1:ntrt) {
        trt.name <- levels(trt)[j]
        X.trt <- X[!cvind & trt==trt.name,]
        if (nrow(X.trt)>0) {
          rownames(X.trt) <- clin[!cvind & trt==trt.name,"Model"]
          
          X.trt.trt.temp = X.trt[, 1:(ntrt-1)]
          X.trt.bio.temp = X.trt[, -(1:(ntrt-1))]
          Xt.trt.temp = matrix(NA, ncol = ncol(X.trt.trt.temp) * ncol(X.trt.bio.temp), nrow = nrow(X.trt))
          colindex = 1
          for(it in 1:ncol(X.trt.trt.temp)){
            for(ix in 1:ncol(X.trt.bio.temp)){
              Xt.trt.temp[,colindex] = X.trt.trt.temp[,it] * X.trt.bio.temp[,ix]
              colindex = colindex + 1
            }
          }
          X.trt.int <- cbind(X.trt, Xt.trt.temp)
          
          pred <- predict(cvlasso, as.matrix(X.trt.int), s=lamopt)
          names(pred) = row.names(X.trt)
          for (k in 1:length(pred)) {
            pred.trt.outcome[names(pred)[k],trt.name] <- pred[k]
          }
        } else {
          pred.trt.outcome[names(pred)[k],trt.name] <- NA
        }
        
      }
      best.trts <- rep("",length(folds[[i]]))
      names(best.trts) <- folds[[i]]
      for (j in 1:length(best.trts)) {
        best.trts[j] <- colnames(pred.trt.outcome)[which(pred.trt.outcome[folds[[i]][j],]==
                                                           max(pred.trt.outcome[folds[[i]][j],],na.rm=T))]
      }
      # Store predicted outcome for held out mice
      predY <- rep(NA,length(Y))
      for (line in folds[[i]]) {
        for (tr in levels(trt)) {
          predY[clin[,"Model"]==line & trt==tr] <- pred.trt.outcome[line,tr]
        } 
      }
      replace.ind <- !is.na(predY)
      if (t==1) {
        pred.cent.BAR[replace.ind] <- predY[replace.ind]
      }
      if (t==2) {
        pred.cent.TTD[replace.ind] <- predY[replace.ind]
      }
      # Evaluate the decision rule using the ith fold outcomes
      decis_BAR <- rep(NA,length(folds[[i]]))
      decis_TTD <- rep(NA,length(folds[[i]]))
      if (t==1) {
        for (j in 1:length(folds[[i]])) {
          decis_BAR[j] <- Y[clin[,"Model"]==names(best.trts)[j] & trt==best.trts[j]]
          decis_TTD[j] <- Y2[clin[,"Model"]==names(best.trts)[j] & trt==best.trts[j]]
        }
        value_by_fold_BAR[i] <- mean(decis_BAR)
        value_by_fold_TTD[i] <- mean(decis_TTD)
      }
      if (t==2) {
        for (j in 1:length(folds[[i]])) {
          decis_BAR[j] <- Y2[clin[,"Model"]==names(best.trts)[j] & trt==best.trts[j]]
          decis_TTD[j] <- Y[clin[,"Model"]==names(best.trts)[j] & trt==best.trts[j]]
        }
        value_by_fold_BAR[i] <- mean(decis_BAR)
        value_by_fold_TTD[i] <- mean(decis_TTD)
      }
      
      
      if (verbose) {print(paste("fold",i,"done"))}
    }         # End loop over CV folds
    # Find the optimal value assuming oracle knowledge
    if (t==1) {
      max_values_BAR <- sapply(lines, function(x) max(Y[clin[,"Model"]==x],na.rm=TRUE))
      max_values_TTD <- sapply(lines, function(x) max(Y2[clin[,"Model"]==x],na.rm=TRUE))
      obs_values_BAR <- Y
      obs_values_TTD <- Y2
    }
    if (t==2) {
      max_values_BAR <- sapply(lines, function(x) max(Y2[clin[,"Model"]==x],na.rm=TRUE))
      max_values_TTD <- sapply(lines, function(x) max(Y[clin[,"Model"]==x],na.rm=TRUE))
      obs_values_BAR <- Y2
      obs_values_TTD <- Y
    }
    # Fit LASSO to full data and compute variable importance
    #rf.full <- randomForest(X,Y,importance=TRUE)
    X.trt.temp = X[, 1:(ntrt-1)]
    X.bio.temp = X[, -(1:(ntrt-1))]
    Xt.temp = matrix(NA, ncol = ncol(X.trt.temp) * ncol(X.bio.temp), nrow = nrow(X))
    colindex = 1
    for(it in 1:ncol(X.trt.temp)){
      for(ix in 1:ncol(X.bio.temp)){
        Xt.temp[,colindex] = X.trt.temp[,it] * X.bio.temp[,ix]
        colindex = colindex + 1
      }
    }
    bio_index = rep(1:ncol(X.bio.temp), ncol(X.trt.temp))
    X.temp.int <- cbind(X, Xt.temp)
    cvlasso  = cv.glmnet(as.matrix(X.temp.int), Y, alpha=1, type.measure = "mse", nfolds = 5)
    lamopt = cvlasso$lambda.min
    lasso_imp_raw = coef(cvlasso, s=lamopt)[-c(1:(1+ncol(X)))] ## remove the intercept
    lasso_imp = tapply(lasso_imp_raw, INDEX=bio_index, FUN=function(x) max(abs(x)))
    #importance[[t]] <- rf.full$importance[,1]
    importance[[t]] <- lasso_imp
    if (verbose) {print("importance computed")}
    
    value[[t]][[1]][1] <- mean(value_by_fold_BAR)
    value[[t]][[1]][2] <- sd(value_by_fold_BAR)
    value[[t]][[1]][3] <- mean(obs_values_BAR)
    value[[t]][[1]][4] <- mean(max_values_BAR)
    value[[t]][[1]][5] <- mean(max_values_BAR)/sqrt(log(ntrt))
    value[[t]][[2]][1] <- mean(value_by_fold_TTD)
    value[[t]][[2]][2] <- sd(value_by_fold_TTD)
    value[[t]][[2]][3] <- mean(obs_values_TTD)
    value[[t]][[2]][4] <- mean(max_values_TTD)
    value[[t]][[2]][5] <- mean(max_values_TTD)/sqrt(log(ntrt))
  }           # End loop over outcome variables
  
  ret <- sapply(c("value","importance","pred.cent.BAR","pred.cent.TTD"),function(x) NULL)
  ret$value <- value
  ret$importance <- importance
  ret$pred.cent.BAR <- pred.cent.BAR
  ret$pred.cent.TTD <- pred.cent.TTD
  return(ret)
}

