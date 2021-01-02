# collect.results.overall.R

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q week R CMD BATCH"  # code to submit job to cluster


setwd(dir1)

load("full.data.rda")

num_genes = c(50, 100, 500, 1000)
trains = c("BAR", "TTD")
evals = c("BAR", "TTD")

collect.results.overall = function(method) {
  
  templist = list(overall = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                              "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))))
    
    ntrt = length(unique(dat$Treatment))
    tumor_type = "overall"
    
    for (j in 1:length(num_genes)) {
      
      for (k in 1:length(trains)) {
        
        for (m in 1:length(evals)) {
          
          if (trains[k] == "BAR") outcome = "BAR"
          if (trains[k] == "TTD") outcome = "Surv"
          
          file = sprintf("%s_%s_%s.%d.biomarkers_%s.csv", tumor_type, outcome, tumor_type, num_genes[j], method)
          #dat = read.csv(file)
          dat = NA
          try(expr = (dat = read.csv(file)))
          if (!is.na(dat)) {
            if (evals[m] == "BAR") {
              temp = c(dat$mean.response, sqrt(dat$var.response), dat$observed.resp, dat$optimal.resp, dat$optimal.resp / sqrt(log(ntrt - dat$c1 - 1)), dat$c1, dat$c2)
              names(temp) = c("mean_value", "sd_value", "obs_value", "optimal_value", "optimal_adj_value", "c1", "c2")
            }
            if (evals[m] == "TTD") {
              temp = c(dat$mean.survival, sqrt(dat$var.survival), dat$observed.surv, dat$optimal.surv, dat$optimal.surv / sqrt(log(ntrt - dat$c1 - 1)), dat$c1, dat$c2)
              names(temp) = c("mean_value", "sd_value", "obs_value", "optimal_value", "optimal_adj_value", "c1", "c2")
            }
            
            templist[[1]][[j]][[k]][[m]] = temp
          }
          
          if (is.na(dat)) {
            if (evals[m] == "BAR") {
              temp = rep(NA, 7)
              names(temp) = c("mean_value", "sd_value", "obs_value", "optimal_value", "optimal_adj_value", "c1", "c2")
            }
            if (evals[m] == "TTD") {
              temp = rep(NA, 7)
              names(temp) = c("mean_value", "sd_value", "obs_value", "optimal_value", "optimal_adj_value", "c1", "c2")
            }
            
            templist[[1]][[j]][[k]][[m]] = temp
            
          }
          
          rm(dat)
          
          
        }  # end loop through outcomes evaluated on
      }  # end through outcomes trained on
    }  # loop through number of genes 
  
  return(templist)
  
}  # end collect function

methods = c("ql1", "ql2", "owllinear", "owlkernel", "qlrf1", "qlrf2", "owllinearsmooth", "owlkernelsmooth")
for (i in 1:length(methods)) {
  results = collect.results.overall(methods[i])
  assign(sprintf("value_%s_full", methods[i]), results)
  savecmd = sprintf("save(value_%s_full, file = 'value_%s_full.rda')", methods[i], methods[i])
  eval(parse(text = savecmd))
}

# collect value results of RF/LASSO analyses using screened genes
collect.values.overall.rf = function(method,ngeneslist) {
  
  value_rf <- sapply(ngeneslist,function(x) NULL)
  
  for (r in 1:length(ngeneslist)) {
    file=sprintf("%s/results_%s_overall_%s.rda",dir2,method,ngeneslist[r])
    load(file)
    value_rf[[r]] <- temp$value
  }
  
  return(value_rf)  
  
}  # end collect.value.rf function

# collect importance results of RF/LASSO analyses using screened genes
collect.imps.overall.rf = function(method,ngeneslist) {
  
  importance_rf <- sapply(ngeneslist,function(x) NULL)
  
  for (r in 1:length(ngeneslist)) {
    file=sprintf("%s/results_%s_overall_%s.rda",dir2,method,ngeneslist[r])
    load(file)
    value_rf[[r]] <- temp$importance
  }
  
  return(importance_rf)  
  
}  # end collect.imps.rf function

# collect predicted values results of RF/LASSO analyses using screened genes
collect.predvals.overall.rf = function(ngeneslist) {
  
  pred_vals_rf <- sapply(ngeneslist,function(x) NULL)
  for (r in 1:length(ngeneslist)) {
    pred_vals_rf[[r]] <- sapply(c("BAR","TTD"),function(x) NULL)
  }
  
    for (r in 1:length(ngeneslist)) {
      file=sprintf("%s/results_rf_overall_%s.rda",dir2,ngeneslist[r])
      load(file)
      pred_vals_rf[[r]][[1]] <- temp$pred.cent.BAR
      pred_vals_rf[[r]][[2]] <- temp$pred.cent.TTD
    }
  
  
  return(pred_vals_rf)  
  
}  # end collect.predvals.rf function

ngeneslist <- c("50","100","500","1000")
methods=c("rf","lasso")
for (method in methods) {
  results <- collect.values.overall.rf(method,tumorlist,ngeneslist)
  print(results)
  assign(sprintf('value_%s_full',method),results)
  savecmd = sprintf("save(value_%s_full, file = 'value_%s_full.rda')",method,method)
  eval(parse(text = savecmd))
  if (method=="rf") {
    pred_vals_rf_full <- collect.predvals.overall.rf(tumorlist,ngeneslist)
    save(pred_vals_rf_full,file="pred_vals_rf_full.rda")
  }
  imp <- collect.imp.overall.rf(method,tumorlist,ngeneslist)
  assign(sprintf('importance_%s_full',method),imp)
  savecmd <- sprintf("save(importance_%s_full,file='importance_%s_full.rda')",method,method)
  eval(parse(text=savecmd))
}

