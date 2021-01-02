# collect.results.R

# collect results of PDX analyses from separate csv files into lists

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q week R CMD BATCH"  # code to submit job to cluster

setwd(dir1)
load("split.cm.data.rda")

library(stringi)

tumor_types = c("BRCA", "CM", "CRC", "NSCLC", "PDAC")
num_genes = c(50, 100, 500, 1000)
trains = c("BAR", "TTD")
evals = c("BAR", "TTD")

# collect results of analyses using screened genes
collect.results = function(method) {
  
  templist = list(BRCA = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                     "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                     "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                     "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
  CM = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
              "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
              "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
              "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
  CRC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
               "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
               "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
               "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
  NSCLC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                 "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                 "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                 "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
  PDAC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))))
  
  for (i in 1:length(tumor_types)) {
    
    ntrt = length(unique(split.cm.data[[i]]$Treatment))
    
    for (j in 1:length(num_genes)) {
      
      for (k in 1:length(trains)) {
        
        for (m in 1:length(evals)) {
          
          if (trains[k] == "BAR") outcome = "BAR"
          if (trains[k] == "TTD") outcome = "Surv"

          file = sprintf("%s_%s_%s.%d.biomarkers_%s.csv", tumor_types[i], outcome, tumor_types[i], num_genes[j], method)
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
          
          templist[[i]][[j]][[k]][[m]] = temp

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

          templist[[i]][[j]][[k]][[m]] = temp

          }

          rm(dat)
          
        }  # end loop through outcomes evaluated on
      }  # end through outcomes trained on
    }  # loop through number of genes 
  }  # loop through tumor types
  
  return(templist)  

}  # end collect.results function

methods = c("ql1", "ql2", "qlrf1", "qlrf2", "ql1smooth", "ql2smooth", "qlrf1smooth", "qlrf2smooth",  "owllinear", "owlkernel", "owllinearsmooth", "owlkernelsmooth", "sl4", "sl6", "sl8", "sl16")
for (i in 1:length(methods)) {
  results = collect.results(methods[i])
  print(results)
  assign(sprintf("value_%s", methods[i]), results)
  savecmd = sprintf("save(value_%s, file = 'value_%s.rda')", methods[i], methods[i])
  eval(parse(text = savecmd))
}

# collect results analyses using deep learning reduced data
collect.results.dl = function(method) {
  
  templist = list(BRCA = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                              "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
                  CM = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                            "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                            "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                            "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
                  CRC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                             "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                             "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                             "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
                  NSCLC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                               "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                               "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                               "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))),
                  PDAC = list("50" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)), 
                              "100" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "500" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL)),
                              "1000" = list(train_BAR = list(eval_BAR = NULL, eval_TTD = NULL), train_TTD = list(eval_BAR = NULL, eval_TTD = NULL))))
  
  for (i in 1:length(tumor_types)) {
    
    ntrt = length(unique(split.cm.data[[i]]$Treatment))
    
    for (j in 1:length(num_genes)) {
      
      for (k in 1:length(trains)) {
        
        for (m in 1:length(evals)) {
          
          if (trains[k] == "BAR") outcome = "BAR"
          if (trains[k] == "TTD") outcome = "Surv"
          
          # biomarkers.feat.new references deep learning reduced data
          file = sprintf("%s_%s_%s.%d.biomarkers.feat.new_%s.csv", tumor_types[i], outcome, tumor_types[i], num_genes[j], method)          
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
            
            templist[[i]][[j]][[k]][[m]] = temp
            
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
            
            templist[[i]][[j]][[k]][[m]] = temp
            
          }
          
          rm(dat)
          
        }  # end loop through outcomes evaluated on
      }  # end through outcomes trained on
    }  # loop through number of genes 
  }  # loop through tumor types
  
  return(templist)  
  
}  # end collect.results.dl function

methods = c("qldl1", "qldl2", "qlrfdl1", "qlrfdl2", "ql1smoothdl", "ql2smoothdl", "qlrf1smoothdl", "qlrf2smoothdl", "owllineardl", "owlkerneldl", "owllinearsmoothdl", "owlkernelsmoothdl")
for (i in 1:length(methods)) {
  results = collect.results.dl(methods[i])
  print(results)
  assign(sprintf("value_%s", methods[i]), results)
  savecmd = sprintf("save(value_%s, file = 'value_%s.rda')", methods[i], methods[i])
  eval(parse(text = savecmd))
}

# collect value results of RF/LASSO analyses using screened genes
collect.values.rf = function(method,tumorlist,ngeneslist) {
  
  value_rf <- sapply(tumorlist,function(x) NULL)
  for (q in 1:length(tumorlist)) {
    value_rf[[q]] <- sapply(ngeneslist,function(x) NULL)
  }
  
  for (q in 1:length(tumorlist)) {
    for (r in 1:length(ngeneslist)) {
      file=sprintf("%s/results_%s_%s_%s.rda",dir2,method,tumorlist[q],ngeneslist[r])
      load(file)
      value_rf[[q]][[r]] <- temp$value
    }
  }
  
  return(value_rf)  
  
}  # end collect.value.rf function

# collect importance results of RF/LASSO analyses using screened genes
collect.imps.rf = function(method,tumorlist,ngeneslist) {
  
  importance_rf <- sapply(tumorlist,function(x) NULL)
  for (q in 1:length(tumorlist)) {
    value_rf[[q]] <- sapply(ngeneslist,function(x) NULL)
  }
  
  for (q in 1:length(tumorlist)) {
    for (r in 1:length(ngeneslist)) {
      file=sprintf("%s/results_%s_%s_%s.rda",dir2,method,tumorlist[q],ngeneslist[r])
      load(file)
      value_rf[[q]][[r]] <- temp$importance
    }
  }
  
  return(importance_rf)  
  
}  # end collect.imps.rf function

# collect predicted values results of RF/LASSO analyses using screened genes
collect.predvals.rf = function(tumorlist,ngeneslist) {
  
  pred_vals_rf <- sapply(tumorlist,function(x) NULL)
  for (q in 1:length(tumorlist)) {
    pred_vals_rf[[q]] <- sapply(ngeneslist,function(x) NULL)
    for (r in 1:length(ngeneslist)) {
      pred_vals_rf[[q]][[r]] <- sapply(c("BAR","TTD"),function(x) NULL)
    }
  }
  
  for (q in 1:length(tumorlist)) {
    for (r in 1:length(ngeneslist)) {
      file=sprintf("%s/results_rf_%s_%s.rda",dir2,tumorlist[q],ngeneslist[r])
      load(file)
      pred_vals_rf[[q]][[r]][[1]] <- temp$pred.cent.BAR
      pred_vals_rf[[q]][[r]][[2]] <- temp$pred.cent.TTD
    }
  }
  
  return(pred_vals_rf)  
  
}  # end collect.predvals.rf function

tumorlist <- c("BRCA","CM","CRC","NSCLC","PDAC")
ngeneslist <- c("50","100","500","1000")
methods=c("rf","rf_dl","lasso","lasso_dl")
for (method in methods) {
  results <- collect.values.rf(method,tumorlist,ngeneslist)
  print(results)
  assign(sprintf('value_%s',method),results)
  savecmd = sprintf("save(value_%s, file = 'value_%s.rda')",method,method)
  eval(parse(text = savecmd))
  if (method=="rf") {
    pred_vals_rf <- collect.predvals.rf(tumorlist,ngeneslist)
    save(pred_vals_rf,file="pred_vals_rf.rda")
  }
  if (method=="rf" | method=="lasso") {
    imp <- collect.imps.rf(method,tumorlist,ngeneslist)
    assign(sprintf('importance_%s',method),imp)
    savecmd <- sprintf("save(importance_%s,file='importance_%s.rda')",method,method)
    eval(parse(text=savecmd))
  }
}
