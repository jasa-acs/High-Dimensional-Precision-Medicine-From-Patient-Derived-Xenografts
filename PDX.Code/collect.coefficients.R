# collect.coefficients.R

# collect coefficients and variable importance scores of PDX analyses from csv files into lists

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q week R CMD BATCH"  # code to submit job to cluster

setwd(dir1)
load("split.cm.data.rda")

tumor_types = c("BRCA", "CM", "CRC", "NSCLC", "PDAC")
num_genes = c(50, 100, 500, 1000)
trains = c("BAR", "TTD")

collect.coefficients = function(method) {
  
  templist = list(BRCA = list("50" = list(train_BAR = NULL, train_TTD = NULL), "100" = list(train_BAR = NULL, train_TTD = NULL), 
                              "500" = list(train_BAR = NULL, train_TTD = NULL), "1000" = list(train_BAR = NULL, train_TTD = NULL)), 
                  CM = list("50" = list(train_BAR = NULL, train_TTD = NULL), "100" = list(train_BAR = NULL, train_TTD = NULL), 
                            "500" = list(train_BAR = NULL, train_TTD = NULL), "1000" = list(train_BAR = NULL, train_TTD = NULL)),
                  CRC = list("50" = list(train_BAR = NULL, train_TTD = NULL), "100" = list(train_BAR = NULL, train_TTD = NULL), 
                             "500" = list(train_BAR = NULL, train_TTD = NULL), "1000" = list(train_BAR = NULL, train_TTD = NULL)),
                  NSCLC = list("50" = list(train_BAR = NULL, train_TTD = NULL), "100" = list(train_BAR = NULL, train_TTD = NULL), 
                               "500" = list(train_BAR = NULL, train_TTD = NULL), "1000" = list(train_BAR = NULL, train_TTD = NULL)),
                  PDAC = list("50" = list(train_BAR = NULL, train_TTD = NULL), "100" = list(train_BAR = NULL, train_TTD = NULL), 
                              "500" = list(train_BAR = NULL, train_TTD = NULL), "1000" = list(train_BAR = NULL, train_TTD = NULL)))
  
  for (i in 1:length(tumor_types)) {
    
    for (j in 1:length(num_genes)) {
      
      for (k in 1:length(trains)) {
        
        if (trains[k] == "BAR") outcome = "BAR"
        if (trains[k] == "TTD") outcome = "Surv"
        
        file = sprintf("%s_%s_%s.%d.biomarkers_%s.coefs.csv", tumor_types[i], outcome, tumor_types[i], num_genes[j], method)
        dat = read.csv(file)
        
        templist[[i]][[j]][[k]] = dat
        
      }  # end through outcomes trained on
    }  # loop through number of genes 
  }  # loop through tumor types
  
  return(templist)
  
}  # end collect.coefficients function

methods = c("ql1", "ql2", "ql1smooth", "ql2smooth", "qlrf1", "qlrf2", "qlrf1smooth", "qlrf2smooth", "owllinear", "owllinearsmooth")
for (i in 1:length(methods)) {
  results = collect.coefficients(methods[i])
  assign(sprintf("importance_%s", methods[i]), results)
  savecmd = sprintf("save(importance_%s, file = 'importance_%s.rda')", methods[i], methods[i])
  eval(parse(text = savecmd))
}


