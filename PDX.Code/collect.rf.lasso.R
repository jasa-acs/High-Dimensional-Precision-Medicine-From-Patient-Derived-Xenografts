library(stringi)

tumor_types = c("BRCA", "CM", "CRC", "NSCLC", "PDAC")
num_genes = c(50, 100, 500, 1000)

rflist = list(BRCA = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
                CM = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
                CRC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
                NSCLC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
                PDAC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL))


for (i in 1:length(tumor_types)) {
  
  for (j in 1:length(num_genes)) {
    
    tempname = sprintf("results_rf_rna_%s_%s.rda", tumor_types[i], num_genes[j])
    load(tempname)
    
    rflist[[i]][[j]] = temp$value
  }
}


lassolist = list(BRCA = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
              CM = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
              CRC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
              NSCLC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL),
              PDAC = list("50" = NULL, "100" = NULL, "500" = NULL, "1000" = NULL))


for (i in 1:length(tumor_types)) {
  
  for (j in 1:length(num_genes)) {
    
    tempname = sprintf("results_lasso_rna_%s_%s.rda", tumor_types[i], num_genes[j])
    load(tempname)
    
    lassolist[[i]][[j]] = temp$value
  }
}

value_rf = rflist
value_lasso = lassolist
save(value_rf, file = 'value_rf.rda')
save(value_lasso, file = 'value_lasso.rda')




