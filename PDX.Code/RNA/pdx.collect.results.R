# pdx.collect.results.R

# collect results and produce tables

library(xtable)

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
setwd(dir1)

load("trts.by.cancer.rda")

# collect function 
# takes as argument tumor_type and outcome
# reports results for maximizing outcome within tumor_type using num_genes genes
collect.results = function(tumor_type, outcome, num) {
  
  # matrix to store results
  #res = matrix(NA, ncol = 7, nrow = 7)
  res = matrix(NA, ncol = 7, nrow = 5)
  res[1, ] = c("Method", "Outcome", "(c_1, c_2)", "Untreated", "Observed", "Value", "Optimal")  
  res[2:nrow(res), 1] = c("QL", "", "QL (pseudo)", "")
  res[2:nrow(res), 2] = c("BAR", "Surv", "BAR", "Surv")
  res[2:nrow(res), 4] = rep("0", nrow(res) - 1)
  
  file1 = paste(tumor_type, "_", outcome, "_", tumor_type, ".", num, ".", "biomarkers", "_ql1.csv", sep = "")
  file2 = paste(tumor_type, "_", outcome, "_", tumor_type, ".", num, ".", "biomarkers", "_ql2.csv", sep = "")
  #file3 = paste(tumor_type, "_", outcome, "_", tumor_type, ".", num, ".", "biomarkers", "_owl.csv", sep = "")
  ql1 = read.csv(file1)
  ql2 = read.csv(file2)
  #owl = read.csv(file3)
  
  # q-learning without pseudo-values results
  res[2 ,3] = sprintf("(%d, %d)", ql1$c1, ql1$c2)
  res[2, 5:7] = c(sprintf("%.4f", round(ql1$observed.resp, 4)), sprintf("%.4f", round(ql1$mean.resp, 4)), sprintf("%.4f", round(ql1$optimal.resp, 4)))
  res[3, 5:7] = c(sprintf("%.4f", round(ql1$observed.surv, 4)), sprintf("%.4f", round(ql1$mean.survival, 4)), sprintf("%.4f", round(ql1$optimal.surv, 4)))
  
  # q-learning with pseudo-values results
  res[4 ,3] = sprintf("(%d, %d)", ql2$c1, ql2$c2)
  res[4, 5:7] = c(sprintf("%.4f", round(ql2$observed.resp, 4)), sprintf("%.4f", round(ql2$mean.resp, 4)), sprintf("%.4f", round(ql2$optimal.resp, 4)))
  res[5, 5:7] = c(sprintf("%.4f", round(ql2$observed.surv, 4)), sprintf("%.4f", round(ql2$mean.survival, 4)), sprintf("%.4f", round(ql2$optimal.surv, 4)))
  
  # owl
  #res[6 ,3] = sprintf("(%d, %d)", owl$c1, owl$c2)
  #res[6, 5:7] = c(sprintf("%.2f", round(owl$observed.resp, 2)), sprintf("%.2f", round(owl$mean.resp, 2)), sprintf("%.2f", round(owl$optimal.resp, 2)))
  #res[7, 5:7] = c(sprintf("%.2f", round(owl$observed.surv, 2)), sprintf("%.2f", round(owl$mean.survival, 2)), sprintf("%.2f", round(owl$optimal.surv, 2)))
  
  if (outcome == "Surv") out.capt = "Log time to double"
  if (outcome == "BAR") out.capt = "Best average response"
  caption = sprintf("Results when maximizing %s in tumor type %s using %d genes", out.capt, tumor_type, num)
  label = sprintf("tab_%s_%s_%d", tumor_type, outcome, num)
  rownames(res) = NULL
  
  return(xtable(res, caption = caption, label = label))
  
}

outcomes = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)

sink('pdx.tables.txt')

for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(outcomes)) {
    outcome = outcomes[j]
    for (k in 1:length(nums)) {
      num = nums[k]
      print.xtable(collect.results(tumor_type, outcome, num), include.rownames = F)
    }
  }
}

sink()

