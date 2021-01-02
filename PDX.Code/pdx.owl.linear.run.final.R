# pdx.owl.linear.run.final.R

# repeats linear OWL on full data to estimate coefficients without cross-validation

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
dumpfile = "-o /netscr/luckett/dump.txt"  # specify file here to suppress emails
setwd(dir1)

load("trts.by.cancer.rda")

# linear OWL
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_owllinear.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2
      
      # owl with linear decision function
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.owl.linear.final.R')"
      R[2] = sprintf("pdx.owl.linear.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)", 
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owllinearfinal.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}

# linear OWL with random forest smoothed outcomes
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_owllinearsmooth.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2
      
      # owl with linear decision function
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.owl.linear.smooth.final.R')"
      R[2] = sprintf("pdx.owl.linear.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)", 
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owllinearsmoothfinal.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}




