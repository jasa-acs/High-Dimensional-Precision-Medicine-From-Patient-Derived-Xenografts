# pdx.run.ql.rfsmooth.R

# runs analysis using Q-learning on PDX data

dir1 = "/mnt/home/yzhan234/PDX/PDX.CodeAndData/RNA/"
dir2 = "/mnt/home/yzhan234/PDX/PDX.CodeAndData/RNA/"  # scratch directory
tosubmit = "bwsubmit r"  # code to submit job to cluster
setwd(dir1)

load("trts.by.cancer.rda")

# RF q-learning
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.rna.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1smooth_rna.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
      
      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf2.smooth.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2smooth_rna.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}
