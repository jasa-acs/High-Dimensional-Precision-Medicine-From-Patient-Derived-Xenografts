# pdx.owl.kernel.run.R

# runs analysis using outcome weighted learning with Gaussian kernel on PDX data

dir1 = "/mnt/home/yzhan234/PDX/PDX.CodeAndData/RNA/"
dir2 = "/mnt/home/yzhan234/PDX/PDX.CodeAndData/RNA/"  # scratch directory
tosubmit = "bwsubmit r"  # code to submit job to cluster
setwd(dir1)

load("trts.by.cancer.rda")

# owl
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.rna.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      c2string = "6:18"
      if (tumor_type == "CM") c2string = "6:12"
      # owl with kernel decision function
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.owl.kernel.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.owl.kernel(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1s = c(0, 1), c2s = %s, outstring = '_owlkernel_rna.csv')",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c2string)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_owlkernel_rna.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}






