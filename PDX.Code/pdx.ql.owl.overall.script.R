# pdx.ql.owl.overall.script.R

# runs analysis using Q-learning cross-cancer data

dir1 = "~/PDX/"  # main directory
dir2 = "/netscr/luckett/"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q day R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)
tumor_type = "overall"

# Q-learning
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]

    # q-learning without pseudo-values
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.ql1.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.ql1(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_ql1.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))

    # q-learning with pseudo-values
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.ql2.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.ql2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_ql2.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}


# random forest Q-learning
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]

    # RF q-learning without pseudo-values
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.ql.rf1.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.ql.rf1(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))

    # RF q-learning with pseudo-values
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.ql.rf2.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.ql.rf2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}


# linear OWL
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]

    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.owl.linear.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.owl.linear(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1s = c(0), c2s = NA)",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_owllinear.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

# Gaussian OWL
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]

    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.owl.kernel.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.owl.kernel(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1s = c(0), c2s = NA,)",
                   tumor_type, outcome, gene.data.file, dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_owlkernel.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}


# linear OWL with random forest smoothed outcomes
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]
      
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.owl.linear.smooth.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.owl.linear.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1s = c(0), c2s = NA)", 
                   tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_owllinearsmooth.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

# Gaussian OWL with random forest smoothed outcomes
for (j in 1:length(nums)) {
  gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
  for (k in 1:length(outs)) {
    outcome = outs[k]
    
    R = rep(0, 2)
    pathtoRcode = paste(dir1, "pdx.owl.kernel.smooth.R", sep = "")
    R[1] = sprintf("source('%s')", pathtoRcode)
    R[2] = sprintf("pdx.owl.kernel.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1s = c(0), c2s = NA,)", 
                    tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s_%s_%s_submit_owlkernelsmooth.R", tumor_type, outcome, nums[j])
    pathtoRfile = paste(dir2, pathtoRfile, sep = "")
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}
