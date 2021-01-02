# pdx.ql.owl.final.script.R

# script to submit all final Q-learning and OWL analyses on full data to save coefficients and variable importance scores

dir1 = "~/PDX/"  # main directory
dir2 = "/netscr/luckett/"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q day R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

load("trts.by.cancer.rda")

# Q-learning
outs = c("Surv", "BAR")
nums = c("50", "100", "500", "1000")
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%s.biomarkers_ql1.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql1.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %s, c2 = %s)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql1.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%s.biomarkers_ql2.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql2.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %s, c2 = %s)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql2.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# random forest Q-learning
outs = c("Surv", "BAR")
nums = c("50", "100", "500", "1000")
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%s.biomarkers_qlrf1.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %f, c2 = %f)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%s.biomarkers_qlrf2.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # RF q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf2.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %f, c2 = %f)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# Q-learning with random forest smoothed outcomes
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_ql1smooth.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.smooth.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql1.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql1.smooth.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_ql2smooth.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql2.smooth.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql2.smooth.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# random forest Q-learning with random forest smoothed outcomes
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_qlrf1smooth.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1.smooth.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # read resulting tuning parameters
      file = sprintf("%s_%s_%s.%d.biomarkers_qlrf2smooth.csv", tumor_type, outcome, tumor_type, nums[j])
      params = read.csv(file)
      c1 =  params$c1; c2 = params$c2

      # RF q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf2.smooth.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2.smooth.final.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

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
      pathtoRcode = paste(dir1, "pdx.owl.linear.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.owl.linear.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_owllinearfinal.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# linear OWL with random forest smoothed outcomes
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)
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
      pathtoRcode = paste(dir1, "pdx.owl.linear.smooth.final.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.owl.linear.smooth.final(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', c1 = %d, c2 = %d)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c1, c2)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_owllinearsmoothfinal.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

