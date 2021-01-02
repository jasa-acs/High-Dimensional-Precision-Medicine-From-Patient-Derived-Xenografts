# pdx.ql.owl.sl.pre.script.R

# script to submit jobs for all Q-learning and OWL analyses as preparation for super-learning

dir1 = "~/PDX/"  # main directory
dir2 = "/netscr/luckett/"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q day R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

load("trts.by.cancer.rda")

# Q-learning
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.sl.R", sep = "")
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
      pathtoRcode = paste(dir1, "pdx.ql2.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql1.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# random forest Q-learning
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.sl.R", sep = "")
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
      pathtoRcode = paste(dir1, "pdx.ql.rf2.sl.R", sep = "")
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
}

# Q-learning with RF smoothed outcomes
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql1smooth.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql2.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql2smooth.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}
#
# random forest Q-learning with random forest smoothed outcomes
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1smooth.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # RF q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2smooth.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# linear Q-learning with deep learning reduced data
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql1(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qldl1.csv')",
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_dlql1.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql2.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qldl2.csv')",
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_dlql2.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# random forest Q-learning with deep learning reduced data
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qlrfdl1.csv')",
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_dlqlrf1.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # RF q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf2.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qlrfdl2.csv')",
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_dlqlrf2.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# Q-learning with random forest smoothed outcomes and deep learning reduced data
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', strip = F, outstring = '_ql1smoothdl.csv')",
                     tumor_type, outcome, gene.data.file, nums[j],  dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql1smoothdl.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql2.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', strip = F, outstring = '_ql2smoothdl.csv')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_ql2smoothdl.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

# random forest Q-learning with random forest smoothed outcomes and deep learning reduced data
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]

      # RF q-learning without pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qlrf1smoothdl.csv')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf1smoothdl.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))

      # RF q-learning with pseudo-values
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.ql.rf1.smooth.sl.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.ql.rf2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qlrf2smoothdl.csv')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_qlrf2smoothdl.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
    }
  }
}

