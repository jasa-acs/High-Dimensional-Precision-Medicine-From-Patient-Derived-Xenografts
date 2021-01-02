# pdx.ql.owl.sl.script.R

# script to submit jobs for all super-learning training

dir1 = "~/PDX/"  # main directory
dir2 = "/netscr/luckett/"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q day R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

load("trts.by.cancer.rda")


# super-learning 4
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.super.learning.4.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.super.learning.4(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_sl4.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
      
    }
  }
}

# super-learning 6
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.super.learning.6.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.super.learning.6(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_sl6.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
      
    }
  }
}

# super-learning 8
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.super.learning.8.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.super.learning.8(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_sl8.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
      
    }
  }
}

# super-learning 16
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file1 = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    gene.data.file2 = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      R = rep(0, 2)
      pathtoRcode = paste(dir1, "pdx.super.learning.16.R", sep = "")
      R[1] = sprintf("source('%s')", pathtoRcode)
      R[2] = sprintf("pdx.super.learning.16(cancer.type = '%s', outcome = '%s', gene.data.file1 = '%s', gene.data.file2 = '%s', numgenes = '%s', input_dir = '%s', output_dir = '%s')",
                     tumor_type, outcome, gene.data.file1, gene.data.file2, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("%s_%s_%s_submit_sl16.R", tumor_type, outcome, nums[j])
      pathtoRfile = paste(dir2, pathtoRfile, sep = "")
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("%s %s", tosubmit, pathtoRfile))
      
    }
  }
}
