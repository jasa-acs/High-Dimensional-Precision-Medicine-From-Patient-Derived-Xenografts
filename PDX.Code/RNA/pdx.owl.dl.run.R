# pdx.owl.dl.run.R

# all OWL DL analyses together

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
dumpfile = "-o /netscr/luckett/dump.txt"  # specify file here to suppress emails
setwd(dir1)

load("trts.by.cancer.rda")

# owl
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      c2string = "6:18"
      if (tumor_type == "CM") c2string = "6:12"
      # owl with linear decision function
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.dl.owl.linear.R')"
      R[2] = sprintf("pdx.dl.owl.linear(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1s = c(0, 1), c2s = %s)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c2string)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owllineardl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 10 %s -q week R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}

outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      c2string = "6:18"
      if (tumor_type == "CM") c2string = "6:12"
      # owl with Gaussian kernel
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.dl.owl.kernel.R')"
      R[2] = sprintf("pdx.dl.owl.kernel(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', c1s = c(0, 1), c2s = %s)",
                     tumor_type, outcome, gene.data.file, dir1, dir1, c2string)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owlkerneldl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 10 %s -q week R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}

outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      c2string = "6:18"
      if (tumor_type == "CM") c2string = "6:12"
      # owl with linear decision function
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.dl.owl.linear.smooth.R')"
      R[2] = sprintf("pdx.dl.owl.linear.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s',
                     numgenes = %s, input_dir = '%s', output_dir = '%s', c1s = c(0, 1), c2s = %s)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c2string)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owllinearsmoothdl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 -q week %s R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}

outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      c2string = "6:18"
      if (tumor_type == "CM") c2string = "6:12"
      # owl with Gaussian kernel
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.dl.owl.kernel.smooth.R')"
      R[2] = sprintf("pdx.dl.owl.kernel.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s',
                     numgenes = %s, input_dir = '%s', output_dir = '%s', c1s = c(0, 1), c2s=  %s)",
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1, c2string)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owlkernelsmoothdl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 -q week %s R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}

