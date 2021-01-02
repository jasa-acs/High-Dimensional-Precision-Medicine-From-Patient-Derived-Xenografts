# pdx.run.qlsmooth.R

# runs analysis using Q-learning on PDX data with deep learning reduced data

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
dumpfile = "-o /netscr/luckett/dump.txt"  # specify file here to suppress emails
setwd(dir1)

load("trts.by.cancer.rda")

# q-learning
outs = c("Surv", "BAR")
nums = c(50, 100, 500, 1000)  # number of genes to keep -- needs to be a vector so we can index
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(nums)) {
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.feat.new.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      # q-learning without pseudo-values
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.ql1.smooth.R')"
      R[2] = sprintf("pdx.ql1.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', outstring = '_ql1smoothdl.csv')", 
                     tumor_type, outcome, gene.data.file, nums[j],  dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_ql1smoothdl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
      
      # q-learning with pseudo-values
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.ql2.smooth.R')"
      R[2] = sprintf("pdx.ql2.smooth(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', numgenes = %s, input_dir = '%s', output_dir = '%s', outstring = '_ql2smoothdl.csv')", 
                     tumor_type, outcome, gene.data.file, nums[j], dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_ql2smoothdl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}
