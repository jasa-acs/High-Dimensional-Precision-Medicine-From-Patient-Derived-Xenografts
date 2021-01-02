# pdx.run.ql.dl.R

# runs analysis using Q-learning on PDX data

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
      
#      if (tumor_type == "CRC") {

      # q-learning without pseudo-values
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.ql1.R')"
      R[2] = sprintf("pdx.ql1(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qldl1.csv')", 
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_dlql1.R", tumor_type, outcome, nums[j])
      #pathtoRfile = sprintf("%s_%s_%s_submit_ql1.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
      
      # q-learning with pseudo-values
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.ql2.R')"
      R[2] = sprintf("pdx.ql2(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s', strip = F, outstring = '_qldl2.csv')", 
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_dlql2.R", tumor_type, outcome, nums[j])
      #pathtoRfile = sprintf("%s_%s_%s_submit_ql2.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
      #}
    }
  }
}
