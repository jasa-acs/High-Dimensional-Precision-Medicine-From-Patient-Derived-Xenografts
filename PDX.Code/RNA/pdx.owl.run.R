# pdx.owl.run.R

# runs analysis using outcome weighted learning on PDX data

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
    gene.data.file = paste(tumor_type, ".", nums[j], ".biomarkers.csv", sep = "")
    for (k in 1:length(outs)) {
      outcome = outs[k]
      
      # q-learning without pseudo-values
      R = rep(0, 2)
      R[1] = "source('~/PDX/pdx.owl.R')"
      R[2] = sprintf("pdx.owl(cancer.type = '%s', outcome = '%s', gene.data.file = '%s', input_dir = '%s', output_dir = '%s')", 
                     tumor_type, outcome, gene.data.file, dir1, dir1)
      CMD = paste(R, collapse = "\n")
      pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_owl.R", tumor_type, outcome, nums[j])
      #pathtoRfile = sprintf("%s_%s_%s_submit_owl.R", tumor_type, outcome, nums[j])
      write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
      system(sprintf("bsub -M 10 %s -q week R CMD BATCH %s", dumpfile, pathtoRfile))
    }
  }
}






