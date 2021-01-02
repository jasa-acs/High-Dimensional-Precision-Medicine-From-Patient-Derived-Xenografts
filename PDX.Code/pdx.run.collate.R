# pdx.run.collate.R

# collates screening results for PDX analyses

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q week R CMD BATCH"  # code to submit job to cluster
setwd(dir1)

load("trts.by.cancer.rda")

# collate screening results
nums = 'c(50, 100, 500, 1000)'  # number of genes to keep
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  R = rep(0, 2)
  R[1] = "source('~/PDX/pdx.collate.R')"
  R[2] = sprintf("pdx.collate(tumor_type = '%s', nums = %s, input_dir = '%s', output_dir = '%s')", tumor_type, nums, dir2, dir1)
  CMD = paste(R, collapse = "\n")
  pathtoRfile = sprintf("/netscr/luckett/%s_submit_collate.R", tumor_type)
  write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
  system(sprintf("%s %s", tosubmit, pathtoRfile))
}
