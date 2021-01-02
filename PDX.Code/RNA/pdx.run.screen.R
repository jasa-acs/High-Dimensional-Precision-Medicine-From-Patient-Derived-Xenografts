# pdx.run.screen

# runs all gene screening for PDX analyses

dir1 = "~/PDX/rna_only"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
dumpfile = "-o /netscr/luckett/dump.txt"  # specify file here to suppress emails
setwd(dir1)

source("pdx.import.fn.R")
import.pdx()
load("trts.by.cancer.rda")

# supervised screening for prognostic biomarkers
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  for (j in 1:length(trts.by.cancer[[i]])) {
    trt_type = trts.by.cancer[[i]][j]
    R = rep(0, 2)
    R[1] = "source('~/PDX/pdx.bdc1.R')"
    R[2] = sprintf("pdx.bdc1(tumor_type = '%s', trt_type = '%s', input_dir = '%s', output_dir = '%s')", tumor_type, trt_type, dir1, dir2)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("/netscr/luckett/%s_%s_submit_pre_rna.R", tumor_type, trt_type)
    #pathtoRfile = sprintf("%s_%s_submit_pre.R", tumor_type, trt_type)
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
  }
}

# supervised screening for predictive biomarkers
for (i in 1:length(trts.by.cancer)) {
  tumor_type = names(trts.by.cancer)[i]
  trts = trts.by.cancer[[i]]
  # all combinations of two trt's
  trt.comb = combn(trts, 2)
  for (j in 1:dim(trt.comb)[2]) {
    trt_type1 = trt.comb[1, j]
    trt_type2 = trt.comb[2, j]
    R = rep(0, 2)
    R[1] = "source('~/PDX/pdx.bdc2.R')"
    R[2] = sprintf("pdx.bdc2(tumor_type = '%s', trt_type1 = '%s', trt_type2 = '%s', input_dir = '%s', output_dir = '%s')", tumor_type, trt_type1, trt_type2, dir1, dir2)
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("/netscr/luckett/%s_%s_%s_submit_pre_rna.R", tumor_type, trt_type1, trt_type2)
    #pathtoRfile = sprintf("%s_%s_%s_submit_pre.R", tumor_type, trt_type1, trt_type2)
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("bsub -M 4 %s R CMD BATCH %s", dumpfile, pathtoRfile))
  }
}
