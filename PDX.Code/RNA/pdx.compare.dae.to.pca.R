# pdx.compare.dae.to.pca

# Directory where RF code is stored
codedir <- "~/school/research/pdx/code"
# Directory where data object split.cm.data and biomarker datasets are stored
datadir <- "~/school/research/pdx/data"
# Directory where final data objects should be sent
outdir <- "~/school/research/pdx/results"

DAE_errors <- matrix(nrow=length(tumorlist),ncol=length(ngeneslist))
rownames(DAE_errors) <- tumorlist
colnames(DAE_errors) <- ngeneslist
PCA_errors <- matrix(nrow=length(tumorlist),ncol=length(ngeneslist))
rownames(PCA_errors) <- tumorlist
colnames(PCA_errors) <- ngeneslist
for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    summ <- read.csv(sprintf("%s/%s.%s.biomarkers.rna.summary.csv",outdir,tumorlist[q],ngeneslist[r]))
    chosen.row <- which(summ$AE_error==min(summ$AE_error))
    DAE_errors[q,r] <- summ[chosen.row,"AE_error"]
    PCA_errors[q,r] <- summ[chosen.row,"PCA_error"]
  }
}

save(DAE_errors,file=sprintf("%s/DAE_errors.rda",outdir))
save(PCA_errors,file=sprintf("%s/PCA_errors.rda",outdir))
