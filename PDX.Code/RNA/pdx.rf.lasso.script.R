# pdx.rf.lasso.script.R

# script to submit jobs for all standard (RF and LASSO) analyses

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q week R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

load("trts.by.cancer.rda")

tumorlist <- names(trts.by.cancer)
ngeneslist <- c("50","100","500","1000")

##################### MAIN ANALYSIS ########################

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.rf.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = rf.value(%s,clin,verbose=F,rna_only=F,deeplearn=F,fulldata=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_rf_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_rf.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.lasso.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = lasso.value(%s,clin,verbose=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_lasso_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_lasso.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

##################### DEEP LEARNING ########################

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.feat.new.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.rf.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%s]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = rf.value(%s,clin,verbose=F,rna_only=F,deeplearn=T,fulldata=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_rf_dl_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_rf_dl.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.feat.new.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.lasso.dl.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%s]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = lasso.value.DL(%s,clin,verbose=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_lasso_dl_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_lasso_dl.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

##################### RNA ONLY ########################

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.rf.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = rf.value(%s,clin,verbose=F,rna_only=T,deeplearn=F,fulldata=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_rf_rna_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_rf_rna.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.lasso.rna.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = lasso.value.rna(%s,clin,verbose=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_lasso_rna_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_lasso_rna.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

##################### RNA ONLY, DEEP LEARNING ########################

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.rna.feat.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.rf.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = rf.value(%s,clin,verbose=F,rna_only=T,deeplearn=T,fulldata=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_rf_rna_dl_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_rf_rna_dl.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

for (q in 1:length(tumorlist)) {
  for (r in 1:length(ngeneslist)) {
    infile <- sprintf("%s/%s.%s.biomarkers.rna.feat.csv",dir1,tumorlist[q],ngeneslist[r])
    R <- rep(0,5)
    R[1] <- sprintf("source('%s/pdx.lasso.dl.R')",dir1)
    R[2] <- sprintf("load('%s/split.cm.data.rda')",dir1)
    R[3] <- sprintf("clin = split.cm.data[[%d]][,c('logSurvScaled','RespScaled','Treatment','Model')]",q)
    R[4] <- sprintf("temp = lasso.value.DL(%s,clin,verbose=F)",infile)
    R[5] <- sprintf("save(temp,file='%s/results_lasso_rna_dl_%s_%s.rda')",dir1,tumorlist[q],ngeneslist[r])
    CMD = paste(R, collapse = "\n")
    pathtoRfile = sprintf("%s/%s_%s_submit_lasso_rna_dl.R",dir2,tumorlist[q], ngeneslist[r])
    write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
    system(sprintf("%s %s", tosubmit, pathtoRfile))
  }
}

##################### FULL DATA ########################


for (r in 1:length(ngeneslist)) {
  infile <- sprintf("%s/overall.%s.biomarkers.csv",dir1,ngeneslist[r])
  R <- rep(0,8)
  R[1] <- sprintf("source('%s/pdx.rf.R')",dir1)
  R[2] <- sprintf("source('%s/create.full.data.R')",dir1)
  R[3] <- sprintf("load('%s/split.cm.data.rda')",dir1)
  R[4] <- sprintf("load('%s/trts.by.cancer.rda')",dir1)
  R[5] <- sprintf("load('%s/genes.by.cancer.rda')",dir1)
  R[6] <- "clin = create.full.data(split.cm.data,trts.by.cancer,genes.by.cancer,use5=T,use4=F,no.genes=T)"
  R[7] <- sprintf("temp = rf.value(%s,clin,verbose=F,rna_only=F,deeplearn=F,fulldata=T)",infile)
  R[8] <- sprintf("save(temp,file='%s/results_rf_overall_%s.rda')",dir1,ngeneslist[r])
  CMD = paste(R, collapse = "\n")
  pathtoRfile = sprintf("%s/overall_%s_submit_rf.R",dir2,ngeneslist[r])
  write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
  system(sprintf("%s %s", tosubmit, pathtoRfile))
}

for (r in 1:length(ngeneslist)) {
  infile <- sprintf("%s/overall.%s.biomarkers.csv",dir1,ngeneslist[r])
  R <- rep(0,8)
  R[1] <- sprintf("source('%s/pdx.lasso.R')",dir1)
  R[2] <- sprintf("source('%s/create.full.data.R')",dir1)
  R[3] <- sprintf("load('%s/split.cm.data.rda')",dir1)
  R[4] <- sprintf("load('%s/trts.by.cancer.rda')",dir1)
  R[5] <- sprintf("load('%s/genes.by.cancer.rda')",dir1)
  R[6] <- "clin = create.full.data(split.cm.data,trts.by.cancer,genes.by.cancer,use5=T,use4=F,no.genes=T)"
  R[7] <- sprintf("temp = lasso.value(%s,clin,verbose=F)",infile)
  R[8] <- sprintf("save(temp,file='%s/results_lasso_overall_%s.rda')",dir1,ngeneslist[r])
  CMD = paste(R, collapse = "\n")
  pathtoRfile = sprintf("%s/overall_%s_submit_lasso.R",dir2,ngeneslist[r])
  write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
  system(sprintf("%s %s", tosubmit, pathtoRfile))
}

