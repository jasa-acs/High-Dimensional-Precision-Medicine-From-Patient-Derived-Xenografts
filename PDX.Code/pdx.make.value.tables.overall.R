# pdx.make.value.tables.overall.R

# takes the value objects from combined cancer analyses and turns them into usable tables

dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
tosubmit = "bsub -M 4 -o /netscr/luckett/dump.txt -q day R CMD BATCH"  # code to submit job to cluster

setwd(dir1)

value_table_R <- function(datadir,outdir,method) {
  R <- rep(0,25)
  R[1] <- sprintf("load('%s/value_%s.rda')",datadir,method)
  
  R[2] <- sprintf("tumornames <- names(value_%s)",method)
  R[3] <- sprintf("ngenelist <- names(value_%s[[1]])",method)
  R[4] <- "ng <- length(ngenelist)"
  R[5] <- sprintf("value_%s_tables <- sapply(tumornames,function(x) NULL)",method)
  
  R[6] <- "for (i in 1:length(tumornames)) {"
  R[7] <-   sprintf("value_%s_tables[[i]] <- matrix(NA,ncol=8,nrow=(ng*4+1))",method)
  R[8] <-   sprintf("value_%s_tables[[i]][1,] <- c('Lsup','Trained On','Evaluated On','Mean Value','SD Value','Obs. Value','Optimal Value','Adj. Opt. Value')",method)
  R[9] <-   "for (j in 1:length(ngenelist)) {"
  R[10] <-     "lsuprows <- (2+(ng*(j-1))):(1+ng*j)"
  R[11] <-     sprintf("value_%s_tables[[i]][lsuprows,1] <- c(ngenelist[j],rep('',ng-1))",method)
  R[12] <-     sprintf("value_%s_tables[[i]][lsuprows,2] <- c('BAR','','TTD','')",method)
  R[13] <-     sprintf("value_%s_tables[[i]][lsuprows,3] <- rep(c('BAR','TTD'),2)",method)
  R[14] <-   "}"
  R[15] <-   sprintf("nums <- 4:ncol(value_%s_tables[[i]])",method)
  R[16] <-   "for (j in 1:length(ngenelist)) {"
  R[17] <-     "for (k in 1:2) {"
  R[18] <-       "for (l in 1:2) {"
  R[19] <-        "row <- 1+(j-1)*ng+(k-1)*2+l"
  R[20] <-        sprintf("value_%s_tables[[i]][row,nums] <- value_%s[[i]][[j]][[k]][[l]][1:5]",method,method)
  R[21] <-      "}"
  R[22] <-     "}"
  R[23] <-  "}"
  R[24] <- "}"
  
  R[25] <- sprintf("save(value_%s_tables,file='%s/value_%s_tables.rda')",method,outdir,method)
  
  return(paste(R,collapse = "\n"))
}

methods = c("ql1_full", "ql2_full", "owllinear_full", "owlkernel_full", "qlrf1_full", "qlrf2_full", "owllinearsmooth_full", "owlkernelsmooth_full","rf_full","lasso_full")
for (method in methods) {
  CMD <- value_table_R(datadir=dir1,outdir=dir1,method=method)
  pathtoRfile = sprintf("/netscr/luckett/%s_make_value_table.R", method)
  write.table(CMD, file = pathtoRfile, quote = F, col.names = F, row.names = F)
  system(sprintf("%s %s", tosubmit, pathtoRfile))
}
