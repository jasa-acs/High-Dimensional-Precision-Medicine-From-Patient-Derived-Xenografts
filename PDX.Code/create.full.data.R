create.full.data <- function(split.cm.data,trts.by.cancer,genes.by.cancer,use4=FALSE,use5=TRUE,no.genes=FALSE) {
  if (!no.genes) {
    if (use5) {
      intersect.trts <- trts.by.cancer[[1]]
      for (i in 2:5) {
        intersect.trts <- intersect(intersect.trts,trts.by.cancer[[i]])
      }
      intersect.genes <- genes.by.cancer[[1]]
      for (i in 2:5) {
        intersect.genes <- intersect(intersect.genes,genes.by.cancer[[i]])
      }
      use.types <- names(split.cm.data)
    }
    if (use4) {
      int.trts.list <- sapply(names(trts.by.cancer),function(x) NULL)
      for (i in 1:5) {
        ind <- (1:5)[-i]
        int.trts.list[[i]] <- trts.by.cancer[[ind[1]]]
        for (j in 2:4) {
          int.trts.list[[i]] <- intersect(int.trts.list[[i]],trts.by.cancer[[ind[j]]])
        }
      }
      type.to.remove <- which(sapply(int.trts.list,length)==max(sapply(int.trts.list,length)))
      intersect.trts <- int.trts.list[[type.to.remove]]
      ind <- (1:5)[-type.to.remove]
      intersect.genes <- genes.by.cancer[[ind[1]]]
      for (j in 2:4) {
        intersect.genes <- intersect(intersect.genes,genes.by.cancer[[ind[j]]])
      }
      use.types <- names(split.cm.data)[-type.to.remove]
    }
    
    cols.to.keep <- c("Model","Treatment","TumorType","RespScaled","logSurvScaled",intersect.genes)
  }
  if (no.genes) {
    cols.to.keep <- c("Model","Treatment","TumorType","RespScaled","logSurvScaled")
    if (use5) {
      intersect.trts <- trts.by.cancer[[1]]
      for (i in 2:5) {
        intersect.trts <- intersect(intersect.trts,trts.by.cancer[[i]])
      }
      use.types <- names(split.cm.data)
    }
    if (use4) {
      int.trts.list <- sapply(names(trts.by.cancer),function(x) NULL)
      for (i in 1:5) {
        ind <- (1:5)[-i]
        int.trts.list[[i]] <- trts.by.cancer[[ind[1]]]
        for (j in 2:4) {
          int.trts.list[[i]] <- intersect(int.trts.list[[i]],trts.by.cancer[[ind[j]]])
        }
      }
      type.to.remove <- which(sapply(int.trts.list,length)==max(sapply(int.trts.list,length)))
      intersect.trts <- int.trts.list[[type.to.remove]]
      use.types <- names(split.cm.data)[-type.to.remove]
    }
  }
  dat <- NULL
  for (i in 1:length(use.types)) {
    int.rows <- split.cm.data[[use.types[i]]]$Treatment%in%intersect.trts
    int.cols <- colnames(split.cm.data[[use.types[i]]])%in%cols.to.keep
    dat <- rbind(dat,split.cm.data[[use.types[i]]][int.rows,int.cols])
  }
  
  return(dat)
}

