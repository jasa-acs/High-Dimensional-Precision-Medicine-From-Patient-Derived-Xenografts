# pdx.import.R

# imports pdx data

# Function imports PDX datasets and performs unsupervised screening of treatments and genetic features
# perc1: remove RNA-seq features with MAD in the lowest perc1 percentile
# cut2: remove RNA-seq features with mean expression < cut2
# perc2: remove copy number features with MAD in the lowest perc2 percentile
# cut4, cut5: remove mutation features when < cut4 or > cut5 of lines have the mutation
# output.interim: if true, outputs interim (post-unsupervised screening) datasets for RNA, CN, and MUT
import.pdx <- function(input_dir = "./", output_dir = "./", perc1 = 0.2, cut2 = 1, perc2 = 0.2, cut4 = 0.1, cut5 = 0.9, output.interim = F) {
  library(stringr)
  library(stringi)
  
  # copy number data
  cn.file = "nm.3954-S2_copy_number_common_samples.txt"
  
  # mutation data
  mut.file = "nm.3954-S2_mutation_and_cn2_common_samples.txt"
  
  # curve metrics data
  cm.file = "nm.3954-S2_PCT_curve_metrics_common_samples.txt"
  
  # raw data
  rd.file = "nm.3954-S2_PCT_raw_data_common_samples.txt"
  
  # rna seq data
  rna.file = "nm.3954-S2_rnaseq_common_samples.txt"
  
  cn.file = paste(dir1, cn.file, sep = "")
  mut.file = paste(dir1, mut.file, sep = "")
  cm.file = paste(dir1, cm.file, sep = "")
  rd.file = paste(dir1, rd.file, sep = "")
  rna.file = paste(dir1, rna.file, sep = "")
  
  # read in data
  cn.data = read.table(file = cn.file, header = T, sep = "\t",strip.white=T)
  mut.data = read.table(file = mut.file, header = T, sep = "\t",strip.white=T)
  cm.data = read.table(file = cm.file, header = T, sep = "\t", quote = "", strip.white=T)
  rd.data = read.table(file = rd.file, header = T, sep = "\t",strip.white=T, quote = "")
  rna.data = read.table(file = rna.file, header = T, sep = "\t",strip.white=T)
  
  # create (character) list of unique gene names
  genelist.cn <- as.character(cn.data$X[-c(1,2)])
  genelist.rna <- as.character(rna.data$X)
  genelist.mut <- as.character(mut.data$Gene)
  genelist <- unique(c(genelist.cn,genelist.rna,genelist.mut))
  
  # change row names so that copy number and RNA-seq for the same gene will have different names
  cn.data$X = noquote(paste(cn.data$X, ".cn", sep = ""))
  rna.data$X = noquote(paste(rna.data$X, ".rna", sep = ""))
  
  # transpose so that each row is a sample and each column is a gene copy number
  cn.names = cn.data$X
  cn.data = as.data.frame(t(cn.data[ , -1]))
  colnames(cn.data) = cn.names
  #cn.data$factor = factor(row.names(cn.data))
  
  # transpose so that each row is a sample and each column is a gene RNA-seq measurement
  rna.names = rna.data$X
  rna.data = as.data.frame(t(rna.data[ , -1]))
  colnames(rna.data) = rna.names
  #rna.data$factor = factor(row.names(rna.data))
  
  # initial unsupervised screening of RNA-seq and copy number
  
  # screen RNA-seq
  
  # remove RNA-seq data for genes when MAD is less than perc1 percentile
  rna.mad = apply(rna.data, 2, mad)
  cut1 = quantile(rna.mad, perc1)
  rna.data = rna.data[ , which(rna.mad > cut1)]
  
  # remove RNA-seq data for genes with less than cut2 mean expression
  rna.means = apply(rna.data, 2, mean)
  rna.data = rna.data[ , which(rna.means > cut2)]
  
  # screen copy number
  
  # remove copy number data for genes when MAD is less than perc2 percentile 
  cn.mad = apply(cn.data, 2, mad)
  cut3 = quantile(cn.mad, perc2)
  cn.data = cn.data[ , which(cn.mad > cut3)]
  
  # RNA-seq and copy number data together
  gene.data = matrix(0,nrow=nrow(cn.data),ncol=ncol(cn.data)+ncol(rna.data)+length(genelist))
  colnames(gene.data) <- c(colnames(cn.data),colnames(rna.data),noquote(paste(genelist, ".mut", sep = "")))
  rownames(gene.data) <- str_replace_all(rownames(cn.data),"[[:punct:]]","")
  gene.data[,1:ncol(cn.data)] <- as.matrix(cn.data)
  gene.data[,(ncol(cn.data)+1):(ncol(cn.data)+ncol(rna.data))] <- as.matrix(rna.data)
  
  # reformat mutation data lineID to match RNAseq and copy number
  mut.data$Sample <- str_replace_all(mut.data$Sample,"[[:punct:]]","")
  
  # merge on mutation data
  # if a gene is present in the mutation data for a PDX line, then that PDX line has a mutation for that gene;
  #    if the PDX line has no mutation, that gene does not appear in the mutation data at all
  # function takes a row of the mutation data and outputs the corresponding row and column of the merged genetic
  #    dataset (i.e. the specific cell that needs to be set equal to 1 to indicate mutation present)
  findcells <- function(mut.row) {
    row.of.merge <- which(rownames(gene.data)==mut.row["Sample"])
    col.of.merge <- paste(mut.row["Gene"],".mut",sep="")
    return(c(row.of.merge,col.of.merge))
  }
  mut.cells <- apply(mut.data,1,findcells)
  mut.cells[1,] <- as.numeric(mut.cells[1,])
  for (i in 1:ncol(mut.cells)) {
    gene.data[as.numeric(mut.cells[1,i]),mut.cells[2,i]] <- 1
  }
  
  # initial unsupervised screening of mutation data
  # remove mutation data on genes when less than cut4 or more than cut5 have mutations
  mut.temp = gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  mut.means = apply(mut.temp, 2, mean)
  mut.new = mut.temp[ , which(mut.means > cut4 & mut.means < cut5)]
  gene.data = cbind(gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], mut.new)
  
  # create new outcome variables in cm.data
  cm.data$NonPD <- as.numeric(cm.data$ResponseCategory!="PD")
  cm.data$CRorPR <- as.numeric(grepl("CR",cm.data$ResponseCategory,fixed=TRUE)|grepl("PR",
                                                                                     cm.data$ResponseCategory,fixed=T))
  # other outcome variables are BestAvgResponse and TimeToDouble
  
  # add cancer type to clinical data
  rd.data$Model <- str_replace_all(rd.data$Model,"[[:punct:]]","")
  cm.data$Model <- str_replace_all(cm.data$Model,"[[:punct:]]","")
  # function takes PDX line identifier and outputs the tumor type corresponding to that PDX line; this allows us
  #    to sapply over cm.data$Model below
  # note: some rows in rd.data have Tumor.type missing, so for a given PDX line, we must find the first row in rd.data
  #    that actually states the tumor type
  find.cancer.type <- function(model) {
    rowind <- which(rd.data$Model==model[1] & rd.data$Tumor.Type!="")[5]
    return(rd.data[rowind,'Tumor.Type'])
  }
  cm.data$TumorType <- sapply(cm.data$Model,find.cancer.type)
  # remove blank level of TumorType that is copied over thanks to blank rows in rd.data
  cm.data$TumorType <- as.character(cm.data$TumorType)
  cm.data$TumorType <- factor(cm.data$TumorType)
  
  # create data frames of gene data for each cancer type (for screening)
  temp.cm = cm.data[!duplicated(cm.data$Model), ]
  tumors = as.character(temp.cm$TumorType)
  brca.gene.data = gene.data[which(tumors == "BRCA"), ]
  pdac.gene.data = gene.data[which(tumors == "PDAC"), ]
  crc.gene.data = gene.data[which(tumors == "CRC"), ]
  nsclc.gene.data = gene.data[which(tumors == "NSCLC"), ]
  cm.gene.data = gene.data[which(tumors == "CM"), ]
  
  # screen mutations for BRCA
  brca.mut.temp = brca.gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  brca.mut.means = apply(brca.mut.temp, 2, mean)
  brca.mut.new = brca.mut.temp[ , which(brca.mut.means > cut4 & brca.mut.means < cut5)]
  brca.gene.data = cbind(brca.gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], brca.mut.new)
  
  dim(brca.mut.temp)
  dim(brca.mut.new)
  
  # screen mutations for PDAC
  pdac.mut.temp = pdac.gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  pdac.mut.means = apply(pdac.mut.temp, 2, mean)
  pdac.mut.new = pdac.mut.temp[ , which(pdac.mut.means > cut4 & pdac.mut.means < cut5)]
  pdac.gene.data = cbind(pdac.gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], pdac.mut.new)
  
  dim(pdac.mut.temp)
  dim(pdac.mut.new)
  
  # screen mutations for CRC
  crc.mut.temp = crc.gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  crc.mut.means = apply(crc.mut.temp, 2, mean)
  crc.mut.new = crc.mut.temp[ , which(crc.mut.means > cut4 & crc.mut.means < cut5)]
  crc.gene.data = cbind(crc.gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], crc.mut.new)
  
  dim(crc.mut.temp)
  dim(crc.mut.new)
  
  # screen mutations for NSCLC
  nsclc.mut.temp = nsclc.gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  nsclc.mut.means = apply(nsclc.mut.temp, 2, mean)
  nsclc.mut.new = nsclc.mut.temp[ , which(nsclc.mut.means > cut4 & nsclc.mut.means < cut5)]
  nsclc.gene.data = cbind(nsclc.gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], nsclc.mut.new)
  
  dim(nsclc.mut.temp)
  dim(nsclc.mut.new)
  
  # screen mutations for CM
  cm.mut.temp = cm.gene.data[ , (ncol(cn.data) + ncol(rna.data) + 1):ncol(gene.data)]
  cm.mut.means = apply(cm.mut.temp, 2, mean)
  cm.mut.new = cm.mut.temp[ , which(cm.mut.means > cut4 & cm.mut.means < cut5)]
  cm.gene.data = cbind(cm.gene.data[ , 1:(ncol(cn.data) + ncol(rna.data))], cm.mut.new)
  
  dim(cm.mut.temp)
  dim(cm.mut.new)
  
  # scale best average response
  cm.data$RespScaled = cm.data$BestAvgResponse / sd(cm.data$BestAvgResponse)
  
  # take log of survival and scale
  cm.data$logSurv = log(cm.data$TimeToDouble)
  cm.data$logSurvScaled = cm.data$logSurv / sd(cm.data$logSurv)
  
  # create list of post-screening gene datasets to work with loop below
  genedata.list <- list(brca.gene.data,cm.gene.data,crc.gene.data,nsclc.gene.data,pdac.gene.data)
  
  # create split clinical + post-unsupervised screening genetic data
  cancer.types <- unique(cm.data$TumorType)
  split.cm.data <- split(cm.data,cm.data$TumorType)
  untrt.level <- which(levels(cm.data$Treatment)=="untreated")
  # copy.genes function takes the PDX line identifier and returns the relevant row from the screened gene data
  #   to a new dataset; this allows us to sapply over cm.data$Model below
  copy.genes <- function(mouseid,genedata) {
    return(genedata[rownames(genedata)==mouseid,])
  }
  for (i in 1:length(split.cm.data)) {
    # set untreated as reference level of treatment
    split.cm.data[[i]]$Treatment <- relevel(split.cm.data[[i]]$Treatment,ref=untrt.level)
    # re-initialize treatment factor so that levels correspond to only treatments used for that tumor type
    split.cm.data[[i]]$Treatment <- factor(split.cm.data[[i]]$Treatment)
    split.cm.data[[i]] <- cbind(split.cm.data[[i]],t(sapply(split.cm.data[[i]][,"Model"],
                                                            function(x) copy.genes(x,genedata.list[[i]]))))
  }
  
  # remove anomalous PDX line that was only given one treatment
  split.cm.data$NSCLC <- split.cm.data$NSCLC[split.cm.data$NSCLC$Model!="X1795",]
  # remove anomalous PDX line that was not given untreated
  split.cm.data$PDAC <- split.cm.data$PDAC[split.cm.data$PDAC$Model!="X1959",]
  # remove treatments with insufficient sample size (present for <90% of PDX lines within tumor type)
  trt.remove.list.brca <- "LFW527 + everolimus"
  trt.remove.list.cm <- c("LEE011 + binimetinib","CLR457")
  trt.remove.list.crc <- c("WNT974", "abraxane+gemcitabine", "INC424 + binimetinib", "INC424", "gemcitabine-50mpk", 
                           "figitumumab\" + binimetinib", "figitumumab\"", "BKM120 + LDE225", "BKM120 + binimetinib",
                           "abraxane + gemcitabine", "abraxane", "trametinib", "binimetinib-3.5mpk",
                           "LFW527 + binimetinib")
  split.cm.data$BRCA <- split.cm.data$BRCA[!(split.cm.data$BRCA$Treatment%in%trt.remove.list.brca),]
  split.cm.data$CM <- split.cm.data$CM[!(split.cm.data$CM$Treatment%in%trt.remove.list.cm),]
  split.cm.data$CRC <- split.cm.data$CRC[!(split.cm.data$CRC$Treatment%in%trt.remove.list.crc),]
  # re-initialize treatment factor for tumor types that just had levels removed
  split.cm.data$BRCA$Treatment <- as.character(split.cm.data$BRCA$Treatment)
  split.cm.data$BRCA$Treatment <- factor(split.cm.data$BRCA$Treatment)
  split.cm.data$CM$Treatment <- as.character(split.cm.data$CM$Treatment)
  split.cm.data$CM$Treatment <- factor(split.cm.data$CM$Treatment)
  split.cm.data$CRC$Treatment <- as.character(split.cm.data$CRC$Treatment)
  split.cm.data$CRC$Treatment <- factor(split.cm.data$CRC$Treatment)
  
  # FINAL PRODUCTS: 
  # 1) a list of datasets, split by tumor type, containing all relevant clinical data and post-unsupervised
  #      screening genetic data
  # 2) list (i.e. names) of all genetic features left after unsupervised screening
  # 3) list (i.e. names) of all treatments left after unsupervised screening
  # save final product to directory
  save(split.cm.data,file = paste(dir2, "split.cm.data.rda", sep = ""))
  trts.by.cancer <- vector('list',5)
  genes.by.cancer <- vector('list',5)
  names(trts.by.cancer) <- names(split.cm.data) -> names(genes.by.cancer)
  for (i in 1:length(trts.by.cancer)) {
    trts.by.cancer[[i]] <- as.character(levels(split.cm.data[[i]]$Treatment))
    genestart <- which(grepl(".cn",colnames(split.cm.data[[i]])))[1]
    genes.by.cancer[[i]] <- colnames(split.cm.data[[i]])[genestart:ncol(split.cm.data[[i]])]
  }
  save(trts.by.cancer,file = paste(dir2, "trts.by.cancer.rda", sep = ""))
  save(genes.by.cancer,file = paste(dir2, "genes.by.cancer.rda", sep = ""))
  
  # output unsupervised screened genetic datasets
  if (output.interim) {
    cn.columns <- c(which(colnames(split.cm.data[[1]])%in%c("Model","Treatment")),
                    which(grepl(".cn",colnames(split.cm.data[[1]]))))
    rna.columns <- c(which(colnames(split.cm.data[[1]])%in%c("Model","Treatment")),
                     which(grepl(".rna",colnames(split.cm.data[[1]]))))
    mut.columns <- c(which(colnames(split.cm.data[[1]])%in%c("Model","Treatment")),
                     which(grepl(".mut",colnames(split.cm.data[[1]]))))
    write.csv(split.cm.data[[1]][,cn.columns],paste(dir2, "unsup_",names(split.cm.data)[1],"_cn.csv",sep=""))
    write.csv(split.cm.data[[1]][,rna.columns],paste(dir2, "unsup_",names(split.cm.data)[1],"_rna.csv",sep=""))
    write.csv(split.cm.data[[1]][,cn.columns],paste(dir2, "unsup_",names(split.cm.data)[1],"_mut.csv",sep=""))
  }

  return(NULL)
}

