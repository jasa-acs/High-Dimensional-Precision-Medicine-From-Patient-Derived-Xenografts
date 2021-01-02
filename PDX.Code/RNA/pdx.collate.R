# pdx.collate.R

# collates results from supervised screening

library(stringr)
library(stringi)

# collate function
# takes as input a tumor type
# nums is the number of genes that should be kept -- can be a vector
# input_dir is where results of prognostic/redictive screening are (scratch directory)
# output_dir is where collated results should be sent (main directory)
pdx.collate = function(tumor_type, nums, input_dir, output_dir) {
  
  setwd(input_dir)
  
  # removes the suffix from file names
  sub = function(string)  return(stri_sub(string, 1, -5))

  csv.list = list.files(pattern = "*.csv")
  csv.list = csv.list[grep(tumor_type, csv.list)] 
  csv.list = csv.list[grep("screen", csv.list)]
  if (length(csv.list) == 0) print("No files with screening results in this directory!")
  
  # list of data frames
  df.list = list()

  # loop through files and import
  for (i in 1:length(csv.list)) {
  
    df.name = paste(sub(csv.list[i]), ".df")
    temp = read.csv(csv.list[i], header = T)
  
    temp$bdc.p.value = p.adjust(temp$bdc.p.value, method = "BY")
    temp$cca.p.value = p.adjust(temp$cca.p.value, method = "BY")
  
    temp = temp[order(temp$gene.name), ]
  
    df.list[[i]] = temp
  
    assign(df.name, temp)
  
    rm(temp)
  
  }

  gene.names = NULL
  for (i in 1:length(csv.list)) {
  
    gene.names = c(gene.names, as.character(df.list[[i]]$gene.name))
  
  }

  gene.names = unique(gene.names)
  gene.names = gene.names[order(gene.names)]

  # new data frame to store minimum p-values for each gene
  new.df = data.frame(gene = gene.names, pval = NA, cor = NA, cca.pval = NA, cca.cor = NA)

  # fill in p-values
  for (i in 1:dim(new.df)[1]) {
  
    cur.gene = as.character(new.df$gene[i])
  
    pvalues = NULL
    correlations = NULL
    cca.pvalues = NULL
    cca.cors = NULL
  
    for (j in 1:length(df.list)) {
    
      temp.df = df.list[[j]]
      ind = which(as.character(temp.df$gene.name) == cur.gene)
      pvalues = c(pvalues, temp.df$bdc.p.value[ind])
      correlations = c(correlations, temp.df$correlation[ind])
      cca.pvalues = c(cca.pvalues, temp.df$cca.p.value[ind])
      cca.cors = c(cca.cors, temp.df$cca.stat[ind])
    
    }
  
    new.df$pval[i] = min(pvalues, na.rm = T)
    new.df$cca.pval[i] = min(cca.pvalues, na.rm = T)
    temp.inds = which(pvalues == min(pvalues, na.rm = T))
    new.df$cor[i] = max(abs(correlations[temp.inds]))
    temp.inds2 = which(cca.pvalues == min(cca.pvalues, na.rm = T))
    new.df$cca.cor[i] = max(abs(cca.cors[temp.inds2]), na.rm = T)
  
  }

  setwd(output_dir)
  load('split.cm.data.rda')
  load('full.data.rda')
  
  for (k in 1:length(nums)) {
  
    final = new.df[order(new.df$pval), ]
    
    # number of genes with two zero p-values
    num_zero = length(final$gene[which(final$pval == 0  & final$cca.pval == 0)])
    
    num_genes = nums[k]
    
    # if there are more zero p-values than genes requested
    if (num_genes < num_zero) {
      final = final[which(final$pval == 0 & final$cca.pval == 0), ]
      final = final[order(-final$cca.cor), ]
      final = final[1:num_genes, ]
    }
  
    # if there are more genes requested than zero p-values
    if (num_genes >= num_zero) {
      perc = num_genes / dim(final)[1]
      final = final[which(final$pval <= quantile(final$pval, perc, na.rm = T)), ]
    }

    # create file name for results
    filename = paste(tumor_type, ".", num_genes, ".results.rna.csv", sep = "")
  
    # this file contains gene names and p-values and correlation statistics
    write.csv(final, filename)
  
    # get gene data
    if (tumor_type == "BRCA") biomarkers = split.cm.data$BRCA[ , 18:ncol(split.cm.data$BRCA)]
    if (tumor_type == "PDAC") biomarkers = split.cm.data$PDAC[ , 18:ncol(split.cm.data$PDAC)]
    if (tumor_type == "NSCLC") biomarkers = split.cm.data$NSCLC[ , 18:ncol(split.cm.data$NSCLC)]
    if (tumor_type == "CRC") biomarkers = split.cm.data$CRC[ , 18:ncol(split.cm.data$CRC)]
    if (tumor_type == "CM") biomarkers = split.cm.data$CM[ , 18:ncol(split.cm.data$CM)]
    if (tumor_type == "overall") biomarkers = dat[ , 7:ncol(dat)]
  
    # select names of gene variables to use
    vars.to.keep = NULL
    for (i in 1:length(final$gene)) {
      #vars.to.keep = c(vars.to.keep, paste(final$gene[i], ".cn", sep = ""), paste(final$gene[i], ".rna", sep = ""), paste(final$gene[i], ".mut", sep = ""))
      vars.to.keep = c(vars.to.keep, paste(final$gene[i], ".rna", sep = ""))
    }
    biomarkers = biomarkers[ , which(colnames(biomarkers) %in% vars.to.keep)]
  
    filename2 = paste(tumor_type, ".", num_genes, ".biomarkers.rna.csv", sep = "")
  
    # this file contains gene data for screened-in genes
    write.csv(biomarkers, filename2)
  
  }  # end loop through length of num

}

