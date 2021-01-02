# Commented out because this segment of code appears to do nothing
# unsupervised screening for prognostic biomarkers
#dir_in = "C:/Users/jgxchen/Dropbox/documents/UNC/forClass/1Projects/Michael112016_JASAPDX/data/owl"
#dir_code = "C:/Users/jgxchen/Dropbox/documents/UNC/forClass/1Projects/Michael112016_JASAPDX/code/deepL"
#setwd(dir_in)

### CHANGE DIRECTORIES HERE ###
dir_in = "/nas02/home/j/g/jgxchen/PDX_line/data/Screening_results"
dir_code = "/nas02/home/j/g/jgxchen/PDX_line/code/deepL"
setwd(dir_in)
cancer_type = c("BRCA", "CM", "CRC", "NSCLC", "PDAC")
genenum = c(50, 100, 500, 1000)

for (i in 1:length(cancer_type)) {
  for (j in 1:length(genenum)) {
    py = rep(0, 8)
    py[1] = "import sys" ## input
    py[2] = "import os"
    pathtodata = sprintf("'%s'", dir_in)
    pathtocode = sprintf("'%s'", dir_code)
    input = sprintf("%s.%i.biomarkers.rna", cancer_type[i], genenum[j])
    py[3] = sprintf("pathname=%s", pathtodata)
    py[4] = sprintf("sys.path.append(%s)", pathtocode)
    py[5] = sprintf("os.chdir(pathname)")
    py[6] = sprintf("from functions_bio import *")
    py[7] = sprintf("infile='%s'", input)
    py[8] = sprintf("autoencoder(infile)")
    #py[7] = sprintf("os.system(example_bio.py %s)", input)
    CMD = paste(py, collapse = "\n")
    pathtopyfile = sprintf("%s/DL_%s_%i.py", dir_code, cancer_type[i], genenum[j])
    write.table(CMD, file = pathtopyfile, quote = F, col.names = F, row.names = F)
    system(sprintf("bsub -o output.txt python CMD BATCH %s", pathtopyfile))
  }
}

#### longleaf submission
### CHANGE DIRECTORIES HERE ###
dir_in = "/nas/longleaf/home/jgxchen/PDX_line/data/Screening_results"
dir_code = "/nas/longleaf/home/jgxchen/PDX_line/code/deepL"
setwd(dir_code)
py_file = NULL
for (i in 1:length(cancer_type)) {
  for (j in 1:length(genenum)) {
    py = rep(0, 8)
    py[1] = "import sys" ## input
    py[2] = "import os"
    pathtodata = sprintf("'%s'", dir_in)
    pathtocode = sprintf("'%s'", dir_code)
    input = sprintf("%s.%i.biomarkers", cancer_type[i], genenum[j])
    py[3] = sprintf("pathname=%s", pathtodata)
    py[4] = sprintf("sys.path.append(%s)", pathtocode)
    py[5] = sprintf("os.chdir(pathname)")
    py[6] = sprintf("from functions_bio import *")
    py[7] = sprintf("infile='%s'", input)
    py[8] = sprintf("autoencoder(infile)")
    #py[7] = sprintf("os.system(example_bio.py %s)", input)
    CMD = paste(py, collapse = "\n")
    #pathtopyfile = sprintf("%s/DL_%s_%i.py", dir_code, cancer_type[i], genenum[j])
    pathtopyfile = sprintf("DL_%s_%i.py", cancer_type[i], genenum[j])
    py_file = c(py_file, sprintf("python %s", pathtopyfile))
    write.table(CMD, file = pathtopyfile, quote = F, col.names = F, row.names = F)
    #system(sprintf("bsub -o output.txt python CMD BATCH %s", pathtopyfile))
  }
}

write.table(paste(py_file, collapse = "\n"), file="subdeepL.sl", 
            quote = F, col.names = F, row.names = F)

for(i in 1:length(py_file)){
  py = rep(0, 9)
  py[1] = "#!/bin/bash"
  py[2] = ""
  py[3] = "#SBATCH -n 1"
  py[4] = "#SBATCH -p gpu"
  py[5] = "#SBATCH -t 10-00:00:00"
  py[6] = "#SBATCH --qos gpu_access"
  py[7] = "#SBATCH --gres=gpu:1"
  py[8] = ""
  py[9] = sprintf(py_file[i])
  CMD = paste(py, collapse = "\n")
  pathtopyfile = sprintf("subdeep%i.sl", i)
  write.table(CMD, file = pathtopyfile, quote = F, col.names = F, row.names = F)
}



######### rename the feature outputs
### CHANGE DIRECTORY HERE ###
feat_dir <- "C:/Users/jgxchen/Dropbox/documents/UNC/forClass/1Projects/Michael112016_JASAPDX/data/DL_results"
setwd(feat_dir)
fileList <- list.files(pattern=".csv")
for(i in fileList){
  dattmp = read.csv(i,  skip = 0)
  datout = dattmp[,-1]
  colnames(datout) <- paste("feat", 1:ncol(datout), sep="")
  tempname = gsub(".csv", "",i)
  write.csv(datout, file = sprintf("./rename/%s.new.csv", tempname), row.names=F)
}
