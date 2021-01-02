dir1 = "~/PDX"  # main directory
dir2 = "/netscr/luckett"  # scratch directory
setwd(dir1)

source('~/PDX/pdx.collate.R')
pdx.collate(tumor_type = 'overall', nums = c(50, 100, 500, 1000), input_dir = dir2, output_dir = dir1)