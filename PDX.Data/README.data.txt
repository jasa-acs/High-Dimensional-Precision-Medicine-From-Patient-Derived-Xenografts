PDX DATA FILES:

RAW DATA SETS:

nm.3954-S2_copy_number_common_samples.txt: copy number data

nm.3954-S2_mutation_and_cn2_common_samples.txt: mutation data

nm.3954-S2_PCT_curve_metrics_common_samples.txt: response data

nm.3954-S2_PCT_raw_data_common_samples.txt: raw response data

nm.3954-S2_rnaseq_common_samples.txt: RNA_seq data

SCREENED GENE DATA SETS:

tumortype.numgenes.biomarkers.csv: genomic data for given tumor type
and number of genes; these are produced by the function collate.R

tumortype.numgenes.biomarkers.feat.new.csv: features from 
deep autoencoder for given tumor type and number of genes

PROCESSED DATA SETS:

split.cm.data.rda: R list of all data after unsupervised data 
(including both outcomes and gene data); this is produced by 
pdx.import.fn.R

trts.by.cancer.rda: R list of treatment names for each cancer; 
this is produced by pdx.import.fn.R

genes.by.cancer.rda: R list of genes remaining after 
unsupervised screening; this is produced by pdx.import.fn.R

pred_vals_rf.rda: R list of predicted values from random forest
by cancer type

pred_vals_rf_full.rda: R list of predicted values from 
random forest for combined cancer types

full.data.rda: clinical and genomic data for combined cancer types 
