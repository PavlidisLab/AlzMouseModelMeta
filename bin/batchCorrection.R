#' Created: 2015-12-10
#' Updated: 2016-01-06
#' Batch correction
#' ----------------------------------------------------- #
#' INPUT:
#'  - dataset (string): e.g. 'GSE36237'
#'  - object_dir(string): a directory contains the '*_objects.Rdata' from the explore_data_general.R script
#'      e.g. object_dir = '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/GSE36237/results/'
#'  
#' USAGE:
# source('batchCorrection.R')
# batchCorrection('GSE36237', object_dir = '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/GSE36237/results/')
#' ----------------------------------------------------- #
#' 


 
source("ComBat.R")
source("helper_functions.R")

## set defaults
## image size (pch)
width = 1600
height = 1000

batchCorrection <- function(dataset, object_dir){
  #' given a dataset name, the R objects from the data pre-process (explore_data_general.R) are loaded
  #' Outliers have been previously removed and unwanted sample are removed (e.g. select only the hippocampus 
  #' samples, not the cortex samples etc.)
  #' 
  #' INPUT:
  #'    -dataset: experiment id (string). e.g. GSE1556
  #'    -object_dir: dir for R objects, e.g. '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/GSE1556/results/'
  ####### sanity check ########
  if(class(dataset) != "character"){
    stop(paste0("Input dataset is a ", class(dataset), ". A string is required, e.g. 'GSE1234'."))
  }
  if(!dir.exists(object_dir)){
    stop(paste0("Error: Dir for R object file ", object_dir, " is not found."))
  }
  
  ###### input dataset must match with input object####
  #dataset<-'GSE36237'
  tmp <- grep('GSE', unlist(strsplit(object_dir, split = '\\/|\\_')), value =T)
  if(as.character(tmp) != as.character(dataset)){
    stop(paste0(dataset, "does not match with ", object_dir))
  }
  
  ###### get the dir for the R object #####
  f_object <- grep('_objects.Rdata', list.files(object_dir,full.names = T), value =T)
  
  if(length(f_object) == 0){  ##try to look for results subfolder
    object_dir <- paste0(object_dir,'/results/')
    f_object <- grep('_objects.Rdata', list.files(object_dir,full.names = T), value =T)
  }
  
  if(length(f_object) == 0){
    stop(paste0("Error: no R object file is found in ", object_dir))
  }
  
  ###### LOAD DATA
  print('###### loading R data ######')
  print(f_object)
  load(f_object)
  
  ###### PRE-PROCESS for combat input
  print(paste0('##### process the matrix and design file for ComBat input #####'))
  f_array <- paste0(result_pre, "array_dat_before_batch_correction.tsv")
  f_array_after_batch <- paste0(result_pre, "array_dat_after_batch_correction.tsv")
  f_design <- paste0(result_pre, "design_batch.tsv")
  write.table(array_dat, file = f_array, quote =F, sep ='\t', 
              row.names =F)
  ## including all factors
  response_ls <- intersect(x_col, colnames(array_design))
  response_ls <- c('Sample', 'Batch', setdiff(response_ls, c('Sample', 'Batch')))
  covariate_ls <- NULL
  for(i in response_ls){
      if(length(levels(array_design[,i]))>1){
          covariate_ls <- c(covariate_ls, i)
      }
  }
  
  ## for GSE20547.2 (PD) skip the gender
  if(dataset == 'GSE20547.2'){
      covariate_ls <- setdiff(covariate_ls, 'Gender')
  }
  #covariate_ls <- c('Sample', 'Batch')
  write.table(array_design[, covariate_ls], file = f_design, quote =F, sep ='\t', 
              row.names =F)
  
  ###### COMBAT
  print(paste0('##### ComBat #####'))
  
  df <- ComBat(f_array, f_design, skip=1, write=F, prior.plots=F) # don't output plot
  ## combat replaces first column name with 'geneinfo' and removed row names
  colnames(df) <- colnames(array_dat)
  rownames(df) <- rownames(array_dat)
  
  ## write after batch corrected files
  array_dat <- df
  write.table(array_dat, file = f_array_after_batch, quote =F, sep ='\t', 
              row.names =T)
  
  ##### PLOTS
  plot_pre <- paste0(plot_pre, 'combat_corrected_')
  print(paste0('##### plotting #####'))
  prefix <- " (batch_corrected)"
  plotFunctions(array_dat, plot_pre, dataset, array_design, prefix)
  
  
  ##### MODEL TEST
  ## test the factor and expression matrix
  response_ls <- intersect(x_col, colnames(array_design))
  response_ls <- response_ls[-length(response_ls)] # remove the last element, Sample
  msg_model <- modelTest(array_dat, array_design, response_ls)
  
  ##### LOG AND SAVE DATA
  print(paste0('##### log and save data #####'))
  ## log
  msg_batch <- paste0("\n#-------------------------------#\n## Batch correction\n#-------------------------------#\n", 
                     Sys.Date(),
                     '\nnumber of rows: ', nrow(array_dat),
                     '\nnumber of samples: ', ncol(array_dat))
  
  ## save process messages
  msg <-paste0(msg, msg_batch, msg_model)
  sink(paste0(result_pre, "batch_correction_log.txt"), type="output")
  writeLines(msg)
  sink()
  
  ## save r objects
  f_new_object <- paste0(result_pre, "objects_after_batch_correction.Rdata")
  save(array_design, annotation, array_dat, dataset, platform, plot_pre, result_pre, msg,
       file = f_new_object )
  print(paste0("R object saved: ", f_new_object))
}




