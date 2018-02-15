#' 2015-12-01
#' updated 2017-04-12 # publish_codes
## for AD
#' part 0A: download the CEL files from GEO (gemma files need to be downloaded manually)
#' part 0B: prepare the design files
#' part 0C: quality check of CEL files
#' Part 1: explore data (CEl and gemma seperately, for explore to look at gender, batch effect etc)

#' all tar, gz files must be unzipped

#---------------------------------------------------------------------------#
# PART 0A: DOWNLOAD CEL and OTHER RAW FILES
# annotation files should be downloaded from https://gemma.msl.ubc.ca/home.html
#---------------------------------------------------------------------------#
print('PART 0A: DOWNLOAD CEL and OTHER RAW FILES')

setwd(file.path(here::here(),'bin'))
source('config_wrappers.R')
rm(list=setdiff(ls(),'home_dir'))
disease ='AD'
disease_dir <- paste0(home_dir, disease, '_mouse_model_project/')

source('preprocess/GEO_download_CEL.R')

GEO <- c( 'GSE36237', 'GSE52022', 'GSE14499', 'GSE53480', 'GSE48622') # affy
GEO <- c(GEO, 'GSE50521')  ## affy exon array
GEO <- c(GEO, 'GSE15056', 'GSE63617') ## agilent
GEO <- c(GEO,'GSE64398') ## illumina idat
# GEO <- c(GEO, 'GSE1556')  # no raw: file No supplemental files found.Download this from 


(CEL_dir <- paste0(disease_dir, 'data_and_QC/GEO_data/CEL_raw/'))
dir.create(CEL_dir, recursive = T, showWarnings = F)

## download raw data
downloadCEL(GEO, CEL_dir, meta_data = F)

#####
## unzip files # or manually unzip the files
####
(unzip_files <- grep('GSE', list.dirs(CEL_dir, recursive = T, full.names = T), value = T))
(unzip_files <- grep(paste0(GEO, collapse = '|'), list.dirs(CEL_dir, recursive = T, full.names = T), value = T))

msg_all <- '#!/bin/bash'
for(unzip_f in unzip_files){
    (tar_f <-list.files(unzip_f,full.names = T, pattern = '.tar') )
    #print(tar_f)
    if(length(tar_f)>0){
        print(tar_f)
        
        (msg <- paste0('\ncd ', shQuote(unzip_f),'/',
                       '\ntar -xvf ', shQuote(tar_f),
                       '\ngunzip *.gz'
        ))
    }else{msg =''}
    
    msg_all <- paste0(msg_all, msg)
}
msg_all
(f_out <- paste0(CEL_dir,'/unzip_files.sh'))
sink(file = f_out, type = 'output')
writeLines(msg_all)
sink()

(cmd <- paste0('sh ', shQuote(f_out)))
system(cmd)

(tar_files <-list.files(CEL_dir, recursive = T, full.names = T, pattern = 'tar'))
# file.remove(tar_files)


#---------------------------------------------------------------------------#
# PART 0B: PREPARE GSE63617
#---------------------------------------------------------------------------#
print('PART 0B: PREPARE GSE63617')

# GSE63617 -> split to GSE63617.1 and GSE63617.2 due to samples using different platforms
source("helper_functions.R") ## for saveSampleInfo()

design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
df <- read.delim(design_f, comment.char = '#')
old_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE63617/')


for (gse in c('GSE63617.1', 'GSE63617.2')){
    (file_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/', gse, '/'))
   
    dir.create(file_dir, showWarnings = F)
    ## get the GSM ids
    gsm_id <- as.character(df[which(df$Dataset ==gse), 'ExternalID']%>%droplevels())
    gsm_ls <- grep(paste0(gsm_id, collapse = '|'), list.files(old_dir), value = T)

    ## move files to the new file dir
    file.rename(from = paste0(old_dir, gsm_ls), to = paste0(file_dir, gsm_ls))
    
    ## create corresponding design file
    saveSampleInfo(file_dir, design_f)
}



#---------------------------------------------------------------------------#
# PART 0B: PREPARE GSE64398
#---------------------------------------------------------------------------#
print('PART 0B: PREPARE GSE64398')
# GSE64398 -> split to GSE64398.1 and GSE64398.2, GSE64398.3 due to samples using different platforms
# GSE64398.1 is the Hippocampal samples

source("helper_functions.R") ## for saveSampleInfo()

design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
df <- read.delim(design_f, comment.char = '#')
old_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE64398/')


for (gse in c('GSE64398.1', 'GSE64398.2', 'GSE64398.3')){
    (file_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/', gse, '/'))
    
    dir.create(file_dir, showWarnings = F)
    ## get the GSM ids
    gsm_id <- as.character(df[which(df$Dataset ==gse), 'ExternalID']%>%droplevels())
    gsm_ls <- grep(paste0(gsm_id, collapse = '|'), list.files(old_dir), value = T)
    
    ## move files to the new file dir
    file.rename(from = paste0(old_dir, gsm_ls), to = paste0(file_dir, gsm_ls))
    
    ## create corresponding design file
    saveSampleInfo(file_dir, design_f)
}




#---------------------------------------------------------------------------#
# PART 0C: Quality check and normalized data
#---------------------------------------------------------------------------#
print('PART 0C: Quality check and normalized data')

#*****************************#
#### PART 0C.1.1 Affy platforms
#*****************************#
print('PART 0C.1.1 Affy platforms')

rm(list=setdiff(ls(),'home_dir'))

source('preprocess/data_QC_affy.R')
file_dir <- (paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/'))
design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')


GEO <- c( 'GSE36237', 'GSE52022', 'GSE14499', 'GSE53480', 'GSE48622') # affy
(file_dir_ls <- grep(paste0(GEO, collapse = '|'), list.dirs(file_dir), value=T))

for(i in file_dir_ls){
    dataQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, 
           result_dir =out_dir,chip_img = F, nuse_rle =T, mk_sample_info =T)
}


# hope this isn't important. Do not ignore warnings, they can become errors.
# GEO <- c('GSE50521')  ## affy exon array, NUSE and RLE not available, 
# this will output a warning, but a normalized expression table will be done
# (file_dir_ls <- grep(paste0(GEO, collapse="|"), list.dirs(file_dir), value=T))
# 
# for(i in file_dir_ls){
#     dataQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, 
#            result_dir =out_dir, chip_img = F, nuse_rle =F)
# }

#*****************************#
#### PART 0C.1.2 Affy exon platforms
#*****************************#
print('PART 0C.1.2 Affy exon platforms')

#### Exon affy arrays: GPL6096
rm(list=setdiff(ls(),'home_dir'))

source('preprocess/data_QC_exon.R')
file_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/')
design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')


GEO <- c('GSE50521')  ## affy exon array, NUSE and RLE not available
(file_dir_ls <- grep(paste0(GEO, collapse="|"), list.dirs(file_dir), value=T))

for(i in file_dir_ls){
    dataExonArrayQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, result_dir =out_dir,chip_img = F, rle =T, mk_sample_info =T)
}





#*****************************#
#### PART 0C.2.1 Aligent platforms, green only ## need to retart session
#*****************************#
########
#### # preprocess GSE15056, 1 color, only Cy3; no dye-swaps 
########
print('PART 0C.2.1 Aligent platforms, green only ## need to retart session')

source("preprocess/data_QC_agilent.R")

GEO <- c( 'GSE15056')
file_dir <- (paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/', GEO, '/'))
design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')


saveSampleInfo(file_dir, design_f)
f_target <-paste0(file_dir, "SampleInfo.txt")
aligentDataQC(file_dir, f_target = f_target,agilent_source="agilent.median", channel ='green',
              out_dir = out_dir, width = 1600, height=1000, data_dir =data_dir, return_value =F)



########
#### # preprocess GSE63617
########
rm(list=setdiff(ls(),'home_dir'))

source("preprocess/data_QC_agilent.R")


file_dir <- (paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/'))
design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')


GEO <- c( 'GSE63617.1', 'GSE63617.2')
(file_dir_ls <- grep(paste0(GEO, collapse = '|'), list.dirs(file_dir), value=T))

for(i in file_dir_ls){
    (f_target <-paste0(i, '/SampleInfo.txt'))
    aligentDataQC(file_dir = i, f_target = f_target,agilent_source="agilent.median", channel ='green',
                  out_dir = out_dir, width = 1600, height=1000, data_dir =data_dir, return_value =F)
    }



#*****************************#
#### PART 0C.3.1 Illumnia platforms  --iDat data
#*****************************#
print('PART 0C.3.1 Illumnia platforms  --iDat data')

rm(list=setdiff(ls(),'home_dir'))

file_dir <- (paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/'))
design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all AD datasets
out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/other_QC/')
data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')

source('preprocess/data_QC_illumina_iDat.R')

bgx <- '../configs/AD_mouse_dataset_doc/bgx/MouseRef-8_V2_0_R2_11278551_A.bgx'

GEO_to_do<- c('GSE64398.1')
(file_dir_ls <- grep(paste0(GEO_to_do, collapse="|"), list.dirs(file_dir), value=T))

for(file_dir in file_dir_ls){
    iDatDataQC(file_dir, out_dir = out_dir, data_dir =data_dir, bgx = bgx,
                   design_f=design_f)
}



#---------------------------------------------------------------------------#
# PART 1A: EXPLORE Gemma and CEL
#---------------------------------------------------------------------------#
# to do general data explore for one or multiple datasets with output folder defined, f_dataset_ls , 
# otherwise defaults are
#' default output folder is output_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/'
#' default input list of datasets: f_dataset_ls <- "/home/bzhuang/AD_mouse_model_project/config_files/mouse_datasets_explore.tsv"
#' 
#' for all explore (remove outlier and subset for hippocampus only, no batch correction yet)
#' 
#**********************#
#**1A.a for exploring datasets -gemma
#**********************#
# rm(list=setdiff(ls(),'home_dir'))
# 
# 
# source('config/config_wrappers.R')
# source('explore_data_general.R')
# rm_outlier <- T
# to_subset <- T
# f_dataset_ls <- paste0(disease_dir, 'config_files/gemma_mouse_datasets_explore.tsv')
# disease <- 'AD'
# 
# output_folder <- paste0(disease_dir, 'data_and_QC/gemma_data/explore_data/')
# explore_method <- 'explore_mode'
# 
# unfiltered_exp_folder <- paste0(disease_dir, 'data_and_QC/gemma_data/unfiltered_exp_data/')  ## gemma processed
# 
# start_row <- 1
# #to_do_ls <- 4
# #dataset_todo <-'GSE1556'
# source('explore_data_general.R')
# exploreExpData()



#**********************#
#**1A.b for exploring datasets -CEL
#**********************#
# rm(list=setdiff(ls(),'home_dir'))
# disease='AD'
# 
# source('config/config_wrappers.R')
# 
# 
# rm_outlier <- T
# to_subset <- T
# f_dataset_ls <- paste0(disease_dir, 'config_files/CEL_mouse_datasets_explore.tsv')
# output_folder <- paste0(disease_dir, 'data_and_QC/GEO_data/explore_data/')  ## CEL not processed
# explore_method <- 'explore_mode'
# 
# unfiltered_exp_folder <- paste0(disease_dir, 'data_and_QC/GEO_data/normalized_matrix/')
# model_test <- T  # if T do all test, if F only for subseted data
# #start_row <- 25
# #to_do_ls <- 25
# # dataset_todo <-c('GSE63617.1', 'GSE63617.2')
# #dataset_todo <-c('GSE64398.1','GSE64398.2','GSE64398.3')
# dataset_todo <-c('GSE15056')
# dataset_todo <-c('GSE53480')
# 
# source('explore_data_general.R')
# exploreExpData()
