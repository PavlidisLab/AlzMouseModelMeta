# Huntington's, download CEL files, and data QC 
# 2016-05-10
#---------------------------------------------------------------------------#
# PART 0A: DOWNLOHD CEL FILES
#---------------------------------------------------------------------------#
rm(list=setdiff(ls(),'home_dir'))
 
source('preprocess/GEO_download_CEL.R')
disease = 'HD'
disease_dir <- paste0(home_dir, disease, '_mouse_model_project/')
CEL_dir <- paste0(disease_dir, 'data_and_QC/GEO_data/CEL_raw/')

## total 13 datasets
GEO <- c('GSE10202', 'GSE62210', 'GSE48104', 'GSE9857', 'GSE26317', 
         'GSE25232','GSE9375','GSE9038', 'GSE50379', 
         'GSE32417', 'GSE64386', 'GSE19676', 'GSE7958')


## download raw files
downloadCEL(GEO, CEL_dir, meta_data = F)


## unzip all the .tar, .gz files
(unzip_files <- grep('GSE', list.dirs(CEL_dir, recursive = T, full.names = T), value = T))
(unzip_files <- grep(paste0(GEO, collapse = '|'), list.dirs(CEL_dir, recursive = T, full.names = T), value = T))

msg_all <- '#!/bin/bash'
for(unzip_f in unzip_files){
    (tar_f <-list.files(unzip_f,full.names = T, pattern = '.tar') )
    #print(tar_f)
    if(length(tar_f)>0){
        print(tar_f)
        
        (msg <- paste0('\ncd ', unzip_f,
                       '\ntar -xvf ', tar_f,
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

(cmd <- paste0('sh ', f_out))
system(cmd)
## rm tar files
(tar_files <-list.files(CEL_dir, recursive = T, full.names = T, pattern = 'tar'))
file.remove(tar_files)


#---------------------------------------------------------------------------#
# PART 0C: Quality control
#---------------------------------------------------------------------------#
## if errors like below appeared, restart R and re-run (eg. while wruning GSE9375)
# In qc.affy(unnormalised, ...) :
#   CDF Environment name ' mouse4302cdf ' does not match cdfname ' mgu74av2cdf '

#*****************************#
#### PART 0C.1.1 Affy platforms
#*****************************#

rm(list=setdiff(ls(),'home_dir'))

source('preprocess/data_QC_affy.R')
file_dir <- (paste0(home_dir, '/HD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/'))
design_f=paste0('../configs/HD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv') ## the experimental design file for all HD datasets
out_dir <- paste0(home_dir, '/HD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
data_dir <- paste0(home_dir, '/HD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')


GEO <- c(        "GSE10202",
                 "GSE48104",
                 "GSE7958",
                 "GSE9038",
                 "GSE9375",
                 "GSE9857")  ## affy, non GPL6246 platforms

GEO <- c( 
                 "GSE9375",
                 "GSE9857") 
(file_dir_ls <- grep(paste0(GEO, collapse = '|'), list.dirs(file_dir), value=T))

for(i in file_dir_ls){
    dataQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, result_dir =out_dir,chip_img = F, nuse_rle =T, mk_sample_info =T)
}



GEO <- c("GSE32417",
         "GSE50379",
         "GSE62210")  ## affy GPL6246, NUSE and RLE not available
(file_dir_ls <- grep(paste0(GEO, collapse="|"), list.dirs(file_dir), value=T))

for(i in file_dir_ls){
    dataQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, 
           result_dir =out_dir, chip_img = F, nuse_rle =F)
}


# 
# #*****************************#
# #### Illumnia platforms  -- just use gemma
# #*****************************#
# 
# 
rm(list=setdiff(ls(),'home_dir'))
 
source('preprocess/data_QC_illumina.R')
disease = 'HD'
source('config/config_wrappers.R')

file_dir <- paste0(disease_dir,'/data_and_QC/GEO_data/CEL_raw/')
(file_dir_ls <- grep("GSE", list.dirs(file_dir), value=T))

GEO_to_do <- c("GSE19676",
               "GSE25232",
               "GSE26317",
               "GSE64386")



(file_dir_ls <- grep(paste0(GEO_to_do, collapse="|"), list.dirs(file_dir), value=T))


design_f=paste0(disease_dir, '/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv')
out_dir <- paste0(disease_dir, '/data_and_QC/GEO_data/CEL_QC_3/')
data_dir <- paste0(disease_dir, '/data_and_QC/GEO_data/normalized_matrix/')


for(i in file_dir_ls){
    dataIlluminaQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, result_dir =out_dir)
}




#---------------------------------------------------------------------------#
# PART 1A: EXPLORE Gemma and CEL
#---------------------------------------------------------------------------#

# to do general data explore for one or multiple datasets with output folder defined, f_dataset_ls , 
# otherwise defaults are
# default output folder is output_folder <- '/home/bzhuang/HD_mouse_model_project/data_and_QC/gemma_data/explore_data/'
# default input list of datasets: f_dataset_ls <- "/home/bzhuang/HD_mouse_model_project/config_files/mouse_datasets_explore.tsv"
# 
# for all explore (remove outlier and subset for hippocampus only, no batch correction yet)

# **********************#
# **1A.a for exploring datasets -gemma
# **********************#
rm(list=setdiff(ls(),'home_dir'))
 
disease <- 'HD'
source('config/config_wrappers.R')
rm_outlier <- T
to_subset <- T
model_test <- T 
f_dataset_ls <- paste0(disease_dir, 'config_files/gemma_mouse_datasets_explore.tsv')


output_folder <- paste0(disease_dir, 'data_and_QC/gemma_data/explore_data/')
explore_method <- 'explore_mode'

unfiltered_exp_folder <- paste0(disease_dir, 'data_and_QC/gemma_data/unfiltered_exp_data/')  ## gemma processed

start_row <- 1
#to_do_ls <- 10
#dataset_todo <-'GSE1556'
source('explore_data_general.R')
exploreExpData()



#**********************#
#**1A.b for exploring datasets -CEL
#**********************#
rm(list=setdiff(ls(),'home_dir'))
 
disease <- 'HD'
source('config/config_wrappers.R')
rm_outlier <- T
to_subset <- T
f_dataset_ls <- paste0(disease_dir, 'config_files/CEL_mouse_datasets_explore.tsv')


output_folder <- paste0(disease_dir, 'data_and_QC/GEO_data/explore_data/')  ## CEL not processed
explore_method <- 'explore_mode'

unfiltered_exp_folder <- paste0(disease_dir, 'data_and_QC/GEO_data/normalized_matrix/')
model_test <- T  # if T do all test, if F only for subseted data
start_row <- 1
#to_do_ls <- 18
#dataset_todo <-c('GSE63617.1', 'GSE63617.2')
#dataset_todo <-c('GSE32417')
source('explore_data_general.R')
exploreExpData()
