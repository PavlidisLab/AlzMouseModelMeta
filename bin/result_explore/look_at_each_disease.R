rm(list=ls())
home_dir= 'C:/Users/User/Documents/lab_results/ND_project_combined/'
rm(list=setdiff(ls(),'home_dir'))
setwd(paste0("/bin/mouse_dataset_process/"))




df_re <- read.delim(paste0(home_dir, '/AD_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/early/mixed_model_results_with_expression.tsv'), comment.char = '#')
emfd <- paste0(home_dir, '/AD_mouse_model_project/ermineJ/mixed_model/random_intercept_include_NA_low_exp_rm/2016-08-16/analysis/early_down_regulation_mixed_model.tsv')
emfu <- paste0(home_dir, '/AD_mouse_model_project/ermineJ/mixed_model/random_intercept_include_NA_low_exp_rm/2016-08-16/analysis/early_down_regulation_mixed_model.tsv')


df_gene <- read.delim(paste0(home_dir, '/ND_results/MGI_info/annotation_precessed_omim.tsv'),comment.char = '#')

library(dplyr)
df <- left_join(df_re[, c(1:8)], df_gene[, c(1, 5,7, 13,14,6,8:10,2:4)]) ## mixed model results

em <- read.delim(emfu, comment.char = '#') 
em$regulation ='up'
emd <- read.delim(emfd, comment.char = '#') 
emd$regulation ='down'
em <- rbind(em, emd)
em <- em[, c(1,12,13, 2:11)]  ## erminj j result up and down top 50
rm(emd)


