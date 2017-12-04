

## combine GSE63617.1(3m) and GSE63617.2 (6m), and combine annotations as GSE63617.early


library(dplyr)
library(HelperFunctions)
#---------------------------------------------------------------------------#
## combine GSE64398 (various genotype and timepoint) to early: GSE64398.early
# combine GSE64398 (various genotype and timepoint) GSE64398.late
#---------------------------------------------------------------------------#
rm(list=setdiff(ls(),'home_dir'))
disease ='AD'
    source('config/config_wrappers.R')
    md_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv')
    df_info <- read.delim(md_info, comment.char = '#')

load('/home/bzhuang/ND_project_combined/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/GSE64398.1/results/GSE64398.1_objects_after_batch_correction.Rdata')


## get early samples
phase = 'early'
samples = df_info[which(df_info$Dataset =='GSE64398.1' & df_info$Phase ==phase), ]
samples$phase_geno = paste0(samples$Timepoint, samples$Genotype)
tmp = array_design
tmp$phase_geno = paste0(tmp$Timepoint, tmp$Genotype)
index = union(which(tmp$phase_geno %in% samples$phase_geno), which(tmp$Genotype =='WT')) ## get mouse model and wt index
array_design = array_design[index, ]%>% droplevels()

array_dat = array_dat[, as.character(array_design$Sample)]

dataset = 'GSE64398'
(result_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/results/'))
(plot_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/plots/', dataset, '_'))
dir.create(result_pre, recursive = T)
(result_pre= paste0(result_pre, dataset, '.', phase, '_'))
save(array_design, annotation, array_dat, dataset, platform, plot_pre, result_pre, msg,
     file = paste0(result_pre, "objects_after_batch_correction.Rdata"))

sample_order <- as.character(array_design[with(array_design, order(Genotype)), 'Sample'])
print(setdiff(colnames(array_dat), sample_order))

rm(list=setdiff(ls(),'home_dir'))
for(disease in c('AD')){
    source('config/config_wrappers.R')
    md_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv')
    df_info <- read.delim(md_info, comment.char = '#')
}
load('/home/bzhuang/ND_project_combined/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/GSE64398.1/results/GSE64398.1_objects_after_batch_correction.Rdata')

phase = 'late'
samples = df_info[which(df_info$Dataset =='GSE64398.1' & df_info$Phase ==phase), ]
samples$phase_geno = paste0(samples$Timepoint, samples$Genotype)
tmp = array_design
tmp$phase_geno = paste0(tmp$Timepoint, tmp$Genotype)
index = union(which(tmp$phase_geno %in% samples$phase_geno), which(tmp$Genotype =='WT')) ## get mouse model and wt index
array_design = array_design[index, ]%>% droplevels()

array_dat = array_dat[, as.character(array_design$Sample)]

dataset = 'GSE64398'
(result_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/results/'))
(plot_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/plots/', dataset, '_'))
dir.create(result_pre, recursive = T, showWarnings = F)
(result_pre= paste0(result_pre, dataset, '.', phase, '_'))
save(array_design, annotation, array_dat, dataset, platform, plot_pre, result_pre, msg,
     file = paste0(result_pre, "objects_after_batch_correction.Rdata"))
sample_order <- as.character(array_design[with(array_design, order(Genotype)), 'Sample'])
print(setdiff(colnames(array_dat), sample_order))

#---------------------------------------------------------------------------# 
## combine GSE63617.1(3m) and GSE63617.2 (6m), as GSE63617.early, anotation is the same for both platforms
#---------------------------------------------------------------------------#
## all 3m samples
rm(list=setdiff(ls(),'home_dir'))
load('/home/bzhuang/ND_project_combined/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/GSE63617.1/results/GSE63617.1_objects.Rdata')
array_dat1 = array_dat
array_design1 = array_design

# get 6 months only
load('/home/bzhuang/ND_project_combined/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/GSE63617.2/results/GSE63617.2_objects.Rdata')
array_design = array_design[which(array_design$Timepoint == '6_months'), ]



## combine design
array_design = rbind(array_design1, array_design) %>% droplevels()

## combine dat
df = mapBindAllColumns(array_dat1, array_dat)

## get all samples and reorder
array_dat = df[, as.character(array_design$Sample)]
str(array_design)


dataset = 'GSE63617'
phase ='early'
(result_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/results/'))
(plot_pre= paste0(home_dir, "/AD_mouse_model_project/data_and_QC/all_data/explore_data/analysis_mode/",dataset, '.', phase,'/plots/', dataset, '_'))
dir.create(result_pre, recursive = T, showWarnings = F)
(result_pre= paste0(result_pre, dataset, '.', phase, '_'))
save(array_design, annotation, array_dat, dataset, platform, plot_pre, result_pre, msg,
     file = paste0(result_pre, "objects.Rdata"))

sample_order <- as.character(array_design[with(array_design, order(Genotype)), 'Sample'])
print(setdiff(colnames(array_dat), sample_order))



cat("
#---------------------------------------------------------------------------#
# PART 3.1 LIMMA DEA
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
# PART 3.1A LIMMA DEA for gemma and CEL
#---------------------------------------------------------------------------#\n")
## output put all limma (gemma and CEL) in 1 folder

### limma for Fisher only -CEL and gemma for 
rm(list=setdiff(ls(),'home_dir'))
 

disease_ls = c('AD', 'HD')
phase_ls = c('early', 'late')
q_threshold <- 0.1


for(disease in disease_ls){
    source('config/config_wrappers.R')
    for(phase in phase_ls){
        datadir <-paste0(limma_dir, 'limma_combined_genotype/', phase, '/')
        dir.create(datadir, showWarnings=F, recursive = T)
        object_dir <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
        f_model <- paste0(disease_dir, 'config_files/limma_genotype_combined_',phase, '.tsv')
        start_row <- 1
        combined_genotype =T
        
        #dataset_todo <-c('GSE48622')
        # to_do_ls <-1
        source('limma_DEA.R')
    }
}



###############
rm(list=setdiff(ls(),'home_dir'))
disease ='AD'
phase = 'early'
 

source('config/config_wrappers.R')
q_threshold <- 0.1
datadir <-paste0(limma_dir, 'limma_combined_genotype/', phase, '/')
dir.create(datadir, showWarnings=F, recursive = T)
object_dir <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
f_model <- paste0(disease_dir, 'config_files/limma_genotype_combined_',phase, '.tsv')

combined_genotype =T

dataset_todo <-c('GSE63617.early','GSE64398.early')
# to_do_ls <-1
source('limma_DEA.R')


#############
rm(list=setdiff(ls(),'home_dir'))
disease ='AD'
phase = 'late'
 
source('config/config_wrappers.R')
q_threshold <- 0.1
datadir <-paste0(limma_dir, 'limma_combined_genotype/', phase, '/')
dir.create(datadir, showWarnings=F, recursive = T)
object_dir <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
f_model <- paste0(disease_dir, 'config_files/limma_genotype_combined_',phase, '.tsv')
combined_genotype =T

dataset_todo <-c('GSE64398.late')
# to_do_ls <-1
source('limma_DEA.R')