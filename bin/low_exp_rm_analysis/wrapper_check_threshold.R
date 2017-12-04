
# do not use this, use wrapper
#' 2016-03-15
#' updated 2016-04-11, add GSE52022
## wrapper for rm_low_exp_by_gender_genens.R

#' should remove GSE1556, this dataset has wierd gene-gene correlation
 
source("./low_exp_rm_analysis/rm_low_exp_by_gender_genes.R")


#############################
## to investigate different methods and threshold
#############################
method_t_ls <-c("max", "median", "75quantile", "95quantile") 

r_obj_ls <- c("/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE48622/results/GSE48622_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE52022/results/GSE52022_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.1/results/GSE63617.1_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.2/results/GSE63617.2_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.6m/results/GSE63617.6m_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.15m/results/GSE63617.15m_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE36237/results/GSE36237_objects_after_batch_correction.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE14499/results/GSE14499_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/analysis_mode/GSE1556/results/GSE1556_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/analysis_mode/GSE50521/results/GSE50521_objects.Rdata")


plot_dir <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/explore/select_low_exp_threshold/"
dir.create(plot_dir, showWarnings=F)
width = 1000
height = 800
df_gene<- read.table(header=TRUE, text='
                     GenderGenes    GeneSymbols
                     F    Xist
                     M    Kdm5d
                     M    Rps4y1')
### loop for all datasets and all methods
df_all <- NULL
for (r_obj in r_obj_ls){
    df <- mainExpThreshold (method_t_ls, r_obj, df_gene,
                            plot_dir, width=1000, height=800)
    df_all <- rbind(df_all, df)
    write.table(df_all, file = paste0(plot_dir, "stat_thresholds.tsv"), quote = F, sep ='\t',
                row.names =F)
}
