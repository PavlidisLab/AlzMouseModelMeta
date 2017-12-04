#' @note: the input Rdata has been outlier removed, subsetted, and batch-corrected to filter out low exp probes

#' 2016-04-11: add GSE52022
#' @todo
#' 2016-04-11: need to update emineJ to include background for biological process GO terms only


#---------------------------------------------------------------------------#
# PART 1: filter out low exp probes (prepare for DE)
#---------------------------------------------------------------------------#
## wrapper for filter low exp genes

 
source("./low_exp_rm_analysis/filter_low_exp_and_save.R")

r_obj_ls <- c("/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE48622/results/GSE48622_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE52022/results/GSE52022_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.1/results/GSE63617.1_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE63617.2/results/GSE63617.2_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE36237/results/GSE36237_objects_after_batch_correction.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/analysis_mode/GSE14499/results/GSE14499_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/analysis_mode/GSE1556/results/GSE1556_objects.Rdata",
              "/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/analysis_mode/GSE50521/results/GSE50521_objects.Rdata")

for(i in 1: length(r_obj_ls)){
    print(i)
    r_obj <- r_obj_ls[i]
    dataset <- getGSEID(r_obj)
    print(dataset)
    out_folder <- paste0("/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/explore/analysis_mode/", dataset, "/")
    filterExpArray(r_obj, out_folder)
}




#---------------------------------------------------------------------------#
# PART 2.1 LIMMA DEA
#---------------------------------------------------------------------------#
#-----------------------------#
# PART 2.1A: LIMMA   - gemma
#-----------------------------#

### limma for meta analysis only -gemma
rm(list=setdiff(ls(),'home_dir'))
q_threshold <- 0.05
datadir <- paste0('/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/', Sys.Date(), '/')
dir.create(datadir, showWarnings=F)
object_dir <- '/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/explore/analysis_mode/'
f_model <- '/home/bzhuang/AD_mouse_model_project/config_files/limma_prioritized_mouse_datasets_meta_analysis.tsv'
start_row <- 1
dataset_todo <- c('GSE50521', 'GSE1556')
source('limma_DEA.R')


#-----------------------------#
# PART 2.1B: LIMMA   - CEL
#-----------------------------#
### limma for meta analysis only -CEL
rm(list=setdiff(ls(),'home_dir'))
q_threshold <- 0.05
datadir <- paste0('/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/', Sys.Date(), '/')
dir.create(datadir, showWarnings=F)
object_dir <- '/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/explore/analysis_mode/'
f_model <- '/home/bzhuang/AD_mouse_model_project/config_files/limma_prioritized_CEL_mouse_datasets_meta_analysis.tsv'
start_row <- 1
#to_do <-2
dataset_todo <-c('GSE14499', 'GSE48622', 'GSE63617.1', 'GSE63617.2', 'GSE36237') # all hippo and AD datasets
#dataset_todo <-c('GSE52022')
source('limma_DEA.R')


#---------------------------------------------------------------------------#
# PART 2.2 LIMMA DEA result summary
#---------------------------------------------------------------------------#
# load the top table files of gender gene filtered
# known modifiers are the list from Joerg
# 2016-03-07, files must be in the same order as file labels
rm(list=setdiff(ls(),'home_dir'))
 
source("./meta_analysis/summary_DE_results.R")
files <- c("/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE48622/no_subset/results/GSE48622_limma_toptable_GenotypeAPLP2_KO-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE48622/no_subset/results/GSE48622_limma_toptable_GenotypeAPP_KO-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-04-11/GSE52022/subset/Timepoint_4_months/results/GSE52022_limma_toptable_Genotype5XFAD-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.1/subset/Timepoint_3_months_OrganismPart_Hippocampus/results/GSE63617.1_limma_toptable_GenotypeAD11-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE36237/batch_corrected_subset/Treatment_reference/results/GSE36237_limma_toptable_GenotypeTg2576-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/results/GSE63617.2_limma_toptable_GenotypeAD11-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE14499/subset/Treatment_GFP/results//GSE14499_limma_toptable_GenotypeJ20-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-04-11/GSE52022/subset/Timepoint_8_months/results/GSE52022_limma_toptable_Genotype5XFAD-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE1556/no_subset/results/GSE1556_limma_toptable_GenotypeTg2576-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE50521/no_subset/results/GSE50521_limma_toptable_Genotype5xFAD-vs-baseline_WT.tsv",
           "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/results/GSE63617.2_limma_toptable_GenotypeAD11-vs-baseline_WT.tsv")

file_labels <-c("GSE48622(APLP2 KO)_2m",
                "GSE48622(APP KO)_2m",
                "GSE63617_3m",
                "GSE52022_4m",
                "GSE36237_5m",
                "GSE63617_6m",
                "GSE14499_7m",
                "GSE52022_8m",
                "GSE1556_12m",
                "GSE50521_14-15m",
                "GSE63617_15m")
(time_labels <- c(rep("early", 5), rep("late", 6)))

f_out <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/results_lrm/limma_DE_summary/"
known_modifier <- "/home/bzhuang/neurogem/neurogem_modifiers_mouse_AD.tsv"
df <- mainSummaryDE(files, f_out,file_labels, time_labels, df_return=T, known_modifier = known_modifier)






#---------------------------------------------------------------------------#
# PART 3.1a: META ANALYSIS and JACKKNIFE (rm all NAs)
#---------------------------------------------------------------------------#
# need to copy the up/down regulated gene files to "AD_early_1_5_months","AD_med_5_7_months", "AD_late_6_15_months" folders in folder_pre for input
# (GSE50521 combined 5XFAD_elf2a, 5XFAD and wt and elf2a)
### meta analysis for up and down CEL and gemma (loop thru all folders)
rm(list=setdiff(ls(),'home_dir'))
source('meta.R')
adj_method='BH'
meta_method <-'Fisher'
j_threshold <- 0.1
gene_annotation='/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/gene_annotations.tsv'

folder_pre <- paste0('/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/')

folders <- c("AD_early_1_5_months", "AD_med_5_7_months", "AD_late_6_15_months")
for (i in folders){
    print(i)
    (f_meta <- paste0(folder_pre, i,'/'))
    (f_out <- paste0(folder_pre, i,'_'))
    
    for (keyword in c('up_regulation','down_regulation')){
        (meta_output <- paste0(f_out, keyword, '_meta_genes.tsv'))
        (j_output <- paste0(f_out, keyword, '_jackknife.tsv'))
        (file_list <- getFiles(f_meta, keyword))
        ## write the meta-signature genes (all files)
        metaPval(file_list, gene_annotation=gene_annotation, output = meta_output, method = meta_method, adj_method = adj_method, output_df =F,NA_rm =T, write_df = T)
        ## use the jackknife to output core gene
        jackknifeProcedure(file_list, gene_annotation=gene_annotation, threshold = j_threshold, output=j_output, method = meta_method, adj_method = adj_method, 
                           output_df =F, NA_rm =T, write_df = T)
    }
}


#---------------------------------------------------------------------------#
# PART 3.1b: META ANALYSIS and JACKKNIFE (rm all NAs) and rm GSE1556, only late phase
#---------------------------------------------------------------------------#
# need to copy the up/down regulated gene files to "AD_early_1_5_months", "AD_late_6_15_months" folders in folder_pre for input
### meta analysis for up and down CEL and gemma (loop thru all folders)
rm(list=setdiff(ls(),'home_dir'))
source('meta.R')
adj_method='BH'
meta_method <-'Fisher'
j_threshold <- 0.1
gene_annotation='/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/gene_annotations.tsv'

folder_pre <- paste0('/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/rm_gse1556/')

folders <- c("AD_late_6_15_months")
for (i in folders){
    print(i)
    (f_meta <- paste0(folder_pre, i,'/'))
    (f_out <- paste0(folder_pre, i,'_'))
    
    for (keyword in c('up_regulation','down_regulation')){
        (meta_output <- paste0(f_out, keyword, '_meta_genes.tsv'))
        (j_output <- paste0(f_out, keyword, '_jackknife.tsv'))
        (file_list <- getFiles(f_meta, keyword))
        ## write the meta-signature genes (all files)
        metaPval(file_list, gene_annotation=gene_annotation, output = meta_output, method = meta_method, adj_method = adj_method, output_df =F,NA_rm =T, write_df = T)
        ## use the jackknife to output core gene
        jackknifeProcedure(file_list, gene_annotation=gene_annotation, threshold = j_threshold, output=j_output, method = meta_method, adj_method = adj_method, 
                           output_df =F, NA_rm =T, write_df = T)
    }
}


#---------------------------------------------------------------------------#
# PART 3.2: Summary for META ANALYSIS and JACKKNIFE
#---------------------------------------------------------------------------#
rm(list=setdiff(ls(),'home_dir'))
 
source("./meta_analysis/summary_meta_jack_results.R")

known_modifier <- "/home/bzhuang/neurogem/neurogem_modifiers_mouse_AD.tsv"
phase_ls = c("early", "late", "med")
## final jackkinfe and meta tables(filtered by gender genes and NA removed)
file_folder <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/"
f_out_dir <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/results_lrm/meta_analysis_summary/"
dir.create(f_out_dir, showWarnings=F)
prefix <- "low_exp_filtered_"
mainSummaryMetaJack(file_folder, f_out_dir, prefix, phase_ls, known_modifier, venn_plot=T)

#---------------------------------------------------------------------------#
# PART 3.3: correlation heatmap for META ANALYSIS and JACKKNIFE
#---------------------------------------------------------------------------#
rm(list=setdiff(ls(),'home_dir'))
 
f_df <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/"
source("./meta_analysis/meta_spearman_cor_heatmap.R")





#---------------------------------------------------------------------------#
# PART 4: Plot significant jackknife genes
#---------------------------------------------------------------------------#
rm(list=setdiff(ls(),'home_dir'))
setwd("/home/bzhuang/git/bin/mouse_dataset_process")
source("meta_analysis/plot_sig_genes_helper.R")

file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE14499/", 
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-04-11/GSE52022/subset/Timepoint_4_months/",
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE36237/", 
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE48622/",
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE1556/", 
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.1/", 
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/limma_DE_lrm/2016-03-15/GSE50521/" )

jack_dir <-"/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/"
plot_out <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/results_lrm/Gene_expression_of_selected_genes/"
wrapperPlotGeneExp(jack_dir, plot_out, file_dir_ls, plt_jack =T, plt_cell =F, plt_gender =F,
                   phase_ls = c("early", "late", "med"), return_one_panel =T)



#---------------------------------------------------------------------------#
# PART 5: CREATE BACKGROUND FOR ERMINEJ, meta genes, jack genes
#---------------------------------------------------------------------------#
#see erminej_preprocess


#**************************
## 5A. background with all GO terms
#**************************

rm(list=setdiff(ls(),'home_dir'))
source('ermineJ_preprocess/background_for_pre_filtered_genes.R')
folder_in <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/"
fo_out <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/ermineJ_lrm/background/"
(meta_jack<- paste0(folder_in, grep("regulation.*.tsv", list.files(folder_in), value=T)))


# get the full annotation file with GO terms
folder_anno <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/"
(f_anno_ls <- grep("ermineJ.*.tsv", list.files(folder_anno, full.names=T), value=T))

for(keyword_p in c("early", "med", "late")){
    print(keyword_p)
    for(keyword_m in c("meta", "jackknife")){
        (f_anno <- grep(keyword_p, f_anno_ls, value=T))
        (f_up <- grep(paste0(keyword_p, ".*up.*", keyword_m), meta_jack, value=T))
        (f_down <- grep(paste0(keyword_p, ".*down.*", keyword_m), meta_jack, value=T))
        (f_out<- paste0(fo_out, "ermineJ_", keyword_p, "_", keyword_m,"_gene_anno_filtered.tsv"))
        getGOterms(f_anno, f_up, f_down, f_out, notes =paste0("for ", keyword_p, " ", keyword_m, " genes."))
    }
}


#**************************
## 5B. background with BP GO terms only
#**************************


