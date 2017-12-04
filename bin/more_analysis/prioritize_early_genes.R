# 2016-03-23
#' with the meta and jackknife genes (before filter) grouped by phases
#' 1. filter the lowly expressed genes (by hard threshold), and remove NA genes
#' 2. from meta genes, find the genes have consistent low raw p value ( <0.1 across all studies in the phase) and plot these genes

#---------------------------------------------------------------------------#
# PART 1: REMOVE LOWLY EXPRESSED GENES
#---------------------------------------------------------------------------#
## rm low expr genes
rm(list=setdiff(ls(),'home_dir'))
 
source("./ermineJ_preprocess/get_all_gene_annotation.R")
source("./ermineJ_preprocess//filter_out_low_exp_for_enrichment.R")
source('ermineJ_preprocess/background_for_pre_filtered_genes.R')

fo_folder <- paste0("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma//2016-03-23/low_exp_rm/")
dir.create(fo_folder, recursive=T)
keyword_ls <- c("early.*.meta.*.tsv", "early.*.jackknife.*.tsv", "med.*.meta.*.tsv", "med.*.jackknife.*.tsv","late.*.meta.*.tsv", "late.*.jackknife.*.tsv")
for (keyword in keyword_ls){
    fo_erminej_bg_filter <- grep(keyword, list.files("/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/background/", full.names=T), value=T)
    (f_meta = grep(keyword, list.files("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/", full.names=T), value=T))
    (fo_meta = grep(keyword, list.files("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/"), value=T))
    (fo_meta = paste0(fo_folder, fo_meta))
    FilterLowExprGenesMeta(fo_erminej_bg_filter,
                           f_meta = f_meta, fo_meta =fo_meta,
                           rm_meta_NA = T)
}


#---------------------------------------------------------------------------#
# PART 2: Function to find low p genes and get jackknife rank
#---------------------------------------------------------------------------#
#' find meta genes (df1) that have consistent low p values ( <0.1), and get the jackknife ranks too
getPrioritizedGenes <- function(df1, df_j, f_out_up){
    df_j <- read.delim(df_j, comment.char="#")
    df_j$jackknife_rank <- 1:nrow(df_j)
    df_j <- df_j[, c("geneSymbol", "jackknife_rank")]
    df <- read.delim(df1, comment.char="#")
    max_index <- apply(df[5:7], 1, max) 
    df <- left_join(df, df_j)
    index <- which(max_index < 0.1)
    df<- df[index, ]%>%droplevels
    
    df <- df[order(df$jackknife_rank), ]
    write.table(df, file =f_out_up, quote = F, sep ='\t',
                row.names =F)
}

#**********************#
#**2A. for intermediate phase 5-7m
#**********************#
dir_out_up <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/'
dir.create(dir_out_up)
df1 <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_med_5_7_months_up_regulation_meta_genes.tsv'
df_j <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_med_5_7_months_up_regulation_jackknife.tsv'
f_out_up <- paste0(dir_out_up, 'med_up.tsv')
getPrioritizedGenes(df1, df_j, f_out_up)



df1 <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_med_5_7_months_down_regulation_meta_genes.tsv'
df_j <-'/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_med_5_7_months_down_regulation_jackknife.tsv'
f_out_down <-  paste0(dir_out_up,'med_down.tsv')
getPrioritizedGenes(df1, df_j, f_out_down)

setwd("/home/bzhuang/git/bin/mouse_dataset_process")
source('meta_analysis/plot_sig_genes_helper.R')

file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE14499/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE36237/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE48622/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-02-15/GSE1556/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.1/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-03-07/GSE50521/" )

plot_out_dir <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/plots_med_up/'
f_gene_list <- f_out_up
plotExpWithGeneList(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel=T, with_gene_order=T)

plot_out_dir <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/plots_med_down/'
f_gene_list <- f_out_down
plotExpWithGeneList(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel=T, with_gene_order=T)



#**********************#
#**2A. for early phase 1-5m
#**********************#
dir_out_up <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/'
dir.create(dir_out_up)

df1 <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_early_1_5_months_up_regulation_meta_genes.tsv'
df_j <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_early_1_5_months_up_regulation_jackknife.tsv'
f_out_up <- paste0(dir_out_up, 'early_up.tsv')
getPrioritizedGenes(df1, df_j, f_out_up)


df1 <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_early_1_5_months_down_regulation_meta_genes.tsv'
df_j <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/AD_early_1_5_months_down_regulation_jackknife.tsv'
f_out_down <-  paste0(dir_out_up,'early_down.tsv')
getPrioritizedGenes(df1, df_j, f_out_down)



setwd("/home/bzhuang/git/bin/mouse_dataset_process")
source('meta_analysis/plot_sig_genes_helper.R')

file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE14499/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE36237/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE48622/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-02-15/GSE1556/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.1/", 
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/",
                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-03-07/GSE50521/" )

plot_out_dir <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/plots_early_up/'
f_gene_list <- f_out_up
plotExpWithGeneList(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel=T, with_gene_order=T)

plot_out_dir <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-23/low_exp_rm/prioritized_genes/plots_early_down/'
f_gene_list <- f_out_down
plotExpWithGeneList(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel=T, with_gene_order=T)



