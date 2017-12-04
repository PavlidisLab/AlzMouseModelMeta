#' 2016-03-07
#' updated 2016-4-11
#' This is the wrapper script for non-filtered results and low exp are filtered by hard threshold for ermineJ enrichment
#' output are filtered meta and jackknife files and corresponding backgrounds
#' wrapper for get_all_gene_annotation.R (to aggregate all gene terms for specified platforms)
## wrapper for filter_out_low_exp_for_enrichment, this is for DE of non-filtered probes
 
source("./ermineJ_preprocess/get_all_gene_annotation.R")
source("./ermineJ_preprocess//filter_out_low_exp_for_enrichment.R")




#########################################################
## 1. Get gene annotation for specified platforms (all GO terms)
#########################################################



##************
## 1B. for ermineJ format for early samples
##*************
f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_early_gene_annotations.tsv'
platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/'
platform_ls_early <- c("GPL1261", "GPL7042")
df <- geneAnno(platform_ls_early, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)

##************
## 1C. for ermineJ format for late samples
##*************
platform_ls <- c("GPL1261", "GPL7042", "GPL6096", "GPL81")
f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_late_gene_annotations.tsv'
platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/'
df <- geneAnno(platform_ls, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)

##************
## 1D. for ermineJ format for lmed samples (5-7 months)
##*************
platform_ls <- c("GPL1261", "GPL7042")
f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_med_gene_annotations.tsv'
platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/'
df <- geneAnno(platform_ls, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)



#########################################################
## 2. Get gene annotation for specified platforms (for BP GO terms only)
#########################################################

##************
## 2A. for ermineJ format for early samples, 
##*************
 
source("./ermineJ_preprocess/get_all_gene_annotation.R")
source("./ermineJ_preprocess//filter_out_low_exp_for_enrichment.R")
f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/ermineJ_early_gene_annotations_BP_GO_only.tsv'
platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/'
platform_ls_early <- c("GPL1261", "GPL7042")
df <- geneAnno(platform_ls_early, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)

##************
## 2B. for ermineJ format for late samples
##*************
# platform_ls <- c("GPL1261", "GPL7042", "GPL6096", "GPL81")
# f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/ermineJ_late_gene_annotations_BP_GO_only.tsv'
# platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/'
# df <- geneAnno(platform_ls, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)

##************
## 2C. for ermineJ format for lmed samples (5-7 months)
##*************
platform_ls <- c("GPL1261", "GPL7042")
f_erminej_anno<-'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/ermineJ_med_gene_annotations_BP_GO_only.tsv'
platform_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/'
df2 <- geneAnno(platform_ls, platform_folder, f_erminej_anno, ermineJ_format = T, df_return = T)






###############################################################
## 3. meta genes and jackknife genes, rm low expr genes by hard threshold (input_tables/low_exp_genes_filtered/)
###############################################################
#**************************
## 3A. for early genes
#**************************

fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_early_gene_annotations_low_expr_filtered.tsv"

fo_folder <- paste0("/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/input_tables/low_exp_genes_filtered/")
dir.create(fo_folder, recursive=T)
f_meta = c("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_early_1_5_months_up_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/rm_KO_AD_early_1_5_months_down_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_early_1_5_months_up_regulation_jackknife.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_early_1_5_months_down_regulation_jackknife.tsv")
fo_meta = c(paste0(fo_folder, "AD_early_1_5_months_up_regulation_meta_genes.tsv"),
            paste0(fo_folder, "rm_KO_AD_early_1_5_months_down_regulation_meta_genes.tsv"),
            paste0(fo_folder, "AD_early_1_5_months_up_regulation_jackknife.tsv"),
            paste0(fo_folder, "AD_early_1_5_months_down_regulation_jackknife.tsv"))


FilterLowExprGenesMeta(fo_erminej_bg_filter,
                       f_meta = f_meta, fo_meta =fo_meta,
                       rm_meta_NA = T)

#**************************
## 3B. for late genes
#**************************
fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_late_gene_annotations_low_expr_filtered.tsv"

fo_folder <- paste0("/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/input_tables/low_exp_genes_filtered/")
dir.create(fo_folder, recursive=T)
f_meta = c("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_late_6_15_months_up_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_late_6_15_months_up_regulation_jackknife.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_late_6_15_months_down_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_late_6_15_months_down_regulation_jackknife.tsv")

fo_meta = c(paste0(fo_folder, "AD_late_6_15_months_up_regulation_meta_genes.tsv"),
            paste0(fo_folder, "AD_late_6_15_months_up_regulation_jackknife.tsv"),
            paste0(fo_folder, "AD_late_6_15_months_down_regulation_meta_genes.tsv"),
            paste0(fo_folder, "AD_late_6_15_months_down_regulation_jackknife.tsv"))

FilterLowExprGenesMeta(fo_erminej_bg_filter,
                       f_meta = f_meta, fo_meta =fo_meta,
                       rm_meta_NA = T)



#**************************
## 3C. for med genes
#**************************

fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_med_gene_annotations_low_expr_filtered.tsv"

fo_folder <- paste0("/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/input_tables/low_exp_genes_filtered/")
dir.create(fo_folder, recursive=T)
f_meta = c("/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_med_5_7_months_up_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_med_5_7_months_up_regulation_jackknife.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_med_5_7_months_down_regulation_meta_genes.tsv",
           "/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-03-07/AD_med_5_7_months_down_regulation_jackknife.tsv")

fo_meta = c(paste0(fo_folder, "AD_med_5_7_months_up_regulation_meta_genes.tsv"),
            paste0(fo_folder, "AD_med_5_7_months_up_regulation_jackknife.tsv"),
            paste0(fo_folder, "AD_med_5_7_months_down_regulation_meta_genes.tsv"),
            paste0(fo_folder, "AD_med_5_7_months_down_regulation_jackknife.tsv"))

FilterLowExprGenesMeta(fo_erminej_bg_filter,
                       f_meta = f_meta, fo_meta =fo_meta,
                       rm_meta_NA = T)



###############################################################
## 4. create ermineJ background sets for jackknife and meta with filtered input
###############################################################
#**************************
## 4A. background with all GO terms
#**************************
#see erminej_preprocess
rm(list=setdiff(ls(),'home_dir'))
source('ermineJ_preprocess/background_for_pre_filtered_genes.R')
folder_in <- "/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/input_tables/low_exp_genes_filtered/"
fo_out <- "/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/background/"
dir.create(fo_out)
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
## 4B. background with BP GO terms only
#**************************
#see erminej_preprocess
rm(list=setdiff(ls(),'home_dir'))
source('ermineJ_preprocess/background_for_pre_filtered_genes.R')
folder_in <- "/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/input_tables/low_exp_genes_filtered/"
fo_out <- "/home/bzhuang/AD_mouse_model_project/ermineJ/2016-03-07/background_BP_only/"
dir.create(fo_out)
(meta_jack<- paste0(folder_in, grep("regulation.*.tsv", list.files(folder_in), value=T)))


# get the full annotation file with BP GO terms
folder_anno <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/biological_pathway_only/"
(f_anno_ls <- grep("ermineJ.*.tsv", list.files(folder_anno, full.names=T), value=T))

for(keyword_p in c("early", "med")){
    print(keyword_p)
    for(keyword_m in c("meta", "jackknife")){
        (f_anno <- grep(keyword_p, f_anno_ls, value=T))
        (f_up <- grep(paste0(keyword_p, ".*up.*", keyword_m), meta_jack, value=T))
        (f_down <- grep(paste0(keyword_p, ".*down.*", keyword_m), meta_jack, value=T))
        (f_out<- paste0(fo_out, "ermineJ_", keyword_p, "_", keyword_m,"_gene_anno_filtered_BP_only.tsv"))
        getGOterms(f_anno, f_up, f_down, f_out, notes =paste0("for ", keyword_p, " ", keyword_m, " genes."))
    }
}


#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
#*************************************************************
# 
# ## this is the old way for background filtereing, not in use 
# ###############################################################
# ## NA. early and late, med backgrounds, low expr filtered by hard threshold
# ###############################################################
# #**************************
# ## for early genes
# #**************************
# file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE48622/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE36237/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.1/")
# 
# f_erminej_anno <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_early_gene_annotations.tsv"
# fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_early_gene_annotations_low_expr_filtered.tsv"
# 
# df_filtered <- mainFilterLowExprGenesBackground(file_dir_ls, f_erminej_anno, fo_erminej_bg_filter,filter_method = "median",
#                                                 rm_multi_gene =T, 
#                                                 affy_threshold = 6, agilent_threshold = 0,
#                                                 df_return = T)
# #**************************
# ## for late genes
# #**************************
# f_erminej_anno <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_late_gene_annotations.tsv"
# fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_late_gene_annotations_low_expr_filtered.tsv"
# 
# file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE14499/", 
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-02-15/GSE1556/", 
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-03-07//GSE50521/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/")
# 
# 
# df_filtered_late <- mainFilterLowExprGenesBackground(file_dir_ls, f_erminej_anno, fo_erminej_bg_filter,filter_method = "median",
#                                                      rm_multi_gene =T, 
#                                                      affy_threshold = 6, agilent_threshold = 0,
#                                                      df_return = T)
# 
# #**************************
# ## for med genes (5-7 months), 2016-03-16
# #**************************
# f_erminej_anno <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/platforms/ermineJ_med_gene_annotations.tsv"
# fo_erminej_bg_filter <- "/home/bzhuang/AD_mouse_model_project/ermineJ/background/ermineJ_med_gene_annotations_low_expr_filtered.tsv"
# 
# file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE14499/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE36237/",
#                  "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/")
# 
# 
# df_filtered_late <- mainFilterLowExprGenesBackground(file_dir_ls, f_erminej_anno, fo_erminej_bg_filter,filter_method = "median",
#                                                      rm_multi_gene =T, 
#                                                      affy_threshold = 6, agilent_threshold = 0,
#                                                      df_return = T)
