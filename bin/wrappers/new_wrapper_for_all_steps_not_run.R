cat("
    #---------------------------------------------------------------------------#
    # PART 3.2 LIMMA DEA result summary, and removed samples before limma
    #---------------------------------------------------------------------------#\n")

## load the top table files of non-filtered
# known modifiers are the list from Joerg
# 2016-03-07, files must be in the same order as file labels
# 2016-04-05 add GSE52022 4m samples
rm(list=setdiff(ls(),'home_dir'))


source("./meta_analysis/summary_DE_results.R")
source('meta_analysis/summary_removed_samples.R')

summaryTempFun <- function(df_info, f_out, exclude_ls=NULL){
    ## set default files
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma toptables
    
    (files <- grep('toptable_Genotype', list.files(datadir, recursive=T, full.names=T), value=T))
    if(!is.null(exclude_ls)){
        files <- grep(paste0(exclude_ls, collapse='|'), files, invert = T, value = T)
    }
    known_modifier <- known_modifiers_config
    
    cat("
        ##************
        ## PART 3.2.1 get the bar plot DE summary, limme DE counts, only get the ones in dataset info 
        ##************")
    df <- mainSummaryDE(files, f_out, df_info, df_return=T, known_modifier = known_modifier, bar_text=4,match_all =F)
    
    cat("
        ##************
        ## PART 3.2.2 save the input summary and combine with limma DE summary, and the DE ratio
        ##************")
    (f_DE_summary <- max(grep('limma_DE_FDR_counts', list.files(f_out, full.names=T), value = T))) ## pick the most recent one
    getSampleSummary(datadir, f_out, f_DE_summary=f_DE_summary,plot_DE=T)
    
    cat("
        ##************
        ## PART 3.2.3: Summary of removed samples due to outlier, or subset
        ##************")
    data_folder <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
    f_o<- paste0(f_out,'/removed_samples.tsv') 
    summaryOutlierResult(data_folder, f_o)
}

source('config_wrappers.R')
for(disease in disease_ls){
    
    source('config_wrappers.R')
    
    ## run function for limma models only
    exclude_ls <- exclude_genotype_config    ## genotypes, and others to be excluded from the summary
    #df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_all.tsv')
    df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  #dataset label, phase label, order of the datasets etc.
    f_out<- paste0(disease_dir, 'results/limma_DE_summary/', Sys.Date(), '/selected_genotypes/') ## result and plot outdir
    
    summaryTempFun(df_info, f_out, exclude_ls)
    
    
    ## run function for mixed models genotypes
    exclude_ls <- NULL    ## include all mixed model genotypes
    df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  #dataset label, phase label, order of the datasets etc.
    f_out<- paste0(disease_dir, 'results/limma_DE_summary/', Sys.Date(), '/all_genotypes/') ## result and plot outdir
    
    summaryTempFun(df_info, f_out, exclude_ls)
}



cat("
    #---------------------------------------------------------------------------#
    # PART 3.3.1: Cell population estimation from limma object dir-- on local dirs
    #---------------------------------------------------------------------------#\n")
## the estimate.R file from Ogan is modified so that tables can be exported
rm(list=setdiff(ls(),'home_dir'))
disease='AD'

source('config_wrappers.R')
source("MGP_estimation/estimate_cell_population.R")

file_ls <- paste0(limma_dir, 'prioritized_limma/')
all_dataset <- F  # if to do all datasets, if not define dataset_todo to specify which dataset to do


(r_obj_ls <- list.files(file_ls, full.names=T, recursive=T, pattern= ".Rdata"))
if(!all_dataset){
    dataset_todo <-c('GSE63617.1', 'GSE63617.2')
    r_obj_ls <- grep(paste0(dataset_todo, collapse = '|'), r_obj_ls, value = T)
}


#(r_obj_ls <- r_obj_ls[2]) 
#' have trouble with GSE1556 (AD), GSE31458(PD), GSE8030_8_weeks(PD)
#' trouble : Error in prcomp.default(t(relevantExpr[groups %in% unique(groups)[j]]),  : 
#' cannot rescale a constant/zero column to unit variance 
#(r_obj_ls <- r_obj_ls[c(5:14)]) 

(outDir <-paste0(disease_dir, '/results/Cell_population_estimates_',Sys.Date(), '/'))
dir.create(outDir,showWarnings=F, recursive=T)
genes <- cell_pop_marker_genes_config
genes <- puristOut(genes)

for (i in 1: length(r_obj_ls)){
    (r_ob <- r_obj_ls[i])
    print(r_ob)
    cellPopEstimate(r_ob, genes, outDir)
}


cat("
    #---------------------------------------------------------------------------#
    # PART 3.3.2: Cell population estimation from limma object dir-- on local
    #' for individual studies
    #---------------------------------------------------------------------------#\n")
## recorded done the cell type markers for future heatmaps, save in git/doc/
source('cell_population_markers.R')

cat("
    ##************
    ## PART 3.3.2a  Cell population estimation ## use genotype as factor, p value corrected by 'BH'
    ##************")

rm(list=setdiff(ls(),'home_dir'))

source("MGP_estimation/estimate_cell_population.R")
disease_ls=c('AD', 'HD', 'PD')


all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Genotype'  # design_group = 'Original_genotype'
plot_box = T  ## plot violin plots
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT

threshold = 6 # filter probes exp less than threshold


removeNegatives = F  ## ogan suggested use F, both give similar results

for (disease in disease_ls){
    
    source('config_wrappers.R')
    
    file_ls <- paste0(limma_dir, 'prioritized_limma/')
    
    (r_obj_ls <- list.files(file_ls, full.names=T, recursive=T, pattern= ".Rdata"))
    
    # if(!all_dataset){
    #     dataset_todo <-c('GSE26317', 'GSE10202', 'GSE48104')
    #     r_obj_ls <- grep(paste0(dataset_todo, collapse = '|'), r_obj_ls, value = T)
    # }
    
    
    mixed_model =F
    (outDir <-paste0(disease_dir, '/results/Cell_population_estimates_', design_group,'/',Sys.Date(), '/'))
    dir.create(outDir,showWarnings=F, recursive=T)
    #     genes <- cell_pop_marker_genes_config
    #     genes <- puristOut(genes)
    
    for (i in 1: length(r_obj_ls)){
        (r_ob <- r_obj_ls[i])
        print(r_ob)
        cellPopEstimatePlot(r_ob, outDir, design_group = design_group, mixed_model = mixed_model, plot_box =plot_box,
                            sigTest = sigTest, 
                            wt_only = wt_only,
                            removeNegatives = removeNegatives, threshold = threshold)
    } 
}




cat("
    ##************
    ## PART 3.3.2b  Cell population estimation ## use original genotype as factor
    # for individual studies - on local
    ##************")
# out put in diseaese/results/Cell_population_estimates

rm(list=setdiff(ls(),'home_dir'))
disease_ls=c('AD', 'HD', 'PD')



source("MGP_estimation/estimate_cell_population.R")

all_dataset <- F  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Original_genotype'  # design_group = 'Original_genotype'
plot_box =T
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT
mixed_model =F ## F: not expression from mixed model
threshold = 6 # filter probes exp less than threshold

removeNegatives = F

for(disease in disease_ls){
    source('config_wrappers.R')
    
    file_ls <- paste0(limma_dir, 'prioritized_limma/')
    (r_obj_ls <- list.files(file_ls, full.names=T, recursive=T, pattern= ".Rdata"))
    
    #     if(!all_dataset){
    #         dataset_todo <-c('GSE48622', 'GSE64398.1')
    #         r_obj_ls <- grep(paste0(dataset_todo, collapse = '|'), r_obj_ls, value = T)
    #     }
    
    (outDir <-paste0(disease_dir, '/results/Cell_population_estimates_', design_group,'/',Sys.Date(), '/'))
    dir.create(outDir,showWarnings=F, recursive=T)
    genes <- cell_pop_marker_genes_config
    genes <- puristOut(genes)
    
    for (i in 1: length(r_obj_ls)){
        (r_ob <- r_obj_ls[i])
        print(r_ob)
        cellPopEstimatePlot(r_ob, outDir, design_group = design_group, mixed_model =F, plot_box =plot_box,
                            sigTest = sigTest, 
                            wt_only = wt_only,
                            removeNegatives = removeNegatives,
                            threshold= threshold)
    }
}




cat("
    #---------------------------------------------------------------------------#
    ### PART 3.3.2d  Cell population estimation from mixed model-- on local, 
    # input own markers (obj: genes)
    #---------------------------------------------------------------------------#\n")


rm(list=setdiff(ls(),'home_dir'))

source("MGP_estimation/estimate_cell_population.R")

disease_ls =c('AD')
phase_ls = c('early', 'late')
load('../../doc/cell_type_markers/Wang_et_al_TS3_genes.Rdata')


disease_ls =c('HD')
phase_ls = c('early', 'late')

all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Genotype'  # design_group = 'Original_genotype'
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT

threshold = -10 # filter probes exp less than threshold (input is study corrected value, and low expression filtered)
mixed_model =T
removeNegatives = F  ## ogan suggested use F, 


for (disease in disease_ls){
    for(phase in phase_ls){
        source('config_wrappers.R')
        (file_ls <- paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/',disease, '/', phase,'/'))
        
        
        (r_ob <- list.files(file_ls, full.names=T, recursive=T, pattern= "mixed_model_results_exp_corrected.Rdata"))
        (outDir  <- paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/Wang_et_al_TS3/',disease, '/', phase,'/'))
        dir.create(outDir,showWarnings=F, recursive=T)
        
        cellPopEstimatePlot(r_ob, outDir, design_group = design_group, mixed_model = mixed_model,
                            sigTest = sigTest, 
                            wt_only = wt_only,
                            removeNegatives = removeNegatives, threshold = threshold, genes = genes)
    }
}




cat("
    #---------------------------------------------------------------------------#
    # PART 4.0: preperation for Fisher and JACKKNIFE
    #---------------------------------------------------------------------------#\n")

cat("
    ##************
    ## 4.0.1. make a big annotation with all platforms including GOterms, which platform, gemmaID, NCBIid etc
    ## only use the biological only GO terms
    ##*************")
rm(list=setdiff(ls(),'home_dir'))

source("./meta_analysis/summary_DE_results.R")
source("./ermineJ_preprocess/get_all_gene_annotation.R")

disease_ls <- c('AD', 'HD')

for (disease in disease_ls){
    source('config_wrappers.R')
    df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')
    df_info <- read.delim(df_info, comment.char="#")
    f_anno <- paste0(disease_dir, 'data_and_QC/design_and_annotation/annotations/') # the output folder
    dir.create(f_anno, recursive=T, showWarnings=F)
    f_erminej_anno <- paste0(f_anno, '/gene_annotations.tsv')
    #platform_folder <- paste0(disease_dir, 'data_and_QC/design_and_annotation/platforms/')
    platform_folder <- paste0(disease_dir, 'data_and_QC/design_and_annotation/platforms/biological_pathway_only/')
    (platform_ls <- unique(as.character(df_info$Platform)))
    df <- geneAnno(platform_ls, platform_folder, f_erminej_anno, df_return = T)
}


cat("
    ##************
    ## 4.0.2 link gene p values files from limma DE to the meta and jack folder by disease phase
    ##*************")
rm(list=setdiff(ls(),'home_dir'))


source("./meta_analysis/summary_DE_results.R") ## for getDatasetLabels
info_df <- 'dataset_info_mixed_model.tsv'  ## which dataset info to match
match_all = F ### if T, p value files must match exactly the same as listed in dataset info, if F, only select the ones matched with dataset info

disease_ls <- c('AD', 'HD', 'PD')

for (disease in disease_ls){
    source('config_wrappers.R')
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma toptables
    exclude_ls <- exclude_genotype_config    ## genotypes, and others to be excluded from the summary
    df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/', info_df)
    meta_dir <- paste0(disease_dir, 'meta_analysis/meta_jack/')
    dir.create(meta_dir, recursive=T, showWarnings=F)
    
    (files <- grep('gene_.*_regulation', list.files(datadir, recursive=T, full.names=T), value=T))
    if(!is.null(exclude_ls)){
        files <- grep(paste0(exclude_ls, collapse='|'), files, invert = T, value = T)
    }
    
    df <- getDatasetLabels(files, df_info,match_all = match_all)  #match_all = T
    
    for(phase in levels(df$Phase)){
        ## create folders of each phase
        sub_dir <- paste0(meta_dir, "/", phase, "/")
        dir.create(sub_dir, recursive=T, showWarnings=F)
        file_ls <- df$File[which(df$Phase == phase)]
        file.symlink(file_ls, to = sub_dir)
    }
}


cat("
    #---------------------------------------------------------------------------#
    # PART 4A: Fisher for early and late phases
    # all NAs are kept, lowly expressed genes are not filtered out
    #---------------------------------------------------------------------------#\n")

### Fisher for up and down CEL and gemma (loop thru all folders)
rm(list=setdiff(ls(),'home_dir'))


source('meta.R')
adj_method='BH'
meta_method <-'Fisher'
j_threshold <- 0.1
jack_rm_n = 1  ## remove 1 file a time


disease_ls <- c('AD', 'HD', 'PD')

for (disease in disease_ls){
    source('config_wrappers.R')
    gene_annotation=paste0(disease_dir, 'data_and_QC/design_and_annotation/annotations/gene_annotations.tsv')
    
    folder_pre <- paste0(disease_dir, 'meta_analysis/meta_jack/')
    
    (folders <- intersect(list.dirs(folder_pre, recursive=F, full.names=F), c("early", "late", "intermediate"))) # get the foldernames
    for (i in folders){
        print(i)
        (f_meta <- paste0(folder_pre, i,'/'))
        (f_out <- paste0(folder_pre, i,'_'))
        
        for (keyword in c('up_regulation','down_regulation')){
            (meta_output <- paste0(f_out, keyword, '_meta_genes.tsv'))
            (j_output <- paste0(f_out, keyword, '_jackknife.tsv'))
            (file_list <- getFiles(f_meta, keyword))
            ## write the meta-signature genes (all files)
            metaPval(file_list, gene_annotation=gene_annotation, output = meta_output, method = meta_method, adj_method = adj_method, output_df =F,NA_rm =F, write_df = T)
            ## use the jackknife to output core gene
            jackknifeProcedure(file_list, gene_annotation=gene_annotation, jack_rm_n = jack_rm_n, threshold = j_threshold, output=j_output, method = meta_method, adj_method = adj_method, 
                               output_df =F, NA_rm =F, write_df = T)
        }
    }
}




cat("
    #*****************************************************************************
    #---------------------------------------------------------------------------#
    # PART 4B.1: Summary for Fisher and JACKKNIFE
    # probes not filtered
    #---------------------------------------------------------------------------#
    #*****************************************************************************")
## final jackkinfe and meta tables(expression not filtered)

rm(list=setdiff(ls(),'home_dir'))


source('config_wrappers.R')
source("./meta_analysis/summary_meta_jack_results.R")

known_modifier <- known_modifiers_config

file_folder <- paste0(disease_dir, 'meta_analysis/meta_jack/') # the meta and jack files
f_out_dir <- paste0(disease_dir, 'results/meta_analysis_summary/')
prefix <- "not_filtered_"

dir.create(f_out_dir, showWarnings=F)
(phase_ls <- intersect(list.dirs(file_folder, recursive=F, full.names=F), c("early", "late", "intermediate")))
mainSummaryMetaJack(file_folder, f_out_dir, prefix, phase_ls, known_modifier, venn_plot=T)


cat("
    #---------------------------------------------------------------------------#
    # PART 4B.2: correlation heatmap for Fisher and JACKKNIFE
    # probes not filtered
    #---------------------------------------------------------------------------#\n")
rm(list=setdiff(ls(),'home_dir'))


source('config_wrappers.R')
f_df <- paste0(disease_dir, 'meta_analysis/meta_jack/')
source("./meta_analysis/meta_spearman_cor_heatmap.R")
metaHeatmap(f_df)


cat("
    #---------------------------------------------------------------------------#
    # PART 4B.3: Fisher and JACKKNIFE heatmap for top genes (all probes)
    #---------------------------------------------------------------------------#\n")

rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('helper_functions.R')
## a dir contains all the R objects saved after limma DE
r_object_dir <-paste0(limma_dir, 'prioritized_limma/')
df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv') # for early, intermediate and late files
## where the jack and meta files
(jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))

# plot out dir (a sub folder of 'gene_heatmaps' will be created)
plt_out <- paste0(disease_dir, '/meta_analysis/')

low_exp_rm <- F  # whether to remove lowly expressed probes
threshold <- 50  #top genes threshold
cluster_rows <- F  # for heatmap whether to cluster by rows/probes

source('meta_analysis/get_top_gene_heatmap.R')


cat("
    #---------------------------------------------------------------------------#
    # PART 4B.4: test if celltype markers are overrepresented or under represented 
    #         in the top meta/jack genes (all probes)
    #---------------------------------------------------------------------------#\n")
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('meta_analysis/compare_cell_marker_meta_jack.R')

out_dir <-paste0(disease_dir, '/meta_analysis/cell_marker_comparison/')
dir.create(out_dir,showWarnings=F, recursive=T)

(f_marker_genes <- cell_pop_marker_genes_config)
(jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))
## meta files
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'meta.*tsv', full.names = T))
threshold_method_ls <- c('top', 'percent', 'q_value')
threshold_ls <- c(100, 0.1, 0.1)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'meta_celltype_fisher_test_results.tsv', rm_na = T)

## jack files
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'jackknife.*tsv', full.names = T))
threshold_method_ls <- c('top', 'percent', 'q_value')
threshold_ls <- c(100, 0.1, 0.2)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'jack_celltype_fisher_test_results.tsv', rm_na = T)

## jack files, just the percent method only
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'jackknife.*tsv', full.names = T))
threshold_method_ls <- c('percent')
threshold_ls <- c(0.05)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'jack_celltype_fisher_test_results_top_5_percent.tsv', rm_na = T)


cat("
    #***************************************************************************#
    # PART 4B.5 Plot the meta/jack ranking against each study, calculate spearman correlation (all probes)
    #***************************************************************************#
    
    #' study (for each timepoint and genotype base on the dataset_info_mixed_model.tsv), plot the ranking of meta/jack of the genes
    #' against the ranking in the study (scatter plot with smoothing)
    #' 
    #' calculate the spearman correlation")
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('meta_analysis/compare_indi_to_meta_jack.R')

## plot out
(plot_out <- paste0(disease_dir, '/meta_analysis/ranking_comparisons/'))

## get the meta files and jack files and order them (so that meta and jack are correponding)
(jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))
(meta_f_ls <- list.files(jack_meta_folder, pattern = 'meta_genes.tsv', full.names = T))
(jack_f_ls <- list.files(jack_meta_folder, pattern = 'jackknife.tsv', full.names = T))
(meta_f_ls <- meta_f_ls[order(meta_f_ls)])
(jack_f_ls <- jack_f_ls[order(jack_f_ls)])


for (i in 1:length(meta_f_ls)){
    (f_meta <- meta_f_ls[i])
    (f_jack <- jack_f_ls[i])
    mainCompareRanksMetaIndi(f_meta = f_meta, f_jack=f_jack, plot_out=plot_out)
}


cat("
    #*****************************************************************************
    #---------------------------------------------------------------------------#
    # PART 4C.0: FILTER LOWLY EXPRESSED GENES METAANALYSIS AND JACKKNIFE
    # filter out lowly expressed probes (by hard threshold)
    #---------------------------------------------------------------------------#
    #*****************************************************************************")


rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source("./ermineJ_preprocess/filter_lowly_exp_meta_jack.R")

## meta, jack folders
(jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))
r_object_dir <-paste0(limma_dir, 'prioritized_limma/')# the limma object
df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv') # for early, intermediate and late files
out_dir <- paste0(disease_dir, '/meta_analysis/low_exp_rm/')
affy_threshold = 6
agilent_threshold = -10
filter_method = "median"

makeFilteredJackmeta(jack_meta_folder = jack_meta_folder, r_object_dir=r_object_dir, df_info=df_info, 
                     out_dir=out_dir, 
                     affy_threshold = affy_threshold, agilent_threshold = agilent_threshold, 
                     filter_method = filter_method)



cat("
    #***************************************************************************#
    # PART 4C.1: Summary for Fisher and JACKKNIFE
    # filtered lowly expressed probes
    #***************************************************************************#\n")
## final jackkinfe and meta tables(expression not filtered)

rm(list=setdiff(ls(),'home_dir'))


source('config_wrappers.R')
source("./meta_analysis/summary_meta_jack_results.R")

known_modifier <- known_modifiers_config

file_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/') # the meta and jack files (rm probes)
f_out_dir <- paste0(disease_dir, 'results/meta_analysis_summary/low_exp_rm/')
prefix <- "filtered_"

dir.create(f_out_dir, showWarnings=F)
(phase_ls <- intersect(unique(unlist(strsplit(list.files(file_folder, recursive=F, full.names=F), split='_'))), 
                       c("early", "late", "intermediate")))
mainSummaryMetaJack(file_folder, f_out_dir, prefix, phase_ls, known_modifier, venn_plot=T)


cat("
    #***************************************************************************#
    # PART 4C.2: correlation heatmap for Fisher and JACKKNIFE
    # filtered lowly expressed probes
    #***************************************************************************#\n")
rm(list=setdiff(ls(),'home_dir'))


source('config_wrappers.R')
f_df <- paste0(disease_dir, 'meta_analysis/low_exp_rm/')
source("./meta_analysis/meta_spearman_cor_heatmap.R")
metaHeatmap(f_df)


cat("
    #***************************************************************************#
    # PART 4C.3: Fisher and JACKKNIFE heatmap for top genes (filtered probes)
    #***************************************************************************#\n")

rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('helper_functions.R')

keyword_ls <- c('jackknife', 'meta_genes')
regulation_ls <- c('up', 'down')

## a dir contains all the R objects saved after limma DE
r_object_dir <-paste0(limma_dir, 'prioritized_limma/')
df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv') # for early, intermediate and late files
## where the jack and meta files
(jack_meta_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))

# plot out dir (a sub folder of 'gene_heatmaps' will be created)
plt_out <- paste0(disease_dir, '/meta_analysis/low_exp_rm/')

low_exp_rm <- F  # whether to remove lowly expressed probes # hear the low probes are already removed from the meta genes
threshold <- 50  #top genes threshold
cluster_rows <- F  # for heatmap whether to cluster by rows/probes

source('meta_analysis/get_top_gene_heatmap.R')


cat("
    #***************************************************************************#
    # PART 4C.4: test if celltype markers are overrepresented or under represented 
    #         in the top meta/jack genes (filtered probes)
    #***************************************************************************#\n")
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('meta_analysis/compare_cell_marker_meta_jack.R')

out_dir <-paste0(disease_dir, '/meta_analysis/low_exp_rm/cell_marker_comparison/')
dir.create(out_dir,showWarnings=F, recursive=T)

(f_marker_genes <- cell_pop_marker_genes_config)
(jack_meta_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))
## meta files
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'meta.*tsv', full.names = T))
threshold_method_ls <- c('top', 'percent', 'q_value')
threshold_ls <- c(100, 0.1, 0.1)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'meta_celltype_fisher_test_results.tsv', rm_na = T)

## jack files
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'jackknife.*tsv', full.names = T))
threshold_method_ls <- c('top', 'percent', 'q_value')
threshold_ls <- c(100, 0.1, 0.2)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'jack_celltype_fisher_test_results.tsv', rm_na = T)

## jack files, just the percent method only
(f_mj_ls <- list.files(jack_meta_folder, pattern = 'jackknife.*tsv', full.names = T))
threshold_method_ls <- c('percent')
threshold_ls <- c(0.05)
mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir,
                 f_result_out = 'jack_celltype_fisher_test_results_top_5_percent.tsv', rm_na = T)


cat("
    #***************************************************************************#
    # PART 4C.5 Plot the meta/jack ranking against each study, calculate spearman correlation (filtered probes)
    #***************************************************************************#\n")
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('meta_analysis/compare_indi_to_meta_jack.R')

## plot out
(plot_out <- paste0(disease_dir, '/meta_analysis/low_exp_rm/ranking_comparisons/'))

## get the meta files and jack files and order them (so that meta and jack are correponding)
(jack_meta_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))
(meta_f_ls <- list.files(jack_meta_folder, pattern = 'meta_genes.tsv', full.names = T))
(jack_f_ls <- list.files(jack_meta_folder, pattern = 'jackknife.tsv', full.names = T))
(meta_f_ls <- meta_f_ls[order(meta_f_ls)])
(jack_f_ls <- jack_f_ls[order(jack_f_ls)])


for (i in 1:length(meta_f_ls)){
    (f_meta <- meta_f_ls[i])
    (f_jack <- jack_f_ls[i])
    mainCompareRanksMetaIndi(f_meta = f_meta, f_jack=f_jack, plot_out=plot_out)
}


cat("
    #***************************************************************************#
    # PART 4D Plot the significant genes (see sec. PART)
    #***************************************************************************#\n")



cat("
    #---------------------------------------------------------------------------#
    # PART 5A: ermineJ for Fisher and jackknife
    #---------------------------------------------------------------------------#\n")
## make the sh script to run from meta and jack files with lowly expressed genes removed
## and manually remove KO genes from meta and jack files
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')
source('ermineJ_preprocess/make_ermineJ_sh.R')
maxsize = 500

(input_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))
bg_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/ermineJ_background/')
erminej_dir <- paste0(disease_dir,'/ermineJ/', Sys.Date(),'_geneset_', maxsize, '/')
xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'

y = mkErminejSH(disease,input_folder, bg_folder, erminej_dir, xml=xml, maxsize = maxsize)

cat("
    #---------------------------------------------------------------------------#
    # PART 5B: analyse ermineJ for Fisher and jackknife, look at the top genes in the top pathways
    #---------------------------------------------------------------------------#\n")
rm(list=setdiff(ls(),'home_dir'))


disease_ls <- c('AD','HD', 'PD')
keyword_ls <- c('meta', 'jackknife')
mixed_model_only =F

## only with low expre rm
regulation_ls <- c('up', 'down')
phase_ls <- c('early', 'late')
threshold = 50   ## number of top pathways to look at
fdr_threshold = 0.05  # or sig_paths to look at if more than threshold


for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    ## where the mixed model result files
    (jack_meta_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))
    (ermineJ_folder <- grep('mixed_model|archive', list.dirs(paste0(disease_dir, 'ermineJ/'), recursive = F), value = T,
                            invert = T))
    
    source('result_explore/top_genes_for_ermineJ.R')
}##loop1








cat("
    #***************************************************************************#
    # -- lowly expressed genes removed
    # PART 6.2 compare the MM results to Fisher results and also prepare the ermineJ files
    # for random intercept
    #***************************************************************************#\n")


rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/compare_mm_meta.R')


source('config_wrappers.R')
model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')
model_keyword <- '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        mm_dir = paste0(disease_dir, 'mixed_model/',model,model_keyword,'/') ## mixed model dir(parent dir)
        mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls=phase_ls, compare_fisher =F)
    }## loop2 end
}##loop1

cat("
    #***************************************************************************#
    # -- lowly expressed genes removed
    # PART 6.3 make ermineJ sh
    #***************************************************************************#\n")
## make the sh script to run from meta and jack files with lowly expressed genes removed as background
rm(list=setdiff(ls(),'home_dir'))

source('ermineJ_preprocess/make_ermineJ_sh.R')
maxsize = 500
disease_ls <- c('AD','HD', 'PD')
model_ls <- c('random_intercept')
model_keyword <- '_include_NA_low_exp_rm'
xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (input_folder <- paste0(disease_dir, 'mixed_model/', model, model_keyword, '/'))
        bg_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/ermineJ_background/')
        erminej_dir <- paste0(disease_dir,'/ermineJ/mixed_model/',model,model_keyword, '/', Sys.Date(),'_geneset_', maxsize, '/')
        
        y = mkErminejSH(disease,input_folder, bg_folder, erminej_dir,xml=xml, maxsize = maxsize)
    }## loop2 end
}##loop1


cat("
    #***************************************************************************#
    # -- lowly expressed genes removed
    # PART 6.4 plot significant genes and do model diagnostics of MM
    # with expression values before quantile normalization
    #***************************************************************************#\n")

# with expression values before quantile normalization
rm(list=setdiff(ls(),'home_dir'))

source('mixed_models/top_genes_mixed_model.R')

disease_ls <- c('AD','HD')
model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')
regulation_ls <- c('up', 'down')
model_keyword <- '_include_NA_low_exp_rm'

threshold = 50   ## number of top genes to be plotted
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        data_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/')
        for (phase in phase_ls){ ## loop 3 for phases
            for(regulation in regulation_ls){
                mainTopMMGenes(data_dir, phase, regulation, threshold=threshold, 
                               model =model)
            }
        } ##loop3 end
    }## loop2 end
}##loop1




cat("
    #***************************************************************************#
    # PART 6C.4.2 JACKKNIFE: make ermineJ sh for both biological process only and all process
    #***************************************************************************#\n")
## make the sh script to run from meta and jack files with lowly expressed genes removed as background
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

source('ermineJ_preprocess/make_ermineJ_sh.R')


model_ls <- c('random_intercept')
model_keyword <- '_include_NA_low_exp_rm'
xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
maxsize = 500
iteration = '500000'

## get both biological and all process enrichment
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories


for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (input_folder <- paste0(disease_dir, 'mixed_model_jackknife/', model, model_keyword, '/'))
        
        for(process in process_ls){
            bg_folder <- paste0(input_folder, '/ermineJ_background',process,'/')
            erminej_dir <- paste0(disease_dir,'/ermineJ/mixed_model_jackknife/',model,model_keyword,process, '/', Sys.Date(),'_geneset_', maxsize, '/')
            y = mkErminejSH(disease,input_folder, bg_folder, erminej_dir,xml=xml, maxsize = maxsize, iteration=iteration)
        }
    }## loop2 end
}##loop1




cat("
    #***************************************************************************#
    # PART 6C.4.3 JACKKNIFE: after MM ermineJ results, look at the top genes in the top pathways
    #***************************************************************************#\n")
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

model_ls <- c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm')  
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')
phase_ls <- c('early', 'late')
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories


top_threshold = 200  ## the top genes to annotate
mixed_model_only=F
threshold = 50   ## number of top pathways to look at
fdr_threshold = 0.05  # or sig_paths to look at if more than threshold
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            (jack_meta_folder <- paste0(disease_dir, 'mixed_model_jackknife/', model,model_keyword,'/'))
            for(process in process_ls){
                ## make sure only greps the folder with date (otherwise will grep the jackknife folder)
                (ermineJ_folder <- max(grep('201', 
                                            list.dirs(paste0(disease_dir, 'ermineJ/mixed_model_jackknife/', model,model_keyword,process, '/'), recursive = F),
                                            value =T)))
                source('result_explore/top_genes_for_ermineJ.R')
            }
        }## loop 3 end
    }## loop2 end
}##loop1





# 
# #*********************#
# ###PART 6C.5 JACKKNIFE  run mixed models: random intercept NA removed (previously defined top genes only)
# #*********************#
# #' re run the mixed model for the top genes, with full report (later for heatmap and more)(if the mixed model run
# #' is for fast mode only when intercept for each study is not recorded)
# #' result in mixed_model_include_NA_low_exp_rm_top_genes
#  
# rm(list=setdiff(ls(),'home_dir'))
# source('mixed_models/mixed_model.R')
# 
# disease_ls = c('AD', 'HD', 'PD')
# phase_ls =c('early', 'late')
# regulation_ls = c('up', 'down')
# 
# model_ls <- c('random_intercept')
# 
# 
# #### first take the top up and down genes from mm jack results
# #' and make a new expression.Rdata that only contains these genes
# #' in mixed_model/random_intercept_include_NA_low_exp_rm_top_genes/early(late)
# threshold = 100  ## for up and down genes each
# source('mixed_models/top_up_down_genes_expression.R')
# 
# 
# 
# 
# #' re run the mixed model for the top genes, with full report (later for heatmap and more)(if the mixed model run
# #' is for fast mode only when intercept for each study is not recorded)
# #' result in mixed_model_include_NA_low_exp_rm_top_genes
# model_keyword = '_include_NA_low_exp_rm_top_genes'
# full_report =T
# rm_genes = ''
# NA_filter='0.3'  ## NA_filter='0.3' to include NA
# REML=F
# estimateCI=T  # to run fast estimateCI=F full_report = F
# disease_stage_only =T
# to_plot=T
# tmp_result_output=F
# 
# 
# for (disease in disease_ls){## loop 1 for disease
#     source('config_wrappers.R')
#     for(model in model_ls){ ## loop2 by model
#         (out_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword , '/')) ## where the mixed model result
#         (data_dir <- paste0(disease_dir,'/mixed_model/', model,model_keyword, '/')) ## where the mixed model input expression
#         
#         for (phase in phase_ls){ ## loop3 by phase
#             (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
#             ### get the expression data and MM results
#             x <- mixedModelAll(phase =phase, out_dir =out_dir, 
#                                exprdata = exprdata, model=model,
#                                full_report = full_report,
#                                to_plot= to_plot, rm_genes = rm_genes, NA_filter = NA_filter,
#                                REML =REML,
#                                estimateCI=estimateCI,
#                                disease_stage_only =disease_stage_only,
#                                tmp_result_output=tmp_result_output)
#             
#             
#             assign(x = paste0(model, "_",phase,"_",disease), value = x)
#         }
#     }## loop2 end
# }##loop1
# 

#---------------------------------------------------------------------------#
# PART 7. For each disease (for mixed model results, not jackknife results)
# plot mixed model results: heatmap for top genes for each study (up, down, up and down), and in 1 plot
# Plot the mixed model ranking against DE gene ranking in each study
# after MM ermineJ sh results, mark the top genes in the top pathways (analysis folder)
# compare between intercept and slope model
# gather plots and info of top mixed model(indluding NA) genes in mixed model, meta and jack
##---------------------------------------------------------------------------#\n")
## see removed steps (this part compare results from mixed models of each disease, not the jackknife results of each disease)

