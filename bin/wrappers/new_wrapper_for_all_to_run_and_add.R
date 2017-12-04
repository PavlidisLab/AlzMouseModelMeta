

############# running here ############



cat("
    #**************************#
    # PART 9A.1.2 plot the top genes expression for each ND  ## may need to change (2016-9-023)
    # result in '/ND_results/top_genes_expression/'
    #**************************#\n")
rm(list=setdiff(ls(),'home_dir'))

source('helper_functions.R')
source('config_wrappers.R')

model_ls <- c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm')  
phase_ls <- c('early', 'late')
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')

(plotdir <- paste0(home_dir, '/ND_results/top_genes_expression/', Sys.Date(), '/'))


threshold = 10   ## number of top genes to be plotted
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            (jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/'))
            mm_dir <- jack_meta_folder
            # plot out dir (a sub folder of 'gene_heatmaps' will be created)
            source('mixed_models/plot_up_down_genes_thesis_figure.R')
        }## loop 3 end
    }## loop2 end
}##loop1



cat("
    #**************************#
    # PART 9A.2.1 make ermineJ background (mm and mm jack use the same background) 
    #  used mixed model jack as the background to combine
    # both mixed model and mixed model jackknife use the same background
    #**************************#\n")
## low expression removed , make from jackknife background
rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')

phase_ls =c('early','late')
regulation_ls = c('up', 'down')

meta_jack_ls = 'mixed_model'

out_dir <- paste0(home_dir, '/ND_results/ermineJ/background/')
source('ermineJ_preprocess/combine_bg_for_all_ND.R')


cat("
    #**************************#
    # PART 9A.2.2 make ermineJ sh (for mm only)
    #**************************#\n")
rm(list=setdiff(ls(),'home_dir'))
disease_erminej= 'ALL'
xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
scorecol=2  ## must use columne 2 (rank_avg_percent for rank column)
source('ermineJ_preprocess/make_ermineJ_sh.R')

maxsize=500

(input_folder <- paste0(home_dir, '/ND_results/Disease_comparison/'))
(bg_folder <- paste0(home_dir, '/ND_results/ermineJ/background/'))
(erminej_dir <- paste0(home_dir, '/ND_results/ermineJ/results/', Sys.Date(),'_geneset_', maxsize, '/'))
y = mkErminejSH(disease=disease_erminej,input_folder, bg_folder, erminej_dir, scorecol=scorecol,xml=xml, maxsize=maxsize)


cat("
    #**************************#
    # after ermineJ sh (mm only, and Fisher)
    # PART 9A.2.3 after sh ermineJ look at the top genes in the top pathways for combined disease results
    #**************************#\n")

## for each gene, all plots and mixed model, jack, meta values
rm(list=setdiff(ls(),'home_dir'))


regulation_ls <- c('up', 'down')
phase_ls <- c('early', 'late')
keyword_ls <- c('mixed_model')

threshold = 50   ## number of top pathways to look at
fdr_threshold = 0.05  # or sig_paths to look at if more than threshold
mixed_model_only=T
model_ls <- c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm')  ## only with low exp removed

## for mm
for (model in model_ls){## loop1 for models
    for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
        ## where the mixed model result files
        
        (jack_meta_folder <- paste0(home_dir, '/ND_results/Disease_comparison/'))
        (ermineJ_folder <- max(grep('201', list.dirs(paste0(home_dir, '/ND_results/ermineJ/results/'), recursive = F), value = T)))
        source('result_explore/top_genes_for_ermineJ.R')
    }## loop 2 end
}## loop1 end


#######
## for meta and jack (not mixed model)
#######
keyword_ls <- c('meta', 'jackknife')
mixed_model_only=F
(jack_meta_folder <- paste0(home_dir, '/ND_results/Disease_comparison/'))
(ermineJ_folder <- max(list.dirs(paste0(home_dir, '/ND_results/ermineJ/results/'), recursive = F)))
source('result_explore/top_genes_for_ermineJ.R')



cat("
    #---------------------------------------------------------------------------#
    # PART 9B : compare results from all diseases(for mixed model jackknife results)
    #---------------------------------------------------------------------------#\n")

cat("
    #**************************#
    # PART 9B.1.1 combine all 3 disease (mm jack only): 
    #        intercept and slope including NA)
    #        results are in Disease_comparison_jackknife/
    #**************************#\n")

rm(list=setdiff(ls(),'home_dir'))


disease_ls= c('AD','HD')
model_ls = c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm')  ## only with low exp removed

## MM jackknife
correlation =T
threshold =500
threshold_phase =500
out_dir <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/')
source('result_explore/compare_disease_jackknife.R')

# 
# cat("
#     #**************************#
#     # PART 9B.2.1 make ermineJ background (mm and mm jack use the same background) 
#     #  used mixed model jack as the background to combine
#     #**************************#\n")
# ## low expression removed , make from jackknife background
# # both mixed model and mixed model jackknife use the same background, see part 9A
# ## low expression removed , make from jackknife background
# rm(list=setdiff(ls(),'home_dir'))
# 
# disease_ls <- c('AD', 'HD')
# phase_ls =c('early','late')
# regulation_ls = c('up', 'down')
# 
# meta_jack_ls = 'mixed_model'
# 
# out_dir <- paste0(home_dir, '/ND_results/ermineJ/background/')
# source('ermineJ_preprocess/combine_bg_for_all_ND.R')
# 
# 
# 
# cat("
#     #**************************#
#     # PART 9B.2.2 make ermineJ sh (for mm jack)
#     #**************************#\n")
# rm(list=setdiff(ls(),'home_dir'))
# disease_erminej= 'ALL_jack'
# xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
# scorecol=2  ## must use columne 2 (rank_avg_percent for rank column)
# maxsize = 500
# source('ermineJ_preprocess/make_ermineJ_sh.R')
# (input_folder <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/'))
# (bg_folder <- paste0(home_dir, '/ND_results/ermineJ/background/'))
# (erminej_dir <- paste0(home_dir, '/ND_results/ermineJ/results_jackknife/', Sys.Date(),'_geneset_', maxsize, '/'))
# y = mkErminejSH(disease=disease_erminej,input_folder, bg_folder, erminej_dir, scorecol=scorecol,xml=xml, maxsize=maxsize)
# 
# 
# 
# cat("
#     #**************************#
#     # after ermineJ sh (mm and mm jack)
#     # PART 9B.2.3 after sh ermineJ look at the top genes in the top pathways for combined disease results
#     #**************************#\n")
# 
# ## for each gene, all plots and mixed model, jack, meta values
# rm(list=setdiff(ls(),'home_dir'))
#  
# 
# regulation_ls <- c('up', 'down')
# phase_ls <- c('early', 'late')
# keyword_ls <- c('mixed_model')
# 
# threshold = 50   ## number of top pathways to look at
# fdr_threshold = 0.05  # or sig_paths to look at if more than threshold
# mixed_model_only=T
# model_ls <- c('random_intercept')
# model_keyword_ls <- c('_include_NA_low_exp_rm')  ## only with low exp removed
# 
# 
# ## for mm jack
# for (model in model_ls){## loop1 for models
#     for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
#         ## where the mixed model result files
#         
#         (jack_meta_folder <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/'))
#         (ermineJ_folder <- max(grep('201', list.dirs(paste0(home_dir, '/ND_results/ermineJ/results_jackknife/'), recursive = F), value = T)))
#         source('result_explore/top_genes_for_ermineJ.R')
#     }## loop 2 end
# }## loop1 end
# 
# 


cat("
    #---------------------------------------------------------------------------#
    # PART 10A : compare results from any combo of the 2 diseases () from mixed model (not jackknife)
    # results are in Disease_comparison/AD_PD, AD_HD, HD_PD
    # only for mixed model lowly expressed genes removed
    # define disease_comb_ls, not disease_ls
    #---------------------------------------------------------------------------#\n")
# see removed steps



cat("
    #---------------------------------------------------------------------------#
    # PART 10B : JACKKNIFE MM: compare results from any combo of the 2 diseases ()
    # results are in Disease_comparison_jackknife/AD_PD, AD_HD, HD_PD
    # only for mixed model lowly expressed genes removed
    # define disease_comb_ls, not disease_ls
    #---------------------------------------------------------------------------#\n")

cat("
    #**************************#
    # PART 10B.1 combine all 3 disease for MM jack
    # prepare ermineJ gene rank input
    #**************************#\n")
disease_comb_ls = c('AD', 'HD')
model_ls = c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm','_include_NA_low_exp_rm_adj_cell_pop')  ## only with low exp removed
correlation =F
threshold =500
threshold_phase =500

disease_combs <- as.data.frame(combn(disease_comb_ls, m=2))

## MM Jack
for(i in 1: ncol(disease_combs)){
    disease_ls <- (as.character(disease_combs[, i]))
    print(disease_ls)
    (out_dir <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/', paste0(sort(disease_ls), collapse='_'), '/'))
    source('result_explore/compare_disease_jackknife.R')
}

# cat("
#     #**************************#
#     # PART 10B.2.1 make ermineJ background(mm and mm jack share the same background)
#     #**************************#\n")
# ## low expression removed , make from jackknife background
# ## make ermineJ background
#     # (mm and mm jack share the same background)
# 
# ## low expression removed , make from MM jackknife background
# rm(list=setdiff(ls(),'home_dir'))
# disease_comb_ls = c('AD', 'HD')
# phase_ls =c('early','late')
# regulation_ls = c('up', 'down')
# meta_jack_ls = 'mixed_model'
# process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories
# 
# disease_combs <- as.data.frame(combn(disease_comb_ls, m=2))
# for(i in 1: ncol(disease_combs)){
#     disease_ls <- (as.character(disease_combs[, i]))
#     print(disease_ls)
#     for(process_keyword in process_ls){
#         (out_dir <- paste0(home_dir, '/ND_results/ermineJ/background',process_keyword,'/', paste0(sort(disease_ls), collapse='_'), '/'))
#         print(out_dir)
#         source('ermineJ_preprocess/combine_bg_for_all_ND.R')
#     }
# 
# }
# 
# 
# cat("
#     #**************************#
#     # PART 10B.2.2 make ermineJ sh (mm J)
#     #**************************#\n")
# rm(list=setdiff(ls(),'home_dir'))
# 
# disease_comb_ls = c('AD', 'HD','PD')
# 
# xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
# maxsize = 500
# source('ermineJ_preprocess/make_ermineJ_sh.R')
# disease_combs <- as.data.frame(combn(disease_comb_ls, m=2))
# for(i in 1: ncol(disease_combs)){
#     disease_ls <- (as.character(disease_combs[, i]))
#     print(disease_ls)
#     (disease_erminej= paste0(sort(disease_ls), collapse='_'))
#     scorecol=2  ## must use columne 2 (rank_avg_percent for rank column)
#     
#     ## for mm jack
#     (input_folder <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/', paste0(sort(disease_ls), collapse='_'), '/'))
#     (bg_folder <- paste0(home_dir, '/ND_results/ermineJ/background/', paste0(sort(disease_ls), collapse='_'), '/'))
#     (erminej_dir <- paste0(home_dir, '/ND_results/ermineJ/results_jackknife/', paste0(sort(disease_ls), collapse='_'), '/', Sys.Date(),'_geneset_', maxsize, '/'))
#     y = mkErminejSH(disease=disease_erminej,input_folder, bg_folder, erminej_dir, scorecol=scorecol,xml=xml, maxsize = maxsize)
# }
# 
# 
# 
# 
# for(i in 1: ncol(disease_combs)){
#     disease_ls <- (as.character(disease_combs[, i]))
#     print(disease_ls)
#     (disease_erminej= paste0(sort(disease_ls), collapse='_'))
#     scorecol=2  ## must use columne 2 (rank_avg_percent for rank column)
#     
#     ## for mm jack
#     (input_folder <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife_adj_cell_pop/', paste0(sort(disease_ls), collapse='_'), '/'))
#     (bg_folder <- paste0(home_dir, '/ND_results/ermineJ/background/', paste0(sort(disease_ls), collapse='_'), '/'))
#     (erminej_dir <- paste0(home_dir, '/ND_results/ermineJ/results_jackknife_adj_cell_pop/', paste0(sort(disease_ls), collapse='_'), '/', Sys.Date(),'_geneset_', maxsize, '/'))
#     y = mkErminejSH(disease=disease_erminej,input_folder, bg_folder, erminej_dir, scorecol=scorecol,xml=xml, maxsize = maxsize)
# }
# 
# cat("
#     #**************************#
#     # after ermineJ sh: MM jack
#     # PART 10B.2.3 look at the top genes in the top pathways for combined disease results
#     #**************************#\n")
# 
# 
# ## for each gene, all plots and mixed model, jack, meta values
# rm(list=setdiff(ls(),'home_dir'))
#  
# 
# disease_comb_ls = c('AD', 'HD','PD')
# 
# regulation_ls <- c('up', 'down')
# phase_ls <- c('early', 'late')
# keyword_ls <- c('mixed_model')
# 
# threshold = 50   ## number of top pathways to look at
# fdr_threshold = 0.05  # or sig_paths to look at if more than threshold
# #######
# ## only for mixed model jack
# #######
# mixed_model_only=T
# model_ls <- c('random_intercept')
# model_keyword_ls <- c('_include_NA_low_exp_rm')  ## only with low exp removed
# 
# disease_combs <- as.data.frame(combn(disease_comb_ls, m=2))
# 
# ## for jackknife
# for(i in 1: ncol(disease_combs)){ ## loop 0 for each combo of 2 diseases
#     disease_ls <- (as.character(disease_combs[, i]))
#     print(disease_ls)
#     for (model in model_ls){## loop1 for models
#         for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
#             ## where the mixed model result files
#             (jack_meta_folder <- paste0(home_dir, '/ND_results/Disease_comparison_jackknife/', paste0(sort(disease_ls), collapse='_'), '/'))
#             (ermineJ_folder <- max(list.dirs(paste0(home_dir, '/ND_results/ermineJ/results_jackknife/', paste0(sort(disease_ls), collapse='_'), '/'), recursive = F)))
#             source('result_explore/top_genes_for_ermineJ.R')
#         }## loop 2 end
#     }## loop1 end
# }
# 
# 
# 
