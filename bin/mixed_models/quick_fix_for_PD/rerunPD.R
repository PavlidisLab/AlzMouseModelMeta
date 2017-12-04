#2016-08-20, fix Pink1 expression and rerun analysis
home_dir <- '/home/bzhuang/ND_project_combined/'
#home_dir <- "C:/Users/User/Documents/lab_results/ND_project_combined/"
rm(list=setdiff(ls(),'home_dir'))
disease_ls <- c('PD')
phase_ls <- c('early', 'late')
model_keyword =  '_include_NA_low_exp_rm' ## '_include_NA' or '_include_NA_low_exp_rm' to specify which model folder
model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')


keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')
model_keyword_ls <- c('_include_NA_low_exp_rm')  

cat("
    #***************************************************************************#
    # -- lowly expressed genes removed
    # PART 6B.2 compare the MM results to meta analysis results and also prepare the ermineJ files
    # for random slope and randome intercept
    #***************************************************************************#\n")

setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))

source('mixed_models/compare_mm_meta.R')


for (disease in disease_ls){## loop 1 for disease
    source('config/config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        mm_dir = paste0(disease_dir, 'mixed_model/',model,model_keyword,'/') ## mixed model dir(parent dir)
        mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls=phase_ls)
    }## loop2 end
}##loop1

# cat("
#     #***************************************************************************#
#     # -- lowly expressed genes removed
#     # PART 6B.3 make ermineJ sh
#     #***************************************************************************#\n")
# ## make the sh script to run from meta and jack files with lowly expressed genes removed as background
# 
# setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
# source('ermineJ_preprocess/make_ermineJ_sh.R')
# 
# 
# for (disease in disease_ls){## loop 1 for disease
#     print(disease)
#     source('config/config_wrappers.R')
#     for (model in model_ls){## loop2 for models
#         (input_folder <- paste0(disease_dir, 'mixed_model/', model, model_keyword, '/'))
#         bg_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/ermineJ_background/')
#         erminej_dir <- paste0(disease_dir,'/ermineJ/mixed_model/',model,model_keyword, '/', Sys.Date(),'/')
#         
#         mkErminejSH(disease,input_folder, bg_folder, erminej_dir)
#     }## loop2 end
# }##loop1


cat("
    #***************************************************************************#
    # -- lowly expressed genes removed
    # PART 6B.4 plot significant genes and do model diagnostics of MM
    #***************************************************************************#\n")

# with expression values before quantile normalization

setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source('mixed_models/top_genes_mixed_model.R')


phase_ls <- c('early', 'late')
regulation_ls <- c('up', 'down')
model_keyword <- '_include_NA_low_exp_rm'

threshold = 50   ## number of top genes to be plotted
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config/config_wrappers.R')
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

#---------------------------------------------------------------------------#
# PART 7. For each disease
# plot mixed model results: heatmap for top genes for each study (up, down, up and down), and in 1 plot
# Plot the mixed model ranking against DE gene ranking in each study
# after MM ermineJ sh results, mark the top genes in the top pathways (analysis folder)
# compare between intercept and slope model
# gather plots and info of top mixed model(indluding NA) genes in mixed model, meta and jack
##---------------------------------------------------------------------------#\n")


cat("
    #***************************************************************************#
    # PART 7.1 plot heatmap of top significant genes from results of MM (1 direction)
    # low exp rm or not
    #***************************************************************************#\n")
#
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source('helper_functions.R')

#disease_ls <- c('AD','HD', 'PD')
#model_ls <- c('random_intercept',  'random_slope')

threshold = 50   ## number of top genes to be plotted
low_exp_rm <- F  # whether to remove lowly expressed probes # hear the low probes are already removed from the meta genes
cluster_rows <- F  # for heatmap whether to cluster by rows/probes

for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config/config_wrappers.R')
    r_object_dir <-paste0(limma_dir, 'prioritized_limma/')
    df_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv') # for early, intermediate and late files
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/')
            # plot out dir (a sub folder of 'gene_heatmaps' will be created)
            plt_out <- paste0(jack_meta_folder)
            source('meta_analysis/get_top_gene_heatmap.R')
        }## loop 3 end
    }## loop2 end
}##loop1



cat("
    #***************************************************************************#
    # PART 7.2 plot heatmap of top significant genes from results of MM (up and down in 1 plot)
    # for thesis figure
    #***************************************************************************#\n")
#
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source('helper_functions.R')

threshold = 20   ## number of top genes to be plotted
low_exp_rm <- F  # whether to remove lowly expressed probes # hear the low probes are already removed from the meta genes
cluster_rows <- F  # for heatmap whether to cluster by rows/probes

for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config/config_wrappers.R')
    r_object_dir <-paste0(limma_dir, 'prioritized_limma/')
    df_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv') # for early, intermediate and late files
    #df_info <- '../../doc/PD_mouse_dataset_doc/dataset_info_mixed_model.tsv'
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/')
            # plot out dir (a sub folder of 'gene_heatmaps' will be created)
            plt_out <- paste0(jack_meta_folder)
            source('meta_analysis/get_top_gene_heatmap_up_down_figure.R')
        }## loop 3 end
    }## loop2 end
}##loop1

cat("
    #***************************************************************************#
    # PART 7.3 plot heatmap of top significant genes from results of MM (up and down in 1 plot)
    # for thesis figure-all datasets in 1 plot
    #***************************************************************************#\n")
#
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source('helper_functions.R')

disease_ls <- 'PD'
model_ls


threshold = 20   ## number of top genes to be plotted
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config/config_wrappers.R')
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/')
            mm_dir <- jack_meta_folder
            # plot out dir (a sub folder of 'gene_heatmaps' will be created)
            plt_out <- paste0(jack_meta_folder)
            source('mixed_models/top_gene_models_mm_one_panel.R')
        }## loop 3 end
    }## loop2 end
}##loop1




cat("
    #***************************************************************************#
    # PART 7.4 Plot the mixed model ranking against each study, calculate spearman correlation (filtered probes)
    #***************************************************************************#\n")
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))

source('meta_analysis/compare_indi_to_meta_jack.R')
for(disease in disease_ls){
    source('config/config_wrappers.R')
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            
            ## where the mixed model result files and meta files
            (jack_meta_folder <- paste0(disease_dir, 'meta_analysis/low_exp_rm/'))
            (meta_f_ls <- list.files(jack_meta_folder, pattern = 'meta_genes.tsv', full.names = T))
            (mixed_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/'))
            (mix_ls <- list.files(mixed_folder, pattern = 'mixed_model.tsv', full.names = T))
            (meta_f_ls <- meta_f_ls[order(meta_f_ls)])
            (mix_ls <- mix_ls[order(mix_ls)])
            
            (plot_out <- paste0(mixed_folder, '/ranking_comparisons/'))
            ## loop for each file
            for (i in 1:length(meta_f_ls)){
                (f_meta <- meta_f_ls[i])
                (f_jack <- mix_ls[i])
                mainCompareRanksMetaIndi(f_meta = f_meta, f_jack=f_jack, plot_out=plot_out,plot_meta =F)
            }
        }## loop 3 end
    }## loop2 end
}

cat("
    #***************************************************************************#
    # PART 7.5: after MM ermineJ results, look at the top genes in the top pathways
    #***************************************************************************#\n")

# rm(list=setdiff(ls(),'home_dir'))
# setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
# 
# 
# disease_ls <- c('AD','HD', 'PD')
# model_ls <- c('random_intercept',  'random_slope')
# model_keyword_ls <- c('_include_NA', '_include_NA_low_exp_rm')  
# keyword_ls <- c('mixed_model')
# regulation_ls <- c('up', 'down')
# phase_ls <- c('early', 'late')


disease_ls <- c('AD','HD')
model_ls <- c('random_intercept',  'random_slope')
model_keyword_ls <- c('_include_NA', '_include_NA_low_exp_rm')  
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')
phase_ls <- c('early', 'late')
# 

top_threshold = 200  ## the top genes to annotate
mixed_model_only=F
threshold = 50   ## number of top pathways to look at
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config/config_wrappers.R')
    
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            (jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/'))
            (ermineJ_folder <- max(list.dirs(paste0(disease_dir, 'ermineJ/mixed_model/', model,model_keyword,'/'), recursive = F)))
            
            source('result_explore/top_genes_for_ermineJ.R')
            
        }## loop 3 end
    }## loop2 end
}##loop1
