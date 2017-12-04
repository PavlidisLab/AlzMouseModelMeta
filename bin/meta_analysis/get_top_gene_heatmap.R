#' 2016-01-01
#' update 2016-05-30
#' updated 2016-08-02 modified to accomondate mixed model results
#' ################################################################################################
#' Plot heatmap for experiments of provided genelist after limma DE
#' 
#' load the saved R objects with expression data and annotation
#' plot the heatmap of a given list of genes (each gene is represented by the best logFC probes, including NA)
#' required function topGeneHeatmap() from "helper_functions.R"
#' required filterLowExprGenes() from 'ermineJ_preprocess/filter_out_low_exp_for_enrichment.R'
#' 
#' input: a list of gene symbols, input R object dir and output plot dir
#' output: heatmap plots with gene expression scaled and not scaled for comparison
#' ----------------------------------------------------- #
#' USAGE: to run by range of rows, specified rows, or datasets
# rm(list=setdiff(ls(),'home_dir'))
#  
# 
# ## a dir contains all the R objects saved after limma DE
# r_object_dir <-'/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/all_prioritized_limma/'
# ## where the jack and meta files
# jack_meta_folder <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-04-05/'
# # plot out dir (a sub folder of 'gene_heatmaps' will be created)
# plt_out <- '/home/bzhuang/AD_mouse_model_project/meta_analysis'
# low_exp_rm <- T  # whether to remove lowly expressed probes (not the input genes, so that the expression of 
#                   the gene will be NA?? need to think about this)
# threshold <- 50  #top genes threshold
#  cluster_rows <- T  # for heatmap whether to cluster by rows/probes
# 
# source('meta_analysis/get_top_gene_heatmap.R')
#' ----------------------------------------------------- # 
#' 


#---------------------------------------------------------------------------#
# PART 1: SET SYS DEFAULTS
#---------------------------------------------------------------------------#
 
source("helper_functions.R")
source('ermineJ_preprocess/filter_out_low_exp_for_enrichment.R')  # use the filterLowExprGenes function to filter out lowly exp probes
width = 1600
height = 1000


#---------------------------------------------------------------------------#
# PART 2: SET INPUTS
#---------------------------------------------------------------------------#
(f_ls_all <- list.files(path = r_object_dir, recursive = T, pattern= '.Rdata', full.names=T))
## a df with 'Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype', 'File'(limma objectdir), file_labels(for plt)
info_df <- getRObjectLabels(f_ls_all, df_info)

#---------------------------------------------------------------------------#
# PART 3: LOOPs
#---------------------------------------------------------------------------#

for(keyword in keyword_ls){
  print(paste0('start ', keyword))
  for(regulation in regulation_ls){
    print(paste0('regulation: ', regulation))
    
    ## get the meta or jack files
    (f_gene_list_all <- list.files(jack_meta_folder, full.names=T,
                                   pattern=paste0(regulation, '_regulation_',keyword, '.tsv')))
    for(j in 1:length(f_gene_list_all)){
      
      (f_gene_list <-f_gene_list_all[j])
        
        
        #***************************
        #' loop for each disease phase and correcponding studies
        ## f_gene_list: a meta or jack file
        ## threshold: top # of genes in the f_gene_list
        ## info_df: a df with File(the limma object files), Phase, Extra_phase
        #***************************
        (disease_phase <- grep('early|late|intermediate', unlist(strsplit(f_gene_list, split = "/|_")), value = T))
        
        print(f_gene_list)
        
        ## get the genes from jack or meta file
        df <- read.delim(f_gene_list,comment.char = "#")
        
        ## reorder by the score column
        (score_col <- intersect(c('adj_combined_max_p','Fisher', 'up_pval', 'down_pval'), colnames(df)))
        df <- df[order(df[, score_col]), ]
        
        gene_list <- as.character(df$geneSymbol[1:threshold])
        
        
        ## get the studies corresponding to the time
        (selected_study <- filterContain(info_df, "Phase", disease_phase))  ## TODO: need to update add genotype
        if(nrow(selected_study) == 0){  # for more phases
            selected_study <- filterContain(info_df, "Extra_phase", disease_phase)
        }
        (f_ls <- selected_study$File) # contain a list of object dirs of the specified phase

      
      ###### LOOP for each study and plot heatmap of input gene list
      for (i in 1:length(f_ls)){
        print(f_ls[i])
        
        load(f_ls[i])
        print(selected_study$file_labels[i])
        ## get timepoint and plot prefix
        (timepoint <- as.character(selected_study$Timepoint[i] ))
        (dataset_p <- as.character(selected_study$file_labels[i]))
        (genotype <- as.character(selected_study$Genotype[i]))
        

        ## only use genotype in the heatmap
        array_design <- array_design[, c("Bioassay","ExternalID","Genotype", "Sample", "Dataset")]%>% droplevels()
        
        ##************************##
        ## subset the array by genotype
        ##************************##
        ## only plot the specified Genotype
        #array_design <- filterContain(array_design, 'Genotype', c('WT', genotype))
        array_design <- array_design[which(array_design$Genotype %in% c('WT', genotype)),]%>% droplevels()
        array_dat <- array_dat[, as.character(array_design$Sample)]
        
        
        ##************************##
        ## if lowly expressed probes is removed from heatmap
        ##************************##
        if(low_exp_rm){
            # do not filter agilent (ratio based) for now 2016-05-17
          array_dat <- filterLowExprGenes(file_dir_ls = f_ls[i], affy_threshold = 6, agilent_threshold = -10, 
                                         filter_method ="max", return_array = T)
        }
        
        ##************************##
        ## plot output dir
        ##************************##
        
         plotdir_1 <- paste0(plt_out, '/gene_heatmaps/', Sys.Date(), '/not_clustered/')
        
        (plotdir <- paste0(plotdir_1, '/all_probes/', disease_phase,"_", keyword, "_",
                           regulation, "/"))
        if(low_exp_rm){
          (plotdir <- paste0(plotdir_1, '/low_exp_rm/', disease_phase,"_", keyword, "_",
                             regulation, "/"))
        }
        print(paste0('plotdir is ', plotdir))
        
        dir.create(plotdir, showWarnings=F, recursive=T)
        
        ################################
        ## plot heatmap, scaled by genes
        ## pick the best probe per gene
        ## keep genes with no probes (NA)
        ################################
        (plt_title <- paste0(dataset_p, "_top_", threshold, "_", regulation, "_regulated_", keyword))
        (plot_pre <- plotdir)
        # topGeneHeatmap from helper_function.R
        # plot name is paste0(plot_pre, plt_title, prefix, '.png')

        topGeneHeatmap(array_dat, gene_list, array_design, annotation, plot_pre, 
                       height, width, 
                       prefix = "_scaled",
                       plt_title = plt_title, 
                       scale_row =T,
                       write_df =T,
                       best_probe =T,
                       array_eb_fit = array_eb_fit,
                       legend_breaks =  -2:2, 
                       legend_labels = -2:2)
        
        #### for figure quality
        ## for pub not to show sample names, plots are the same dimentions regardless of sample size
        # width based on the sample size, height bease on the number of probes/genes
        p_w <- 25 * ncol(array_dat)+ 153
        p_h <- 25 * length(gene_list)
        topGeneHeatmap(array_dat, gene_list, array_design, annotation, plot_pre, 
                       height = p_h, 
                       width = p_w, 
                       prefix = "",
                       plt_title = paste0(dataset_p),
                       scale_row =T,
                       show_colnames =F,
                       write_df =F,
                       best_probe =T,
                       auto_w_d= F,
                       array_eb_fit = array_eb_fit,
                       display_full_name = F,
                       legend_breaks =  -2:2, 
                       legend_labels = -2:2,
                       size_r=24)
        
        
        
        ################################
        ## plot heatmap, not scaled , all probes
        ################################
        topGeneHeatmap(array_dat, gene_list, array_design, annotation, plot_pre, plt_title = plt_title, 
                           height, width, prefix = "", scale_row =F,best_probe =F)
        
        ################################
        ## if cluster genes, all probes
        ################################
        if(cluster_rows){
          ## redefine plot dir
            plotdir_1 <- paste0(plt_out, '/gene_heatmaps/clustered/')
            (plotdir <- paste0(plotdir_1, '/all_probes/', disease_phase,"_", keyword, "_",
                               regulation, "/"))
            
            
            if(low_exp_rm){
              (plotdir <- paste0(plotdir_1, '/low_exp_rm/', disease_phase,"_", keyword, "_",
                                 regulation, "/"))
            }
            print(paste0('plotdir is ', plotdir))
            dir.create(plotdir, showWarnings=F, recursive=T)
            
          
            ## start plots
          (plt_title <- paste0(dataset_p, "_top_", threshold, "_", regulation, "_regulated_", keyword))
          (plot_pre <- plotdir)
          ## plot heatmap, scaled by genes, clustered
          topGeneHeatmap(array_dat, gene_list, array_design, annotation, plot_pre, 
                             height, width, 
                             prefix = "_scaled",
                             plt_title = plt_title, 
                             scale_row =T,
                             cluster_rows = T,
                             best_probe =F,
                             legend_breaks =  -2:2, 
                             legend_labels = -2:2) 
          
          ## plot heatmap, not scaled, clustered
          topGeneHeatmap(array_dat, gene_list, array_design, annotation, plot_pre, plt_title = plt_title, 
                             height, width, prefix = "", scale_row =F,  cluster_rows = T, best_probe =F)
        }
      }
    }
  }
}


