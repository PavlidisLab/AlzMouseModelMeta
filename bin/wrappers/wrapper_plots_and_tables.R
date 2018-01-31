
setwd(file.path(here::here(),'bin'))

cat("
    #---------------------------------------------------------------------------#
    # PART 11 thesis tables, figures, rdata for easy comparison
    #---------------------------------------------------------------------------#\n")


cat("
    #**************************#
    # PART 11.0 Rdata for all jacknife results for all disease
    #**************************#\n")
#' make all mixed model jackknife results into a r data, and save to 
#' 'ND_results/mm_results/mm_jack_each_disease.Rdata'

rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

regulation_ls = c('up', 'down')
phase_ls =c('early', 'late')

## for results of mixed models without MGP
variable_prefix =''
input_dir = 'mixed_model_jackknife/random_intercept_include_NA_low_exp_rm/'  # where to grab the mixed_model_results.tsv
output_r ='mm_jack_each_disease.Rdata'  ## output r data name
source('summary_tables/summary_rdata/mm_jack_summary_each_disease.R')

## for results of mixed models with MGP
variable_prefix = '_adj_cell'
input_dir = 'mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/'  # where to grab the mixed_model_results.tsv
output_r ='mm_jack_each_disease_adj_cell.Rdata'  ## output r data name
source('summary_tables/summary_rdata/mm_jack_summary_each_disease.R')




cat("
    #**************************#
    # PART 11.3 plot Heatmaps
    #**************************#\n")
cat("
    #-------------#
    # PART 11.3.1 make the heatmap for top jackknife genes (top genes of the phase)- with cell population corrected
    # and top mixed model genes (not used in thesis)
    # expression are corrected to the corresponding gene list (corrected for study or corrected for celltypes)
    # also plot the top genes of the opposite phase expression
    #-------------#\n")
## wrapper
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('helper_functions.R')
source('mixed_models/plot_up_down_genes_thesis_figure_corrected_value.R')


phase_ls <- c('early','late')
# result_rank_ls <- c('mixed_model', 'mixed_model_jackknife')  ## which rank list to get genes from
result_rank_ls <- c('mixed_model_jackknife')  ## which rank list to get genes from
plot_dir <- paste0(home_dir,'/ND_results/top_gene_heatmaps/')

#legend_ls <- list(c('Genotype', 'Study'), c('Genotype', 'Original_genotype', 'Study'),c('Genotype', 'Model_types'))
legend_ls= list(c('Genotype', 'Model_types'))

model_ls <- c('random_intercept')
# mm_rdata_keyword_ls <- c('_include_NA_low_exp_rm','_include_NA_low_exp_rm_adj_cell_pop')  ## choose which top genes for plotting 
mm_rdata_keyword_ls <- c('_include_NA_low_exp_rm_adj_cell_pop')


regulation_ls <- c('up', 'down')
rdata_keyword = 'mixed_model_results_exp_corrected.Rdata'  # choose which rdata to get from: corrected or not

rm_gene_list = c('Thy1')  ## remove these genes from the gene rank list ## if no gene rm: rm_gene_list = NULL
#rm_gene_list = c('Pink1')  ## remove these genes from the gene rank list ## if no gene rm: rm_gene_list = NULL
#rm_gene_list = NULL

threshold = 20   ## number of top genes to be plotted
top_genes =T ## to plot the top genes, otherwise define gene_list to plot the specified genes
opposite_phase_ls = c(F, T)
rdata_keyword = 'mixed_model_results_exp_corrected.Rdata'

for (disease in disease_ls){## loop 1 for disease
    print(disease)
    geno_f=paste0('../configs/', disease, '_mouse_dataset_doc/dataset_info_genotypes.tsv')
    
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        for(result_rank in result_rank_ls){ ## loop3, jackknife or not
            for(opposite_phase in opposite_phase_ls){ ## loop4 produce which set of heatmaps
                for(mm_rdata_keyword in mm_rdata_keyword_ls){ # loop5: correct for celltype or not
                    ## where the mixed model result files
                    mm_rdata_dir <- paste0(disease_dir, 'mixed_model/', model,mm_rdata_keyword,'/') ## where the expression results
                    mm_dir <- paste0(disease_dir,result_rank,'/', model,mm_rdata_keyword,'/') ## where the ranked results
                    if(opposite_phase){
                        (plt_out <- paste0(plot_dir,disease,'/', result_rank,'/opposite/',model,mm_rdata_keyword,'/', Sys.Date(), '/'))
                    }else{
                        (plt_out <- paste0(plot_dir,disease,'/', result_rank,'/',model,mm_rdata_keyword,'/', Sys.Date(), '/'))
                    }
                    
                    plotThesisHeatGenes(mm_dir,mm_rdata_dir,
                                        plt_out,
                                        phase_ls, 
                                        threshold,
                                        rdata_keyword = rdata_keyword,
                                        rm_gene_list = rm_gene_list,
                                        top_genes =top_genes, 
                                        opposite_phase =opposite_phase, 
                                        gene_list =gene_list,
                                        row_width =25,
                                        legend_ls =legend_ls,
                                        plot_indi = T,
                                        geno_f=geno_f)
                }# loop5: correct for celltype or not
            }## loop4 end
        }## loop 3 end
    }## loop2 end
}##loop1


cat("
    #-------------#
    # PART 11.3.2 make the heatmap for cell type marker genes ## just to explore, not included in thesis
    # input expression is corrected for study only
    #-------------#\n")

#' with model_keyword_ls <- c('_include_NA_low_exp_rm')  ## this correct all the genes if intercept is provided

rm(list=setdiff(ls(),'home_dir'))

source('helper_functions.R')
source('mixed_models/plot_up_down_genes_thesis_figure_corrected_value.R')
source('config_wrappers.R')

legend_ls <- list(c('Genotype', 'Study'), c('Genotype', 'Original_genotype', 'Study'))


model <- c('random_intercept')
model_keyword <- c('_include_NA_low_exp_rm')  
phase_ls <- c('early', 'late')
regulation_ls <- c('up', 'down')
result_rank <- c('mixed_model_jackknife')  ## which rank list to get genes from
opposite_phase = F
top_genes = F

## all cell types (hippocampus and striatum)
cell_type_ls = c('Microglia', 'Astrocyte', 'DentateGranule', 'GabaSTTReln', 'Oligo', 'Pyramidal_Thy1',
                 'ForebrainCholin', 'Spiny')
mm_rdata_keyword_ls <- c('_include_NA_low_exp_rm')  ## choose which top genes for plotting 
rdata_keyword = 'mixed_model_results_exp_corrected.Rdata'  # choose which rdata to get from: corrected for study only

for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    ## the cell marker folder
    if(disease == 'AD'){
        (folder=max(grep('201.*Hippocampus', list.dirs('../config/cell_type_markers/',recursive = T, full.names = T), value = T)))
    }else{
        (folder=max(grep('201.*Striatum', list.dirs('../../config/cell_type_markers/',recursive = T, full.names = T), value = T)))
    }
    
    ## get all the cell type markers
    f_marker_ls <- grep(paste0(cell_type_ls,  collapse = '.tsv|'), list.files(folder, full.names = T), value = T)
    for(f in f_marker_ls){ ## loop 2 for each cell type files
        
        ## get the gene list from the cell type markers
        df <- read.delim(f, comment.char = '#')
        cell_type <- as.character(levels(df$cell_type))
        print(cell_type)
        gene_list <- sort(as.character(df$geneSymbol))
        gene_list <- grep('\\|', gene_list, invert = T, value = T) ## rm multiple genes
        
        
        ## where the mixed model result files
        mm_rdata_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/') ## where the rdata
        mm_dir <- paste0(disease_dir,result_rank,'/', model,model_keyword,'/') ## where the ranked results
        # plot out dir (a sub folder of 'gene_heatmaps' will be created)
        plt_out <- paste0(home_dir,'ND_results/gene_heatmaps_explore/cell_markers/', disease, '/', Sys.Date(), '/', cell_type, '/')
        
        
        
        plotThesisHeatGenes(mm_dir,mm_rdata_dir,
                            plt_out,
                            phase_ls, 
                            threshold,
                            rdata_keyword = rdata_keyword,
                            rm_gene_list = rm_gene_list,
                            top_genes =top_genes, 
                            opposite_phase =opposite_phase, 
                            gene_list =gene_list,
                            row_width =25,
                            legend_ls =legend_ls, cluster_rows =T, cluster_cols =F,plot_indi=F)
        
    } ## loop 2
}##loop1




cat(" # updated 2017-02-20
    #**************************#
    # PART 11.5.1 plot the cell population changes for all diseases for indi study
    #**************************#\n")
#' must run part 3.3.2a and 3.3.2b

rm(list=setdiff(ls(),'home_dir'))
source('MGP_estimation/estimate_cell_population_summary_and_plots.R')
source('config_wrappers.R')
phase_ls =c('early', 'late')


(f_out_dir = paste0(home_dir, '/ND_results/cell_population_plots_indi_studies/', Sys.Date(), '/'))

font_size <- 30
x_angle <- 90  ## rotation of x labels
plot_indi <- F ## whether to plot individual studies in a separate plot
for(disease in disease_ls){
    source('config_wrappers.R')
    df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_all.tsv')
    df <- mainPlotIndiStudyCellEstimates(disease, phase_ls, disease_dir,df_info, original_genotype =F,
                                         write_df =T,
                                         f_out_dir = f_out_dir,
                                         plot_indi= plot_indi,
                                         font_size=font_size,
                                         x_angle=x_angle)
}
# source('summary_tables/cell_pop_rotation_rdata.R')

cat("
    #**************************#
    # PART 11.5.2 plot the cell population changes for all diseases (input in mixed model cell types)
    #**************************#\n")

#' after all the estimate cell populations
#' get a summary table for all and plots

rm(list=setdiff(ls(),'home_dir'))
source('MGP_estimation/estimate_cell_population_summary_and_plots_for_disease.R')
source('config_wrappers.R')

phase_ls =c('early', 'late')
geno_f = '../configs/AD_HD_samples_model_types.tsv'
# font_size <- 30
x_angle <- 0 ## rotation of x labels
one_plot_font_size = 14  ## thesis plot font size
poster_font_size = 24  ## poster font size


## output of the plots and table
in_dir = paste0(disease_dir, '/MGP_estimation/all_sample_estimation/')

out_dir = paste0(disease_dir, '/MGP_estimation/all_sample_estimation/plots/')


cellPopPlots(disease_ls, phase_ls, in_dir, out_dir,x_angle = x_angle,one_plot_font_size = one_plot_font_size,
             geno_f, violin =F, outlier_p = F,plot_type = 'png', poster_font_size =poster_font_size, 
             outlier_rm_from_box = F,poster =T)
# cellPopPlots(disease_ls, phase_ls, in_dir, out_dir,x_angle = x_angle,one_plot_font_size = one_plot_font_size,
#              geno_f, violin =F, outlier_p = F,plot_type = 'svg', poster_font_size =poster_font_size, 
#              outlier_rm_from_box = F,poster =T)


cat("
    #**************************#
    # PART 11.5.3 get the correlation for each jackknife run to other runs
    #**************************#\n")


rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('thesis_stuff/jackknife_run_correlation.R')

#model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder
model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'  ## to specify which model folder
disease_ls =c('AD', 'HD')
out_folder=paste0(home_dir, '/ND_results/jackknife_correlations_for_runs/', Sys.Date(), '/')

for(disease in disease_ls){
    input_folder = paste0(home_dir, '/', disease, '_mouse_model_project/mixed_model_jackknife/')
    jackknifeRunCorr(disease, model_keyword,input_folder, out_folder)
}


cat("
    #**************************#
    # PART 12 thesis plots thesis tables
    #**************************#\n")
# rm(list=setdiff(ls(),'home_dir'))
#  
# mm_dir <- paste0(home_dir, "/ND_results/mm_results/2017-02-02/")
# mm_adj_dir <- paste0(home_dir, "/ND_results/mm_results_cell_adj/2017-02-02/")
# # # make the all_ranks and all_ranks_adj
# source('./thesis_stuff/mk_all_ranks_all_ranks_adj.R')
# save(all_ranks, all_ranks_adj,file = '../../results/ranks_paths_2017-02-02.Rdata')

##############
## summary cell marker rotation estimations and plots
# saved r objects: cell_marker_freq, cell_marker_rotations, cell_marker_msg
##############
#'summary of cell markers and estimation of rotations and variations
#' results in /results/ND_results/DE_genes_markers/
source('thesis_stuff/check_cell_pop_rotation.R')
#+++++++++++++++++++
## get which DE genes are also cell type markers before and after correction
#+++++++++++++++++++
## thesis tables /results/ND_results/DE_genes_markers/'


## see check_genes.R for plots and tables


source('./thesis_stuff/check_genes_helpers.R')
####################################
#### check correlations of ranks early and late
####################################
df_all_el <- NULL
source('config_wrappers.R')
for (disease in disease_ls){
    print('compare cell early and late correclations ')
    df <- compareELCor(all_df = all_ranks_adj, disease, prefix = 'adj')
    df2 <- compareELCor(all_df = all_ranks, disease, prefix = '')
    if(is.null(df_all_el)){
        df_all_el <- rbind(df, df2)
    }else{
        df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
    }
}
## check top genes early late correction
top_t_ls = c(100, 200, 500)
for (disease in disease_ls){
    print('compare cell early and late correclations ')
    for(regulation in c('up', 'down')){
        for(top_t in top_t_ls){
            df2 <- compareELCor(all_df = all_ranks[which(all_ranks[paste0(regulation, '_jack')] <= top_t), ], disease, prefix = paste0('top_', regulation, top_t))
            df <- compareELCor(all_df = all_ranks_adj[which(all_ranks_adj[paste0(regulation, '_jack')] <= top_t), ], disease, prefix = paste0('adj_top_', regulation, top_t))
            df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
        }
    }
}

# ### compare between AD and HD
# for(phase in c('early', 'late')){
#     print('compare AD and HD ')
#     df <- compareADHDCor(all_df = all_ranks_adj, disease,  phase, prefix = 'adj')
#     df2 <- compareADHDCor(all_df = all_ranks, disease,  phase, prefix = '')
#     
#     df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
#     
# }




outdir <- pasete0(home_dir, '/ND_results/early_late_corr/')
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, Sys.Date(), "_early_late_corr.tsv")
writeTable(df_all_el, f_out = f)
################################

threshold <- 0.05
df <- upDownCount(all_ranks, threshold)
df$Count <- paste0(df$Freq, '(', round(df$ratio*100,2), '%)')
df <- df[,c('disease','phase', 'regulation', 'Freq','Count')]
df$correction <- 'Before'
df_before <- df

df <- upDownCount(all_ranks_adj, threshold)

df$Count <- paste0(df$Freq, '(', round(df$ratio*100,2), '%)')
df <- df[,c('disease','phase', 'regulation', 'Freq','Count')]
df$correction <- 'After'
df_after <- df

df <- rbind(df_before, df_after)
df <- left_join(df, df_all_count)
df <- as.data.frame(unclass(df))
df$correction <- factor(df$correction, levels= c('Before', 'After'))

plot_dir <- pasete0(home_dir, '/ND_results/figures/pvalues/')
dir.create(plot_dir, recursive = T,showWarnings = F)
writeTable(df, f_out = paste0(plot_dir, Sys.Date(), '_padj5_count.tsv'))

########
## thesis bar plot DE genes before and after cell type correction
########
source('config_wrappers.R')
for(disease in disease_ls){
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar.png') 
    plotBarDE(df2, one_plot_font_size=10, f_out)
}

for(disease in disease_ls){
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar.svg') 
    plotBarDE(df2, one_plot_font_size=10, f_out)
}





#+++++++++++++++++++
## thesis tables: get which DE genes are also cell type markers before and after correction
#+++++++++++++++++++
## thesis tables /results/ND_results/DE_genes_markers/'

#+++++++++++++++++++
##### get the enriched GO terms and cross disease overlap terms -thesis
#+++++++++++++++++++

#+++++++++++++++++++
# compare my top genes and the reported genes (thesis tables)
#+++++++++++++++++++

#+++++++++++++++++++
# get the gender gene non-expression median (reason for the threshold) # must use the expression.R in the 
# /mixed_model/random_intercept_include_NA/
# source('thesis_stuff/sex_gene_threshold.R')
#+++++++++++++++++++

#+++++++++++++++++++ thesis table
#' summary of all samples and genes counts for input data
#' genes are after filter: input data from mixed_model/random_intercept_include_NA_low_exp_rm/
source('summary_tables/summary_mm_input_data.R')
#+++++++++++++++++++

#+++++++++++++++++++
##### get the enriched GO terms and cross disease overlap terms -thesis
#' 
#+++++++++++++++++++

#+++++++++++++++++++
########## thesis plot
####### summary of each mouse model is in how many studies
#+++++++++++++++++++

#+++++++++++++++++++ 
# 2017-03-01 # thesis plots
# plot the expression of a gene (or gene list) from each study
# with raw expression, study corrected or study and MGP corrected
source('thesis_stuff/plot_indi_gene_exp_per_study.R')
#+++++++++++++++++++ 
f_out_dir <- paste0(home_dir, '/ND_results/gene_expr_before_after_MGP/', Sys.Date(), '/')
dir.create(f_out_dir, showWarnings = F,recursive = T)

#++
disease <- 'AD'
phase <- 'late'
#++
x <- arrayData(disease,phase)
array_dat_raw <- x[[1]]
array_dat_MGP <- x[[2]]
array_dat_study <- x[[3]]
array_design <- x[[4]]

### compare before and after MGP ( 1 gene a plot) ## plot
gene_ls <- c('Trem2','Msmo1')
mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir)
# 
# 
# #++
# disease <- 'HD'
# phase <- 'late'
# #++
# x <- arrayData(disease,phase,kuhn =T)
# array_dat_raw <- x[[1]]
# array_dat_MGP <- x[[2]]
# array_dat_study <- x[[3]]
# array_design <- x[[4]]
# 
# ### compare before and after MGP ( 1 gene a plot) ## plot
# gene_ls <- c('Ddit4l', 'Nrep')
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,one_plot_font_size =9,y_axis_null =F)
# 
# 
# 
