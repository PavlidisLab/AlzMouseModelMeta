
setwd(file.path(here::here(),'bin'))


#****************
## PART 11 thesis tables, figures, rdata for easy comparison ----
#****************


#**************************#
# PART 11.0 [Preperation]Rdata for all jacknife results for all disease-----
#**************************
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





#**************************#
# PART 11.3 plot Heatmaps for top genes----
#**************************

#**************
# PART 11.3.1a [Figure 2]. Top 20 up and down-regulated genes for AD early phase after MGPs correction. ----
# PART 11.3.1b [Figure 3]. Top 20 up and down-regulated genes for AD late phase after MGPs correction. ----
#' OUTPUT: 
#' [early AD]
#' paste0(home_dir,'/ND_results/top_gene_heatmaps/AD/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/*/early/AD_early_Genotype_Model_types.png')
#' [late AD]
#' paste0(home_dir,'/ND_results/top_gene_heatmaps/AD/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/*/AD_late_Genotype_Model_types.png')

# make the heatmap for top jackknife genes MGP corrected
# (top genes of the phase)- with cell population corrected
# and top mixed model genes (not used in thesis)
# expression are corrected to the corresponding gene list (corrected for study or corrected for celltypes)
# also plot the top genes of the opposite phase expression
#**************
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
                    mm_rdata_dir <- paste0(disease_dir, '/mixed_model/', model,mm_rdata_keyword,'/') ## where the expression results
                    mm_dir <- paste0(disease_dir,'/',result_rank,'/', model,mm_rdata_keyword,'/') ## where the ranked results
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



#**************
# PART 11.3.2 [PLOT NOT USED]make the heatmap for cell type marker genes ----
# OUTPUT: paste0(home_dir,'/ND_results/gene_heatmaps_explore/cell_markers/', disease, '/', Sys.Date(), '/', cell_type, '/')
## just to explore, not included in thesis
# input expression is corrected for study only
#**************

#' with model_keyword_ls <- c('_include_NA_low_exp_rm')  ## this correct all the genes if intercept is provided

rm(list=setdiff(ls(),'home_dir'))

## load the cell type markers rdata: variable `mouseMarkerGenes`
# load('../configs/cell_type_markers/2017-02-27/mouseMarkerGenes.Rdata')
load('../configs/mouseMarkerGenesCombined.rda') ## load the markergenes
mouseMarkerGenes = mouseMarkerGenesCombined


source('helper_functions.R')
source('mixed_models/plot_up_down_genes_thesis_figure_corrected_value.R')
source('config_wrappers.R')

legend_ls <- list(c('Genotype', 'Study'), c('Genotype', 'Study'))


model <- c('random_intercept')
model_keyword <- c('_include_NA_low_exp_rm')  
phase_ls <- c('early', 'late')
regulation_ls <- c('up', 'down')
result_rank <- c('mixed_model_jackknife')  ## which rank list to get genes from
opposite_phase = F
top_genes = F

## all cell types (hippocampus and striatum)
cell_type_ls = c('Microglia', 'Astrocyte', 'DentateGranule', 'GabaSSTReln', 'Oligo', 'Pyramidal',
                 'ForebrainCholin', 'Spiny')
mm_rdata_keyword_ls <- c('_include_NA_low_exp_rm')  ## choose which top genes for plotting 
rdata_keyword = 'mixed_model_results_exp_corrected.Rdata'  # choose which rdata to get from: corrected for study only



for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    ## the cell markers for the brain region
    if(disease == 'AD'){
        cell_marker <- mouseMarkerGenes$Hippocampus
    }else{ ## for HD
        cell_marker <- mouseMarkerGenes$Striatum
    }
    
    ## get all the cell type markers that apply
    final_cell_type_ls <- intersect(cell_type_ls, names(cell_marker))
    for(cell_type in final_cell_type_ls){ ## loop 2 for each cell type files
        
        print(cell_type)
        gene_list <- sort(unlist(cell_marker[cell_type]))
        gene_list <- grep('\\|', gene_list, invert = T, value = T) ## rm multiple genes
        
        
        ## where the mixed model result files
        mm_rdata_dir <- paste0(disease_dir, '/mixed_model/', model,model_keyword,'/') ## where the rdata
        mm_dir <- paste0(disease_dir,'/',result_rank,'/', model,model_keyword,'/') ## where the ranked results
        # plot out dir (a sub folder of 'gene_heatmaps' will be created)
        plt_out <- paste0(home_dir,'/ND_results/gene_heatmaps_explore/cell_markers/', disease, '/', Sys.Date(), '/', cell_type, '/')
        
        
        
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



####

#**************************#
# PART 11.5.2a [Figure 1A.] MGPs of neurons in AD mouse models. ----
# PART 11.5.2b [Figure 1B.] MGPs of glial cells in AD mouse models. ----

#' plot the cell population changes (generate a lot of figures but only 2 are used in the paper)

# for all diseases (input in mixed model cell types)
# OUTPUT:
#[Figure 1A] 'files_results/AD_mouse_model_project/MGP_estimation/plots/*[date]/AD_neurons_box.png'
#[Figure 1B] 'files_results/AD_mouse_model_project/MGP_estimation/plots/*[date]/AD_glia_box.png'
#**************************

#' after all the estimate cell populations
#' get a summary table for all and plots

rm(list=setdiff(ls(),'home_dir'))
source('MGP_estimation/estimate_cell_population_summary_and_plots_for_disease.R')
source('config_wrappers.R')
source('helper_functions.R')

phase_ls =c('early', 'late')
geno_f = '../configs/AD_HD_samples_model_types.tsv'
# font_size <- 30
x_angle <- 0 ## rotation of x labels
one_plot_font_size = 14  ## thesis plot font size
poster_font_size = 24  ## poster font size


## output of the plots and table
in_dir = paste0(disease_dir, '/MGP_estimation/')

out_dir = paste0(disease_dir, '/MGP_estimation/plots/')


cellPopPlots(disease_ls, phase_ls, in_dir, out_dir,x_angle = x_angle,one_plot_font_size = one_plot_font_size,
             geno_f, violin =F, outlier_p = F,plot_type = 'png', poster_font_size =poster_font_size, 
             outlier_rm_from_box = F,poster =T)
# cellPopPlots(disease_ls, phase_ls, in_dir, out_dir,x_angle = x_angle,one_plot_font_size = one_plot_font_size,
#              geno_f, violin =F, outlier_p = F,plot_type = 'svg', poster_font_size =poster_font_size, 
#              outlier_rm_from_box = F,poster =T)



#**************************#
# PART 11.5.3 [TABLE NOT USED]get the correlation for each jackknife run to other runs ----
#**************************


rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('thesis_stuff/jackknife_run_correlation.R')

#model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder
model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'  ## to specify which model folder
out_folder=paste0(home_dir, '/ND_results/jackknife_correlations_for_runs/', Sys.Date(), '/')

for(disease in disease_ls){
    input_folder = paste0(home_dir, '/', disease, '_mouse_model_project/mixed_model_jackknife/')
    jackknifeRunCorr(disease, model_keyword,input_folder, out_folder)
}



#**************************#
# PART 12.1 [Preperation] gather ranking results ----
# gene ranks and pvalues
# enrichment results
# save as a Rdata for plots and tables
#**************************

## gather all the jackknife ranking results for before MGP and after MGP (adj),
## mark the genes that are used as cell type markers
## a table for jackknife rankings before MGP (`all_ranks`)
## a table  for jackknife rankings after MGP (`all_ranks_adj`)

rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
(mm_dir <- max(list.dirs(paste0(home_dir, "/ND_results/mm_results/"))))
(mm_adj_dir <- max(list.dirs(paste0(home_dir, "/ND_results/mm_results_cell_adj/"))))
# # make the all_ranks and all_ranks_adj
source('thesis_stuff/mk_all_ranks_all_ranks_adj.R')

## make GO_non_adj and GO_adj for mixed models:
folder_pre = 'ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm'
f_result_out =paste0(home_dir, "/ND_results/ermineJ/", Sys.Date(),'.Rdata')
go_plot_out <- paste0(home_dir, "/ND_results/ermineJ_enrichment/")
source('thesis_stuff/Go_enrichment_summary.R')

save(all_ranks, all_ranks_adj, GO_adj, GO_non_adj, 
     file = paste0(home_dir, "/ND_results/ranks_tables.Rdata"))


#**************
## PART 12.2 [TABLE NOT USED]summary cell marker rotation estimations and plots ----
# saved r objects: cell_marker_freq, cell_marker_rotations, cell_marker_msg
#**************
#'summary of cell markers and estimation of rotations and variations
#' results in /results/ND_results/DE_genes_markers/
#' 
source('config_wrappers.R')

(outdir <- paste0(home_dir, '/ND_results/DE_genes_markers/rotations/', Sys.Date(), '/'))
(plotdir <- paste0(outdir, 'plots/')) ## for heatmaps of cell markers (study corrected value,, per cell type, disease and phase)
(f_rdata_out <- paste0(home_dir, '/ND_results/DE_genes_markers/cell_marker_summary_', Sys.Date(), '.Rdata'))

source('thesis_stuff/check_cell_pop_rotation.R')

#*******************thesis table
##PART 12.3 [Table S3]. Top-ranked cell-type marker genes before and after marker gene profiles correction. ----
#' OUTPUT: paste0(home_dir, '/ND_results/DE_genes_markers/DEmarkers/thesis_DE_markers.tsv'
#' Some formatting needed for a neat table (e.g. change NA to --)
#' get which DE genes are also cell type markers before and after MGP correction
#' 
#*******************

source('./thesis_stuff/check_genes_helpers.R')
(outdir <- paste0(home_dir, '/ND_results/DE_genes_markers/DEmarkers/'))

## need to load all_rank, all_rank_adj
load(paste0(home_dir, "/ND_results/ranks_tables.Rdata"))

threshold = 50

phase_ls =c('early','late')
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, "thesis_DE_markers.tsv")
#if(file.exists(f)){file.remove(f)}
for(disease in disease_ls){
    for(phase in phase_ls){
        print(paste0(disease, '_', phase))
        df_m <- getCellMarkers(disease, phase, threshold)
        noWarnings(writeTable(df_m, f_out = f, file_append =T))
    }
}


#**************************
## PART 12.4.1 [NOT USED]check correlations of ranks early and late ----
# OUTPUT in paste0(home_dir, '/ND_results/early_late_corr/')
#**************************
## need to load all_rank, all_rank_adj
load(paste0(home_dir, "/ND_results/ranks_tables.Rdata"))

df_all_el <- NULL
source('./thesis_stuff/check_genes_helpers.R')

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



(outdir <- paste0(home_dir, '/ND_results/early_late_corr/'))
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, Sys.Date(), "_early_late_corr.tsv")
writeTable(df_all_el, f_out = f)


#**************************
## PART 12.4.2 [Figure S2]. Number of differentially expressed genes (FDR < 0.05) in AD mouse models before and after marker gene profiles correction. ----
# OUTPUT: paste0(home_dir, '/ND_results/figures/pvalues/.[Date]_AD_pvalue_bar_no_dir.png')
#**************************

## need to load all_rank, all_rank_adj
source('./thesis_stuff/check_genes_helpers.R')
source('config_wrappers.R')
load(paste0(home_dir, "/ND_results/ranks_tables.Rdata"))

plot_dir <- paste0(home_dir, '/ND_results/figures/pvalues/')

df_all_count <- NULL
for(phase in c('early', 'late') ){
    for(regulation in c('up', 'down')){
        print('compare cell adj vs non cell adj correclations ')
        df <- compareAdjCor(disease, phase, regulation)
        if(is.null(df_all_count)){
            df_all_count <- df
        }else{
            df_all_count <- rbind(df_all_count,df)
        }
    }
}


threshold <- 0.05
df <- upDownCount(all_ranks, threshold)
df <- df[,c('disease','phase', 'regulation', 'Freq', 'ratio')]
df$correction <- 'Before'
df_before <- df

df <- upDownCount(all_ranks_adj, threshold)


df <- df[,c('disease','phase', 'regulation', 'Freq', 'ratio')]
df$correction <- 'After'
df_after <- df

df <- rbind(df_before, df_after)
df$Count <- paste0(df$Freq, '(', round(df$ratio*100,2), '%)')
df <- left_join(df, df_all_count)
df <- as.data.frame(unclass(df))
df$correction <- factor(df$correction, levels= c('Before', 'After'))

# get no direction
df1 <- aggregate(Freq ~ disease + phase +correction , data = df, FUN = sum)
df2 <- aggregate(ratio ~ disease + phase +correction , data = df, FUN = sum)
df_no_reg_dir <- left_join(df1, df2)
df_no_reg_dir$Count <- paste0(df_no_reg_dir$Freq, '(', round(df_no_reg_dir$ratio*100,2), '%)')


dir.create(plot_dir, recursive = T,showWarnings = F)
writeTable(df, f_out = paste0(plot_dir, Sys.Date(), '_padj5_count.tsv'))
writeTable(df_no_reg_dir, f_out = paste0(plot_dir, Sys.Date(), '_padj5_count_no_dir.tsv'))

## thesis bar plot DE genes before and after cell type correction

## x label is regulation (up or down, facet by phase) - not used
df2 <- filterContain(df, column = 'disease', value = disease)
f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar.png') 
p <- plotBarDE(df2, one_plot_font_size=10, f_out, x_col = 'regulation',
               y_col = 'Freq' ,
               x_label = 'Regulation',
               y_label = 'Count',
               return_p = T)




## Figure S2
df2 <- filterContain(df_no_reg_dir, column = 'disease', value = disease)
f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar_no_dir.png') 
p <- plotBarDE(df2, one_plot_font_size=14, f_out, x_col = 'phase',
               y_col = 'Freq' ,
               x_label = 'Disease Phase',
               y_label = 'DE Gene Count',
               return_p = T, facet_phase =F, save_p = T)



#+++++++++++++++++++
# PART 12.5 [PLOT NOT USED]plot pvalue distribution before and after MGP correction ----
# OUTPUT: paste0(home_dir, '/ND_results/figures/pvalues/.[Date]_AD_pvalue.png')
#+++++++++++++++++++


plot_dir <- paste0(home_dir, '/ND_results/figures/pvalues/')
dir.create(plot_dir, recursive = T,showWarnings = F)

## need to load all_rank, all_rank_adj
load(paste0(home_dir, "/ND_results/ranks_tables.Rdata"))


## plot pvalue distribution
df <- all_ranks[, c('disease', 'phase','pvalue')]%>%droplevels()
df$correction <- 'Before'
df2 <- all_ranks_adj[, c('disease', 'phase','pvalue')]%>%droplevels()
df2$correction <- 'After'
df <- rbind(df,df2)
df$phase <- paste0(df$disease, ' ', df$phase)
df <- as.data.frame(unclass(df))
df$correction <- factor(df$correction, levels= c('Before', 'After'))

for(disease in disease_ls){
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_', disease, '_pvalue.png')  ## or use svg
    plotPvalue(df2, one_plot_font_size =14, f_out)
}




#+++++++++++++++++++
##PART 12.6 [TABLE NOT USED]get the enriched GO terms and the top genes in the gene set ----
#' OUTPUT: paste0(home_dir, '/ND_results/tables/sig_go_terms/2018-.*_sig_go_terms.tsv'
#' 
#' In the paper, top enriched GO terms are specified in the context, without a table
#+++++++++++++++++++

source('thesis_stuff/check_genes_helpers.R')
outdir <- paste0(home_dir, '/ND_results/tables/sig_go_terms/')

# 'GO Term ID':
# 'GO Term Description'
# 'FDR'
# 'Multifunctionality Score'
# 'Number of Genes': number of genes in the gene set, in brackets are number of genes that are top ranked after MGP correction
# 'Top Hits': top hit genes, in brackets are the ranking of the genes

df <- enrichedGOTerms(outdir, GO_adj,phase_ls =c('early', 'late'))


#+++++++++++++++++++ thesis table
# PART 12.7 [Table S4]. Comparisons between top genes and top hits reported in original studies for the late AD phase. ----
#' 
#' Only the table with late AD phase with the top 20 genes is used in Table S4
#' OUTPUT: paste0(home_dir, '/ND_results/tables/compare_to_original_study/AD_late_20_comparison__adj.*.tsv'
#' 
#' Some formatting needed for a neat table: e.g. change 'no expression' to '--' etc
#' 
#' compare my top genes and the reported genes
#' tables are made for both early and late phases.


#+++++++++++++++++++

outdir <- paste0(home_dir, '/ND_results/tables/compare_to_original_study/')
load(paste0(home_dir, "/ND_results/ranks_tables.Rdata"))

source('thesis_stuff/compare_pub_thesis.R')
source('config_wrappers.R')
threshold_ls =c(20,50)
phase_ls = c('early', 'late')

for(phase in phase_ls){
    for(threshold in threshold_ls){
        for(disease in disease_ls){
            
            ## get the reported genes from repo
            (gene_list_dir <- paste0('../configs/original_publication_gene_list/',disease, '/'))
            ## get the most recent results from MGP adj
            keyword='_adj_cell_pop'
            (input_dir <- max(list.dirs(paste0(home_dir,'/ND_results/top_gene_heatmaps/', disease, 
                                               '/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm',keyword, '/'), recursive = F)))
            (f_exp_adj <- grep(paste0(disease, '_', phase, '_inputdata.tsv'), 
                               list.files(input_dir, recursive = T, full.names = T),
                               value = T))
            
            print(paste0('INPUT file: ', f_exp_adj))
            y <- comparePublication(df = all_ranks_adj, threshold, phase, disease, gene_list_dir, prefix ='_adj', 
                                    out_dir = outdir,
                                    write_out = T, f_exp = f_exp_adj)
        }
    }
}





#+++++++++++++++++++ thesis table
#PART 12.8 [Table 1]. Summary of selected gene expression profiling studies for AD ----
# OUTPUT: paste0(home_dir, '/ND_results/tables/input_data_summary_', Sys.Date(), '.tsv')
#' summary of all samples and genes counts for input data for early and late phase
#' Need manual summary of the total counts
#' 
#' # of genes are after filter: input data from mixed_model/random_intercept_include_NA_low_exp_rm/
#+++++++++++++++++++
source('config_wrappers.R')

f_out <- paste0(home_dir, '/ND_results/tables/input_data_summary_', Sys.Date(), '.tsv')
source('summary_tables/summary_mm_input_data.R')





#+++++++++++++++++++ thesis plot
# PART 12.9.1 [Figure 4A]. Expression of Trem2 in the late AD phase before (top row) and after (bottom row) MGPs correction. ----  
# PART 12.9.2 [Figure 4B]. Expressions of Msmo1 in the late AD phase before and after MGPs correction. ----

# OUTPUT: Fig 4A: paste0(home_dir, '/ND_results/gene_expr_before_after_MGP/.[date]/AD_late_Trem2_mgp_Study.png')
#'This plot need to add red boxes to highlight: Studies (data sets) that are highlighted by a red box have 
#' reported the gene as top hit in the original publications. 

# OUTPUT: Fig 4B: paste0(home_dir, '/ND_results/gene_expr_before_after_MGP/.[date]/AD_late_Msmo1_mgp_Study.png')

# plot the expression of a gene (or gene list) from each study (Trem2 and Msmo1 in late AD are selected as examples in the paper)
# with raw expression, study corrected or study and MGP corrected
source('thesis_stuff/plot_indi_gene_exp_per_study.R')
#+++++++++++++++++++ 
f_out_dir <- paste0(home_dir, '/ND_results/gene_expr_before_after_MGP/', Sys.Date(), '/')
dir.create(f_out_dir, showWarnings = F,recursive = T)

#++
disease <- 'AD'
phase <- 'late'
#++
## the design file that have the Study and the Model_types
(geno_f <- list.files(paste0('../configs/',disease, '_mouse_dataset_doc/'), 
                pattern = paste0('dataset_info_genotypes.tsv'),
                full.names = T))


x <- arrayData(disease,phase,kuhn =F, geno_f=geno_f)
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

#+++++++++++++++++++ thesis table
#PART 12.10 [Table 2]. Summary of selected AD mouse model studies. ----
# This table is manually produced.
#+++++++++++++++++++

#+++++++++++++++++++ thesis table
#PART 12.11 [Table S1]. Summary of Alzheimer’s disease mouse models analyzed. ----
# This table is manually produced.
#+++++++++++++++++++

#+++++++++++++++++++ thesis table
#PART 12.12 [Table S2]. Samples and data sets removed from analysis. ----
# This table is manually produced.
#+++++++++++++++++++



#+++++++++++++++++++ thesis table
# Part 12.13.1 [Table S5]. Top 50 up-regulated genes for AD mouse models in the early phase after marker gene profiles correction. ----
# Part 12.13.2 [Table S6]. Top 50 down-regulated genes for AD mouse models in the early phase after marker gene profiles correction. ----
# OUTPUT: paste0(home_dir, '/ND_results/tables/top_genes/top_genes_cell_adj/jackknife_early_top_50_genes_.*.tsv')
# the output table contains both top 50 up and top 50 down-regulated genes
#+++++++++++++++++++

#' run the `wrapper_for_all_results.R` run part 8.6.1
#**************************#
# PART 8.6.1 save the results of top genes (pvalue is raw p, not up or down p)
# mixed model and jackknife results :corrected for cell types
#**************************#

paste0(home_dir, '/ND_results/mm_results/', Sys.Date(), '/')  ## save the r data
threshold =50
PART 6D save the results of top genes
