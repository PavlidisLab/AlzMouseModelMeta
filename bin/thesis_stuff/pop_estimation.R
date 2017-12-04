# to_run for AD and HD added new sample 2017-01-31 (in test folder)


## recorded done the cell type markers for future heatmaps, save in git/doc/
library(markerGeneProfile)
library(ogbox)
source('cell_population_markers.R')

rm(list=setdiff(ls(),'home_dir'))
 
source("estimate_cell_population.R")



#for A1 and A2 markers
##__
## get A1 A2 markers
df <- read.delim('../../doc/A1_A2_non_activated_astrocyte_markers.tsv', comment.char = '#')
df <- filterContain(df, column ='regulation', value = 'up')

#setdiff(df$geneSymbol, all_ranks$geneSymbol)

# mouseMarkerGenes$Striatum$Astrocyte= setdiff(mouseMarkerGenes$Striatum$Astrocyte, c("Slc1a2","Slc1a3"))
# print(paste0('removed "Slc1a2","Slc1a3" from markers'))

(A1_markers <- as.character(df$geneSymbol[which(df$cell_type =='A1_reactive')]))
(A2_markers <- as.character(df$geneSymbol[which(df$cell_type =='A2')]))

intersect(mouseMarkerGenes$Striatum$Astrocyte, union(A1_markers, A2_markers))
intersect(mouseMarkerGenes$Striatum$Microglia_activation, A1_markers)
intersect(mouseMarkerGenes$Striatum$Microglia_deactivation, A2_markers)


intersect(mouseMarkerGenes$Hippocampus$Astrocyte, union(A1_markers, A2_markers))
intersect(mouseMarkerGenes$Hippocampus$Microglia_activation, A1_markers)
intersect(mouseMarkerGenes$Hippocampus$Microglia_deactivation, A2_markers)


mouseMarkerGenes$Striatum$A1 = setdiff(A1_markers, A2_markers)
mouseMarkerGenes$Striatum$A2 = setdiff(A2_markers, A1_markers)
mouseMarkerGenes$Hippocampus$A1 = setdiff(A1_markers, A2_markers)
mouseMarkerGenes$Hippocampus$A2 = setdiff(A2_markers, A1_markers)


## use Lilah's cholin markers
newcholin = c("Ecel1","Slc18a3","Crabp2","Zic3","Trpc3","Slc5a7","Lhx8") ## 2016-11-01
mouseMarkerGenes$Striatum$ForebrainCholin = newcholin
##----



disease_ls =c('HD', 'AD')
phase_ls = c('early', 'late')


disease_ls =c('HD')
phase_ls = c('early', 'late')


# disease_ls =c('AD')
# phase_ls = c('late')

all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Genotype'  # design_group = 'Original_genotype'
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT

threshold = -10 # filter probes exp less than threshold (input is study corrected value, and low expression filtered)
mixed_model =T
removeNegatives = F  ## ogan suggested use F, 


for (disease in disease_ls){
    if(disease == 'AD'){
        genes =mouseMarkerGenes$Hippocampus
    }else{
        genes =    mouseMarkerGenes$Striatum 
    } 
    
    for(phase in phase_ls){
        source('config/config_wrappers.R')
        (file_ls <- paste0(home_dir, '/ND_results/cell_population/all_sample_estimation_test/',disease, '/', phase,'/'))
        
        
        (r_ob <- list.files(file_ls, full.names=T, recursive=T, pattern= "mixed_model_results_exp_corrected.Rdata"))
        (outDir <-paste0(file_ls,Sys.Date(), '/'))
        dir.create(outDir,showWarnings=F, recursive=T)
        
        cellPopEstimatePlot(r_ob, outDir, design_group = design_group, mixed_model = mixed_model,
                            sigTest = sigTest,
                            wt_only = wt_only,
                            removeNegatives = removeNegatives, threshold = threshold,genes = genes)
    }
}

cat("
    #**************************#
    # PART 11.5.2 plot the cell population changes for all diseases (input in mixed model cell types)
    # saved in 'ND_results/cell_population/all_sample_estimation/plots/'
    #**************************#\n")
#' must run part 3.3.2c (study intercept corrected)
#' after all the estimate cell populations
#' get a summary table for all and plots, saved in ''git/results/ND_results/cell_population_plots/'
 
rm(list=setdiff(ls(),'home_dir'))
source('estimate_cell_population_summary_and_plots_for_disease.R')

disease_ls = c('AD', 'HD')
phase_ls =c('early', 'late')
geno_f = '../../doc/AD_HD_samples_model_types.tsv'
# font_size <- 30
x_angle <- 0 ## rotation of x labels
one_plot_font_size = 14  ## thesis plot font size
one_plot_font_size = 24  ## poster font size
in_dir = paste0(home_dir, '/ND_results/cell_population/all_sample_estimation_test/')
out_dir = in_dir

cellPopPlots(disease_ls, phase_ls, in_dir, out_dir,x_angle = x_angle,one_plot_font_size = one_plot_font_size,
             geno_f, violin =F, outlier_p = F,plot_type = 'png', poster =T)