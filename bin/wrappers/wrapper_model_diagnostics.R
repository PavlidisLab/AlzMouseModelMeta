## 2017-02-16: look at the model diagnostic for before and after marker gene profiles (MGP) correction 
#' MGP correct for cell types, degree of freedom will decrease 


 
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/mixed_model.R')
source('mixed_models/mixed_model_adj_celltype.R')

model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder

disease <- c('AD')
source('config/config_wrappers.R')
model <- c('random_intercept')
phase <- c('late')
# phase <- c('early')

genes =c( 'Msmo1','Sqle')

plot_out <- paste0('../../results/ND_results/mixed_model_diagnostics/', disease, '/', Sys.Date(), '/')
dir.create(plot_out, showWarnings = F, recursive = T)


### load the marker dataframes from adj results
data_dir <- paste0(disease_dir,'/mixed_model/random_intercept_include_NA_low_exp_rm_adj_cell_pop/') ## after MGP correction to get df_markers
(marker<- paste0(data_dir, phase, '/mixed_model_results.Rdata'))
load(marker)

if(disease == 'AD'){
    cell_types = c("Astrocyte" ,"DentateGranule", 'GabaSSTReln',"Microglia", 'Oligo', 'Pyramidal_Thy1',
                   "GABAergic",  "Pyramidal")
}else if(disease =='HD'){
    cell_types = c("Astrocyte","Cholinergic" ,"Microglia" ,"Spiny", 'Oligo','ForebrainCholin')
}
## Cholinergic" = ForebrainCholin'
cell_types =intersect(cell_types,colnames(df_markers))
print(paste0('Input cell types are ', paste0(cell_types, collapse = ', ')))


## load array data before and after correction
## input of expressions are the same for before and after
data_dir <- paste0(disease_dir,'/mixed_model/random_intercept_include_NA_low_exp_rm/') ## where the mixed model before MGP correction
(exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
### get the expression data and MM results
print(data_dir)
load(exprdata)


## before correction
(plot_dir1 <- paste0(plot_out, phase,'/beforeMGP/'))
## after correction
(plot_dir2 <- paste0(plot_out, phase,'/afterMGP/'))


for(i in genes){
    print(i)
    #before MGP correction
    mixedModelStudy(array_dat, array_design, i, 
                    model=model,
                    to_plot=T,
                    plot_dir = plot_dir1,
                    full_report = F,
                    REML=FALSE,
                    estimateCI= F,
                    disease_stage_only =T,
                    study_color = NULL)
    
    ## after MGP correction
    anova_r=mixedModelStudyCellType(array_dat, array_design, i, 
                            model=c('random_intercept'),
                            to_plot=T,
                            plot_dir =plot_dir2,
                            full_report = F,
                            REML=FALSE,
                            estimateCI= F,
                            disease_stage_only =T,
                            study_color=NULL,
                            df_markers=df_markers,
                            cell_types=cell_types,
                            return_anova = T)
} # for each gene



