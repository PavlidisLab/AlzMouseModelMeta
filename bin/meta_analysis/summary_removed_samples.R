#' 2016-06-13
#' summary of outliers and subset 
#' 
#' @example 
#  
# source('config/config_wrappers.R')
# source('meta_analysis/summary_removed_samples.R')
# data_folder <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
# f_out <- paste0(disease_dir, 'results/limma_DE_summary/removed_samples.tsv') 
# summaryOutlierResult(data_folder, f_out)

library('HelperFunctions')
summaryOutlierResult <- function(data_folder, f_out){
    #' data_folder is the explore folder (before limma) that contains the sample removed .tsv tables
    f_ls <- list.files(data_folder, recursive =T, pattern ='_removed.tsv', full.names = T)
    df_all <- NULL
    for(i in f_ls){
        df <- read.delim(i, comment.char = '#')
        if(is.null(df_all)){
            df_all <- df
        }else{
            df_all <- noWarnings(rbindAllColumns(df_all, df))
        }
    }
    
    writeTable(df_all, f_out=f_out)
    print(paste0('sample removed summary: ', f_out))
}

