#' 2016-03-15
#' use wrapper: wrapper_filter
#' filter the non-expressed probes by gender genes 
#' (probes lower than the median of the opposite gender probes(i.e. none expressed) are removed)
#' @examples
# r_obj <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/analysis_mode/GSE50521/results/GSE50521_objects.Rdata"
# dataset <- getGSEID(r_obj)
# out_folder <- paste0("/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/explore/analysis_mode/", dataset, "/")
# filterExpArray(r_obj, out_folder)


 
source("./low_exp_rm_analysis/rm_low_exp_by_gender_genes.R")
source("helper_functions.R")

#' @param r_obj: non-filtered array data (R data) and other info saved from explore_data_general.R
#' @param out_folder: output of the new array data and Rdata
#' @return plots and new filtered R data 

filterExpArray <- function(r_obj, out_folder){
    ##--------------------##
    ## define gender genes
    ##--------------------##
    df_gene<- read.table(header=TRUE, text='
                         GenderGenes    GeneSymbols
                         F    Xist
                         M    Kdm5d
                         M    Rps4y1')
    
    ##--------------------##
    ## get the new filtered array_dat, filtered by gender genes and median threshold
    ##--------------------##
    x <- probeExp(r_obj, df_gene)
    df_exp_t <- x[[1]]
    dataset <- x[[2]]
    platform <- x[[3]]
    array_dat_old <- x[[4]]
    df_exp <- x[[5]]
    msg <- x[[6]]
    annotation <- x[[7]]
    array_design <- x[[8]]
    
    y <- findThreshold(method_t="median", df_exp_t, dataset,platform, array_dat = array_dat_old, df_exp, to_plot =F, msg =msg)
    array_dat <- y[[2]]
    msg <- y[[3]]
    
    
    ##--------------------##
    ## plot and save
    ##--------------------##
    ###
    print(paste0("output folder is ", out_folder))
    
    plot_folder <- paste0(out_folder, 'plots/')
    result_folder <- paste0(out_folder, 'results/')
    
    dir.create(plot_folder, showWarnings = FALSE, recursive=T)
    dir.create(result_folder, showWarnings = FALSE, recursive=T)
    
    plot_pre <- paste0(plot_folder, dataset, '_')
    result_pre <- paste0(result_folder, dataset, '_')
    
    
    ##plots 
    plot_pre <- paste0(plot_pre, "none_exp_filtered_")
    prefix <- " (non-expressed probes filtered)"
    ## plot after outlier removed
    plotFunctions(array_dat, plot_pre, dataset, array_design, prefix, width = 1600, height = 1000)
    
    ## save log and robject
    sink(paste0(result_pre, "exp_process_log_low_exp_filtered.txt"), type="output")
    writeLines(msg)
    sink()
    
    ## if the array is already batch corrected, add _objects_after_batch_correction.Rdata to distinguish
    if(all(grepl("after_batch_correction", r_obj))){
        filename <- paste0(result_pre, "objects_after_batch_correction.Rdata")
    }else(filename <- paste0(result_pre, "objects.Rdata"))
    
    save(array_design, annotation, array_dat, dataset, platform, msg,
         file = filename)
}

