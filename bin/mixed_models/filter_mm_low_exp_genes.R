
#' 2016-08-16

#' must run after the mixed_model_results is obtained
#' take previous result, mixed_model_results.Rdata, and filter out genes in the result variable (df_all), and 
#' return a new mixed_model_results.Rdata, with df_all filtered but array_dat, array_design remain the same (not filtered)
#' copy the result of expression.Rdata and mixed_model_before_results.Rdata to the result dir
#' 
#' 
#' @param robj_in: input dir of mixed_model_results.Rdata
#' @param result_dir, output of filtered df_all result (the array dat is not filtered, only the mixed model
#' result df_all is filtered)
#' @examples 
# robj_dir <- "/home/bzhuang/ND_project_combined/PD_mouse_model_project/mixed_model/random_intercept_include_NA/early/"
# result_dir <- "/home/bzhuang/ND_project_combined/PD_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/early/"
# filter_method ='median'
# affy_threshold =6
# 
# filterMMLowExprGenes(robj_dir,result_dir, affy_threshold = affy_threshold,
#                      filter_method = filter_method, return_df = F)
# df <- filterMMLowExprGenes(robj_dir,result_dir, affy_threshold = affy_threshold,
#                            filter_method = filter_method, return_df = T)

library(HelperFunctions)
filterMMLowExprGenes <- function(robj_dir, result_dir, affy_threshold = 6, 
                               filter_method =c("median", "max", "min",  "mean"), return_df = F){
    ## grep the old rdata and copy to result_dir
    (robj <- grep('expression.Rdata|mixed_model_before_results.Rdata', list.files(robj_dir, full.names = T), value = T))
    dir.create(result_dir, recursive = T,showWarnings = F)
    file.copy(robj, result_dir)
    
    
    load(grep('expression.Rdata', list.files(robj_dir, full.names = T), value = T))
    ## load the mixed_model_results.Rdata
    load(grep('mixed_model_results.Rdata', list.files(robj_dir, full.names = T), value = T))
    
    ## get the index based on filter_method
    filter_method <- filter_method[1]
    assign("index_exp", eval(parse(text = paste0("apply(as.matrix(array_dat),1,",filter_method, ")"))))
    
    index_exp <- apply(as.matrix(array_dat),1,function(x) median(x, na.rm = T))
    
    
    ## remove low exp in array dat
    df <- array_dat[which(index_exp > affy_threshold), ]   #(for Affy)
    
    genes <- row.names(df)
    df_all <- df_all[genes, ]%>%droplevels()
    
    # attach expression value
    df_all_expression <- mapBindAllColumns(df_all, df)
    
    msg <- paste0(msg, '\n', Sys.Date(), ': lowly expressed genes are removed when ', 
                     filter_method, ' is less than ',
                     affy_threshold,
                     '\n   After filter gene count: ', nrow(df_all))
    array_dat <- df
    save(array_dat,array_design, df_all, df_all_expression, msg, file=paste0(result_dir, "mixed_model_results.Rdata"))
    save(array_dat_not_qn, array_dat, array_design, phase, msg, plot_dir,all_array_dat_bf_aggregation, all_df_gene_match,
         file=paste0(result_dir,'expression.Rdata'))
    print(result_dir)
    if(return_df){
        return(df_all)
    }
}







