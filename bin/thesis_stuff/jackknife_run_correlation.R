## 2017-01-25
#' checking the correlation of gene rankings for each jackknife procedure 
#'

 
library(HelperFunctions)
library(dplyr)

# model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder
# model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'  ## to specify which model folder
# disease ='HD'
# input_folder = paste0(home_dir, '/', disease, '_mouse_model_project/mixed_model_jackknife/')
# out_folder='../../results/ND_results/jackknife_correlations_for_runs/'
# jackknifeRunCorr(disease, model_keyword,input_folder, out_folder)



jackknifeRunCorr <- function(disease, model_keyword,input_folder, out_folder=input_folder){
    (folder = grep(paste0(model_keyword, '$'), list.dirs(input_folder, recursive = F), value = T))
    ## input files 
    (f_ls <- paste0(folder, '/', list.files(folder, pattern = 'jackknife.*.tsv')))
    
    dir.create(out_folder, showWarnings = F)
    f_out <- paste0(out_folder, '/', disease, model_keyword, '_rank_corr_for_runs_', Sys.Date(), '.tsv')
    writeTable(df=NULL, f_out = f_out, msg = paste0('# ', Sys.Date()))
    
    
    for(f in f_ls){
        (phase <- grep('early$|late$', unlist(strsplit(f, split = '_')), value = T))
        df <- read.delim(f, comment.char = '#')
        (msg <- grep('# run', readLines(f, n=20), value = T))
        
        
        #grep the rank columns
        cols <- grep('run.*.rank', colnames(df), value = T)
        
        ## get correlation for every combo of 2 runs
        (combos <- as.data.frame(t(combn(cols,2))))
        
        s_cor_ls <- NULL
        for(n in 1:nrow(combos)){
            s_cor <- cor(df[, as.character(combos[n, 1])], df[, as.character(combos[n, 2])], method = 'spearman')
            s_cor_ls <- c(s_cor_ls, s_cor)
        }
        out_df <- cbind(data.frame(disease = disease), phase = phase, 
                        combos, data.frame(spearman_correlation = s_cor_ls))
        out_df <- orderCol(out_df, cols = 'spearman_correlation')
        
        ## get for the summary of correlations for each run to other runs
        tmp <- reshape2::melt(out_df[, 3:5],measure.vars = c("V1", "V2"))
        tmp$value <- as.factor(tmp$value)
        my.summary <- function(x, na.rm=TRUE){
            result <- c(n=as.integer(length(x)),
                        Mean=mean(x, na.rm=TRUE), SD=sd(x, na.rm=TRUE), 
                        Min=min(x), Max=max(x),Median=median(x))}
        
        tmp_df <- plyr::ddply(tmp, c("value"), function(x) my.summary(x$spearman_correlation))
        tmp_df <- orderCol(tmp_df, cols = 'Median')
        
        writeTable(df=out_df, f_out = f_out,
                   msg = paste0('\n# INPUT files ', f,'\n', paste0(msg, collapse = '\n')),
                   file_append = T)
        
        writeTable(df=tmp_df, f_out = f_out,
                   msg = paste0('\n# summary of correlations of each run'),
                   file_append = T)
    }
    print(f_out)
}


