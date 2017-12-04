#' 2016-02-23
#' modified: 2016-04-25
#' use wrapper
#' @note
#'  need f_df defined(where the meta and jackknife files are )
# f_df <- "/home/bzhuang/AD_mouse_model_project/low_exp_rm_analysis/meta_analysis_lrm/"
#' plot the rank correlation (spearman cor) of the p values input for each of the meta analysis
#' to see which p values from which study are highly correlated.
#' overall trend: p values from the same study and comparison with the same wild type samples are usually higher correlated.

source("helper_functions.R")


#' @param f_df: the folder contains the meta_genes files
metaHeatmap <- function(f_df, width = 1600, height = 1000){
    f_df_heat <- paste0(f_df, "/spearman_cor_heatmaps/")
    dir.create(f_df_heat, showWarnings=F)
    (f_ls <-grep ("meta_genes.tsv",list.files(f_df, full.names=T), value=T))
    
    for (i in f_ls){
        print(i)
        df <- read.delim(i, comment.char="#")
        df <- na.omit(df[, grep('GSE', colnames(df), value=T)])
        sample_cor_s <- cor(df, method = "spearman")
        
        diag(sample_cor_s) <- NA
        
        ## save plot
        (heatmap_title <- unlist(strsplit(i, split ="/|\\.")))
        (heatmap_title <- paste0(heatmap_title[length(heatmap_title)-1], "_spearman_correlation"))
        f_heatmap <-paste0(f_df_heat, heatmap_title, ".png") 
        png(filename= f_heatmap, width = width, height = height)
        pheatmapNA(sample_cor_s, main = heatmap_title)
        dev.off()
    }
}






