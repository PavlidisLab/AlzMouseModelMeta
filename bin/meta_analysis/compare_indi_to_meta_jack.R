# 2016-05-24
# updated 2016-08-15

## for ranking comparison between rankings in the individual ranking (contained in the meta file) to meta-ranking, jackknife_rangking
## and mixed model ranking
#***************************************************************************#
#---------------------------------------------------------------------------#
#Plot the meta/jack ranking against each study, calculate spearman correlation
#---------------------------------------------------------------------------#

#' study (for each timepoint and genotype base on the dataset_info.tsv), plot the ranking of meta/jack of the genes
#' against the ranking in the study (scatter plot with smoothing)
#' 
#' calculate the spearman correlation
#' must use meta gene file for all individual dataset ranking (for jackknife and mixed model results, which dont contain
#' individual dataset results)
#***************************************************************************#"



##**********************************************************
## FUNCTION
##**********************************************************

#' @param df_target df for meta or jack with indi pvalues for dataset
#' @param prefix: for plot title and save plot
#' @param plot_out: save plot folder
#' @param return_df: return the all df

compareRanksMetaIndi <- function(df_target, prefix,plot_out, return_df=F, plot =T){
    ## grep the target col (Fisher or 'adj_combined_max_p') and GSE p values
    (dataset_ls <- grep('GSE', colnames(df_target), value = T))
    (target_col <- grep('Fisher$|adj_combined_max_p$|up_pval$|down_pval$', colnames(df_target), value = T))
    df_target_t <- df_target[, c(target_col, dataset_ls)]
    ## rm NA and give the ranks for the fisher and dataset respectively
    df_target_t <- na.omit(df_target_t)
    (target_rank <- paste0('rank_', target_col))
    df_target_t[, target_rank] <- rank(df_target_t[, target_col])
    
    ## rank all the dataset in the meta file/jack file
    for(j in 1:length(dataset_ls)){
        (dataset <- dataset_ls[j])
        eval(parse(text= paste0('df_target_t$rank_', dataset,' <- rank(df_target_t[, dataset])')))
    }
    
    ## get the indi pvalue ranks of each GSE
    (rk_ls <- grep('rank.*GSE', colnames(df_target_t), value = T))
    
    ## get the spearman correlation for each dataset
    sp_cor_ls <- NULL
    for (rk_col in rk_ls){
        (sp_cor <- cor(df_target_t[, target_col], df_target_t[, rk_col], method = 'spearman'))
        print(paste0(rk_col, ' ', sp_cor))
        sp_cor_ls<- c(sp_cor_ls, sp_cor)
    }
    
    
    ## plot the smooth line dataset vs meta
    
    if(plot){
        graphics.off()
        png(filename=paste0(plot_out, prefix, '_ranking_cor_comparison.png') , width = 800, height = 600)
        for (i in 1:length(rk_ls)){
            y.loess <- loess(df_target_t[, target_rank] ~ df_target_t[, rk_ls[i]], span=0.3)
            y.predict <- predict(y.loess, df_target_t[, target_rank])
            if(i==1){
                plot(df_target_t[, target_rank], y.predict,type = 'l',
                     xlab = 'Ranks',
                     ylab = 'Dataset Ranks',
                     col = i,
                     ylim =c(0, max(df_target_t[, target_rank])),
                     main = paste0(prefix, '\ncomparison to gene rankings of individual studies'))
            }else{
                lines(df_target_t[, target_rank], y.predict, col=i)
            }
        }
        legend(x="bottomright", legend = paste0(dataset_ls,': ', round(sp_cor_ls, 3)), col=1:length(rk_ls), lty=c(1,1))
        dev.off()
    }
    
    if(return_df){
        return(df_target)
    }
}


##**********************************************************
## MAIN FUNCTION
##**********************************************************
mainCompareRanksMetaIndi <- function(f_meta, f_jack, plot_out,plot_meta =T){
    dir.create(plot_out, recursive = T, showWarnings = F)
    print(f_meta)
    print(f_jack)
    
    ## meta
    df_target <- read.delim(f_meta, comment.char = '#')
    (prefix <- grep('.tsv', unlist(strsplit(f_meta, split ='/')), value = T))
    (prefix <-gsub('.tsv', '', prefix))
    
    ## plot for meta
    print('meta ranks')
    df_meta <- compareRanksMetaIndi(df_target, prefix,plot_out, return_df=T, plot = plot_meta)
    
    ## jack or mixed model
    df_jack <- read.delim(f_jack, comment.char = "#")
    (prefix <- grep('.tsv', unlist(strsplit(f_jack, split ='/')), value = T))
    (prefix <-gsub('.tsv', '', prefix))
    
    ## get the indi studied pvalues
    (dataset_ls <- grep('^GSE', colnames(df_meta), value = T))
    df_jack_t <- df_jack
    # for jack or mixed model
    (index_col <- intersect(colnames(df_jack_t), c('up_pval', 'down_pval', "adj_combined_max_p")))
    df_jack_t <- noWarnings(left_join(df_jack_t[, c("geneSymbol", index_col)], 
                                      df_meta[, c('geneSymbol', dataset_ls)]))
    
    ## rm NA and give the ranks for the fisher and dataset respectively
    df_jack_t <- na.omit(df_jack_t)
    
    ## plot jack or mixed model
    print('jackknife or mixed model')
    compareRanksMetaIndi(df_target = df_jack_t, prefix,plot_out, return_df=F)
}


