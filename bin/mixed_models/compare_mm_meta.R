# 2016-06-21
# updated 2016-08-02
# compare mixed model results to metaanalysis results for both early and late phase
# save the mm files (with up and down, for ermineJ enrichment)

#' @examples
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease ='AD'
# source('config/config_wrappers.R')
# source('mixed_models/compare_mm_meta.R')
# (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
#                           paste0(disease_dir, 'meta_analysis/meta_jack/')))
# mm_dir = paste0(disease_dir, 'mixed_model/random_intercept_include_NA_low_exp_rm') ## mixed model dir(parent dir)
# mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls, threshold = 200)

library(HelperFunctions)
source('helper_functions.R')
#---------------------------------------------------------------------------#
# functions
#---------------------------------------------------------------------------#

## split the result for one sided p-value(not FDR) based on test stat (p-value for up or)
## - to avoid when a gene is up in one study and down in another
## 1-sided pvalues and adj pvalues are added to the toptable

## down regulated genes(for each probe) and BH adjust the 1 sided p value
splitPvalUp <- function(pval, teststat){
    # 1 sided pvalue for upregulation probes
    # return 1 sided p value
    if(teststat>0){
        temp <- pval/2
    }else{
        temp <- 1-pval/2
    }
    return(temp)
}

splitPvalDown <- function(pval, teststat){
    # 1 sided pvalue for downregulation probes
    # return 1 sided p value
    if(teststat < 0){
        temp <- pval/2
    }else{
        temp <- 1-pval/2
    }
    return(temp)
}

########
## compare to Fisher's combined probability method (jack knife and meta)
########

compareMM <- function(jack_meta_folder, mm_dir, threshold = 200, phase_ls =c('early', 'late'), adj_method = 'BH',
                      compare_fisher =T){
    ## input a jack_meta_folder, and a mixed model dir (has subfolder early, late)
    ## output a stat for all phases, up and down regulation, meta and jack files
    ## and plot xy of rank meta and rank mixed model
    
    stat_df <- NULL
    for (phase in phase_ls){ ## loop 1
        print(phase)
        (datadir = paste0(mm_dir,'/', phase, '/'))
        plotdir = paste0(datadir, '/comparison_to_meta/', Sys.Date(),'/')
        dir.create(plotdir, recursive = T, showWarnings = F)
        
        ## log input notes:
        writeTable(df=NULL, 
                   msg = paste0(Sys.Date(), '\nInput meta files: ', jack_meta_folder,'\nInput mixed model files: ', datadir),
                   f_out = paste0(plotdir,'input_info.txt'))
        
        ## load mm r object results
        (expr_Robject <- paste0(datadir,'mixed_model_results.Rdata' ))
        print(expr_Robject)
        load(expr_Robject)
        
        ## remove CI_sig column ( for ermineJ format)
        if('CI_sig' %in% colnames(df_all)){
            df_all <- df_all[, setdiff(colnames(df_all),'CI_sig')]
        }
        
        ## split the pvalue to 1-sided p value
        df_all$up_pval <- NA
        df_all$down_pval <- NA
        for (i in 1:nrow(df_all)){
            df_all$up_pval[i] <- splitPvalUp(df_all$pvalue[i], df_all$t.value[i])  
            df_all$down_pval[i] <- splitPvalDown(df_all$pvalue[i], df_all$t.value[i])  
        }
        
        ## adj p values:
        df_all$up_padj <- p.adjust(df_all$up_pval, adj_method)
        df_all$down_padj <- p.adjust(df_all$down_pval, adj_method)
        df_all$up_rank <- rank(df_all$up_pval, na.last = 'keep')
        df_all$down_rank <- rank(df_all$down_pval, na.last = 'keep')
        
        ## reorder columns:
        save_col <- c("geneSymbol", grep('pval$', colnames(df_all), value = T), 
          grep('padj', colnames(df_all), value = T),
          grep('rank', colnames(df_all), value = T))

        df_all <- df_all[, c(save_col, setdiff(colnames(df_all), save_col))]
        
        ## save all results
        writeTable(df=df_all, msg = paste0('# ', Sys.Date(), '\n# ', mm_dir), 
                   f_out = paste0(datadir, 'mixed_model_results.tsv'))
        print(paste0('Saved ',paste0(datadir, 'mixed_model_results.tsv')))
        
        if('df_all_expression' %in% ls()){
            print('add expression in the output result table')
            df_all_exp <- noWarnings(left_join(df_all, 
                                    df_all_expression[, c('geneSymbol', setdiff(colnames(df_all_expression), colnames(df_all)))]))
            ## save all results
            writeTable(df=df_all_exp, msg = paste0('# ', Sys.Date(), '\n# ', mm_dir), 
                       f_out = paste0(datadir, 'mixed_model_results_with_expression.tsv'))
            }

        
        
        ######
        ## save for ermineJ enrichment
        ######
        mm_df <- df_all
        for(regulation in c('up', 'down')){ ## loop 2
            print(regulation)
            
            ## save the mm files (with up and down, for ermineJ enrichment, down or up pval at column 8)
            need_col <- c("geneSymbol", "Estimate", "Std..Error", "t.value", "CI_2.5", "CI_97.5", "pvalue")
            need_col <- intersect(need_col, colnames(df_all))
            if(regulation == 'up'){
                pval_col <- 'up_pval' 
            }else{
                pval_col <- 'down_pval' 
            }
            
            df_tmp <- df_all[, c(need_col, pval_col, setdiff(colnames(df_all), c(need_col,'up_pval','down_pval')))]
            
            df_tmp <- df_tmp[order(df_tmp[, pval_col]), ] # order by up or down p value
            writeTable(df_tmp, f_out = paste0(mm_dir, '/', phase, '_', regulation,'_regulation_mixed_model.tsv'),
                       msg = paste0('# ', Sys.Date()))
            
            ## assign mixed model results based on regulation up or down
            if(regulation =='up'){
                mm_col <- 'up_pval'
            }else{
                mm_col <- 'down_pval'
            }


            if(compare_fisher){
                
            for(meta_jack in c('meta', 'jackknife')){## loop 3
                ## select the meta or jack file
                (mj_f <- list.files(jack_meta_folder,full.names = T, 
                                    pattern = paste0(phase, '.*', regulation, '.*', meta_jack, '.*.tsv')))
                
                print(mj_f)
                ## check if the file is low_exp_rm input
                if(grepl('low_exp_rm', mj_f)){
                    low_exp_rm = '_low_exp_rm'
                }else{low_exp_rm = ''}
                
                ## read meta or jack files
                df <- read.delim(mj_f, comment.char = '#')
                (index_col<- intersect(colnames(df), c('Fisher',"adj_combined_max_p")))
                rownames(df) <- df$geneSymbol
                
                ## get the results from jack or meta, and from mixed models
                genes <- data.frame(geneSymbol = intersect(df$geneSymbol, mm_df$geneSymbol))
                genes <- noWarnings(left_join(genes, df))
                genes <- noWarnings(left_join(genes, mm_df))
                genes$rank_meta <- rank(genes[, index_col], na.last = 'keep')
                genes$rank_mixed_model <- rank(genes[, mm_col], na.last = 'keep')
                
                
                ## do spearman correlation
                sp_cor_genes <- na.omit(genes[, c(index_col, mm_col)])
                (sp_cor <- cor(sp_cor_genes[, index_col], sp_cor_genes[, mm_col], method = 'spearman'))
                print(paste0(phase, '_',regulation,'_',meta_jack,low_exp_rm, " spearman cor: ", sp_cor))

                #return(list(genes, threshold))
                ## compare the top genes (after re-rank all common genes)
                (top_meta <- as.character(genes$geneSymbol[order(genes[, 'rank_meta'])][1:threshold]))
                #print(top_meta)
                (top_MM <- as.character(genes$geneSymbol[order(genes[, 'rank_mixed_model'])][1:threshold]))

                df_overlap <- data.frame(geneSymbol = intersect(top_MM, top_meta))
                df_overlap <- noWarnings(left_join(df_overlap, genes))
                df_overlap$avg_rank <- (df_overlap$rank_meta + df_overlap$rank_mixed_model)/2
                df_overlap <- df_overlap[order(df_overlap$avg_rank), ]
                writeTable(df_overlap, f_out=paste0(plotdir, phase, '_',regulation,'_',meta_jack,low_exp_rm,  '_overlap_genes.tsv'),
                           msg = paste0('# ', Sys.Date(),
                                        '\n# ', mj_f))
                
                ## plot and save
                plotXY(genes[, c('rank_meta', 'rank_mixed_model')], x= 'rank_meta', y = 'rank_mixed_model',
                       title = paste0(phase, ' ', disease, ', ',regulation, 
                                      '-regulated (compared to ', meta_jack, ')',
                                      '\nspearman correlation = ', round(sp_cor, 3),', ', nrow(genes), ' genes' ))
                ggsave(filename = paste0(plotdir, phase, '_',regulation,'_',meta_jack,'_', low_exp_rm, '.png'), width = 8, height=8)
                
                
                ## get the correlation of the top genes
                # get the top genes of mm
                (plotdir_top <- paste0(plotdir, '/top_', threshold, '/'))
                dir.create(plotdir_top, recursive = T, showWarnings = F)
                genes_top <- genes[which(genes[, 'rank_mixed_model'] <=threshold), ] %>% droplevels()
                
                #top genes: do spearman correlation
                (sp_cor <- cor(genes_top[, index_col], genes_top[, mm_col], method = 'spearman'))
                print(paste0(phase, '_',regulation,'_',meta_jack,low_exp_rm, ": top ", threshold, ": ", sp_cor))
                
                
                # top genes: plot and save
                plotXY(genes_top[, c('rank_meta', 'rank_mixed_model')], x= 'rank_meta', y = 'rank_mixed_model',
                       title = paste0(phase, ' ', disease, ', ',regulation, 
                                      '-regulated (compared to ', meta_jack, ')',
                                      '\nspearman correlation = ', round(sp_cor, 3),', top ', nrow(genes_top), ' genes' ))
                ggsave(filename = paste0(plotdir_top, phase, '_',regulation,'_',meta_jack,'_', low_exp_rm, '.png'), width = 8, height=8)
                
                ## get the stat
                df_t <- data.frame(gene_count = nrow(genes),
                                   disease_phase = phase,
                                   regulation = regulation,
                                   compare_type = meta_jack,
                                   spearman_corr = sp_cor, 
                                   low_exp_rm = low_exp_rm,
                                   top_genes_overlap = paste0(nrow(df_overlap), '/', threshold),
                                   input_meta_jack = mj_f,
                                   input_mixed_model_result = expr_Robject)
                
                
                stat_df <- rbind(stat_df, df_t)
                
            }  ## meta loop 3
        }
        } 
    }
    return(stat_df)
}


#---------------------------------------------------------------------------#
# Main function
#---------------------------------------------------------------------------#

mainCompareMM <- function(jack_meta_folder_ls, mm_dir,phase_ls = c('early', 'late'), threshold =200, adj_method = 'BH',
                          compare_fisher =T){
    print(mm_dir)
    msg <- paste0('# ', Sys.Date(),
                  '\n# input mixed model result folder: ', mm_dir )
    stat_df_all <- NULL
    for(jack_meta_folder in jack_meta_folder_ls){
        stat_df <- compareMM(jack_meta_folder, mm_dir, phase_ls=phase_ls, threshold = threshold, 
                             adj_method = adj_method, compare_fisher = compare_fisher)
        stat_df_all <- rbind(stat_df_all, stat_df)
        msg <- paste0(msg, '\n# input meta-analysis folder: ', jack_meta_folder)
    }
    
    f_save <- paste0(mm_dir,'/mm_comparison_to_meta_analysis.tsv')
    writeTable(stat_df_all, f_out = f_save, msg = msg)
    print(paste0('STAT OUTPUT: ', f_save ))
}

