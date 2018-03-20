## get the jackknife Mm model results and combine all runs in a file (also for ermineJ)
#2016-08-25

#####
library(dplyr)
library(HelperFunctions)


#### helper functions
grepFiles <- function(mm_jack_dir_p, regulation, phase){
    ## get files from all the runs of the mixed model result, and corresponding subset that's been removed
    f_ls <- NULL
    rm_msg <- NULL
    for(f in mm_jack_dir_p){  # for each jackfile
        (f_ls <- c(f_ls, grep(paste0(phase, '.*', regulation), list.files(f, full.names = T), value =T)))
        (f_rm <- grep('rm_sample.txt$', list.files(f, recursive = T, full.names = T), value = T))
        # print(f_rm)
        if(length(f_rm) >0){
            (rm_msg <- c(rm_msg, noWarnings(readLines(f_rm))[2]))
        }
    }
    returnlist=list(f_ls, rm_msg)
    return(returnlist)
}

########main function
compareJackMM <- function(mm_jack_dir,mm_dir, phase, regulation_ls = c('up', 'down'),
                          return_df=F,df_out_dir =NULL, notes = '',file_pre = ''){
    ## save the ranking result from each run and the result from the full MM (no samples removed)
    ## get the jackknife rank (average rank of all the jackknife runs)
    #' return_df and df_out_dir: output a better version of the table for publication
    #' notes: addition info to add to table comments
    #' file_pre: add to table out file name 
    
    ## get the phase and regulation file
    mm_jack_dir_p <- paste0(mm_jack_dir, '/', phase,'/')
    (mm_jack_dir_p <- paste0(grep('done|run', list.dirs(mm_jack_dir_p, recursive = F), value = T), '/'))
    
    ###########
    ## for each regulation, return a df with all rankings from each run of the regulation,
    ###########
    for(regulation in regulation_ls){
        print(regulation)
        all_cols <- c("up_pval", "up_padj", "down_pval", "down_padj","up_rank","down_rank")
        df_all <- NULL
        msg_all <- NULL
        x <- grepFiles(mm_jack_dir_p, regulation, phase)
        (f_ls <- x[[1]])
        (rm_msg <- x[[2]])
        for(i in 1:length(f_ls)){
            # print(regulation)
            msg <- paste0('run_', i, ': ', rm_msg[i])
            print(msg)
            msg_all <- c(msg_all, msg)
            
            (mmf=f_ls[i])
            # print(paste0('INPUT: ', mmf))
            df <- read.delim(mmf, comment.char = '#')
            rk_col <- grep(regulation, all_cols, value = T)
            # print(paste0('rank col is ', rk_col))
            (need_col <- c('geneSymbol','pvalue', rk_col))
            df <- df[, need_col] %>%droplevels()
            colnames(df) <- c('geneSymbol',paste0('run_', i, '_', c('pvalue', rk_col)))
            row.names(df) <- df$geneSymbol
            if(is.null(df_all)){
                df_all <- df
            }else{
                df_all <- noWarnings(left_join(df_all, df))
            }
        } # end for each f_ls
        ## calculate avg rank, max rank
        rank_c <- grep('rank', colnames(df_all), value = T)
        df_all$avg_jackknife_rank <- apply(df_all[, rank_c], 1, function(x) sum(x)/length(rank_c))
        df_all$max_jackknife_rank <- apply(df_all[, rank_c], 1, max)
        df_all$avg_rank_percent <- df_all$avg_jackknife_rank/nrow(df_all)
        
        ### load the orginial mixed model with no sample ommited
        (full_mm <- grep(paste0(phase, '.*', regulation), list.files(mm_dir, full.names = T,recursive = T), value =T))
        print(paste0('Input full mm: ', full_mm))
        df_mm <- read.delim(full_mm, comment.char = '#')
        # df_mm <- df_mm[, need_col] %>%droplevels()
        # colnames(df_mm) <- c('geneSymbol',paste0('all_samples_', c('pvalue', rk_col)))
        df_mm$FDR <-  p.adjust(df_mm$pvalue, method = 'BH')
        colnames(df_mm)[2: ncol(df_mm)] <- paste0('all_samples_', colnames(df_mm)[2: ncol(df_mm)] )
        df_all <- noWarnings(left_join(df_all, df_mm))
        
        ## make the avg_rank_percent as the 2 col for ermineJ
        all_col <- c('geneSymbol', 'avg_rank_percent', setdiff(colnames(df_all), c('geneSymbol', 'avg_rank_percent')))
        df_all <- df_all[, all_col]
        
        ### get msg
        msg_final <- paste0('# ', Sys.Date(), 
                            '# ', notes,
                            '\n# INPUT FILES: \n#   ',
                            paste0(f_ls, collapse = '\n#   '),
                            '\n# INPUT FULL MIXED MODEL: ', full_mm,
                            '\n# ', paste0(msg_all, collapse = '\n# '),
                            '\n# final jackknife rank is the ranking of the average_jackknife_rank')
        
        ### reorder columns and save 
        df_all$final_jackknife_rank <- rank(df_all$avg_jackknife_rank)
        first_col <- grep('geneSymbol|jackknife|all_samples.*rank|avg_rank_percent',colnames(df_all), value = T)
        df_all <- df_all[,c(first_col, setdiff(colnames(df_all), first_col))]
        df_all <- orderCol(df_all, 'final_jackknife_rank')
        
        ## get genenames
        df_all <- noWarnings(left_join(df_all, genenames))
        
        (f_out <- paste0(mm_jack_dir, 'jackknife_', phase, '_', regulation, '_regulation_mixed_model_results.tsv'))
        writeTable(df=df_all, f_out=f_out, msg=msg_final)
        print(paste0('RESULT OUT: ', f_out))
        
        if(return_df){
            
            if(is.null(df_out_dir)){
                (df_out_dir <- paste0(mm_jack_dir, '/',disease,'_rank_tables/'))
            }
            dir.create(df_out_dir, showWarnings = F, recursive = T)
            
            df <- df_all
            colnames(df) <- gsub('all_samples_', '', colnames(df))
            
            need_col <- c( "final_jackknife_rank","geneSymbol", 'Name', "pvalue","FDR",
                           "Estimate",
                           "Std..Error","t.value")
            
            
            (random_effects <- grep('_intercept',colnames(df), value = T))
            (jack <- grep('run',colnames(df), value = T))
            ## add the cell type FE columns
            cell_type <- c('Oligo', 'Astrocyte', "Microglia","DentateGranule",'GabaSSTReln', 'StriatumCholin',
                           'Pyramidal', 'ForebrainCholin', 'Spiny')
            (cell_col_index <- grep(paste0(cell_type, collapse = '|'),colnames(df)))
            if(length(cell_col_index) >0 ){
                cell_col <- c('Oligodendrocytes','Astrocytes',"Microglia", "Dentate_granule_cells" ,'GABAergic_cells','Cholinergic_neurons', 
                              'Pyramidal_cells','Cholinergic_neurons', 'Medium_spiny_neurons')
                names(cell_col) <- cell_type
                colnames(df)[cell_col_index] <- mapValue(colnames(df)[cell_col_index], cell_col)
                colnames(df)[cell_col_index] <- paste0(colnames(df)[cell_col_index], '_FE')
                ## update column names
                df <- df[,c(need_col, random_effects, colnames(df)[cell_col_index],jack)]
            }else{
                ## update column names, no cell correction
                df <- df[,c(need_col, random_effects, jack)]
            }
            
            

            colnames(df) <- gsub('_intercept', '_RE', colnames(df))
            colnames(df) <- gsub('run_', 'jackknife_run', colnames(df))
            colnames(df) <- gsub('Estimate', 'FE_beta_estimate', colnames(df))
            
            
            msg_output <- paste0('#', disease, ': ', phase, '\n# ',notes, 
                                 '\n# ', paste0(msg_all, collapse = '\n# '),
                                 '\n# Final jackknife rank is the ranking of the average ranks of all the jackknife runs',
                                 '\n# FE: fixed effects estimates; RE: random effects estimates',
                                 '\n# pvalue and all the estimates are calculated from all the samples of the disease phase (not from jackknife procedure estimates)')
            
            
            (f_out <- paste0(df_out_dir, disease, '_', phase, '_', regulation, '_jackknife_mixed_model_results', file_pre, '.tsv'))
            writeTable(df=df, f_out=f_out, msg=msg_output)
            print(paste0('Table OUT: ', f_out))
            
            
            # return(df_all)
        }
        
    } #end loop for each regulation
    
}
