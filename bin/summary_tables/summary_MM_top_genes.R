#' 2016-08-26
#' updated 2016-09-06, 2016-10-24
#' summary, run this script
#' to get to genes for each ND, and top pathways of each ND
#' summarizeTopGenes()  ##mixed model
#' summarizeTopPathways()
#' summarizeTopJackGenes()  ## must run after summarizeTopGenes()
#' 

#'@example
# 
# source('summary_tables/summary_MM_top_genes.R')
# f_out <- paste0('../../results/ND_results/ND_summary/doc_tables/top_genes/')
# 
# threshold =10
# disease_ls <- c('AD', 'HD', 'PD')
# summarizeTopGenes(disease_ls, phase ='early', f_out, threshold =10)
# summarizeTopGenes(disease_ls, phase = 'late', f_out, threshold =10)
# 
# source('summary_tables/summary_MM_top_genes.R')
# disease_ls <- c('AD', 'HD', 'PD')
# 
# f_out <- paste0('../../results/ND_results/ND_summary/doc_tables/pathways/')
# summarizeTopPathways(disease_ls, f_out,
#                      phase_ls=c('early', 'late'),
#                      regulation_ls = c('up', 'down'),
#                      threshold =5)

library('dplyr')
library(HelperFunctions)





#########################################################
## mixed model summary of top genes in all diseases, and save an rdata of all the results
#########################################################

summarizeTopGenes <- function(disease_ls, phase, f_out, f_rdata_out=f_out, threshold =10, top_genes =T, rdata_pre="",
                              mm_folder, mgi_info){
    #' for mixed model only
    #' f_dir: full dir to 'mixed_model_results.tsv', if NA, then auto search in the disease dir mixed model results
    #' top_genes: write a table with the top genes of all disease and phase
    #' f_out: outdir for topgenes results
    #' f_rdata_out outdir for mixed model results
    
    ## get gene names
    df_gene_all <- read.delim(mgi_info,comment.char = '#')
    df_gene <- df_gene_all[, c('geneSymbol', 'Name')]%>%droplevels()
    
    
    dir.create(f_out, recursive = T, showWarnings = F)
    dir.create(f_rdata_out, recursive = T, showWarnings = F)
    print(phase)  
    fo_all <- paste0(f_out, phase,'_top_', threshold, '_genes_', Sys.Date(), '.tsv')
    if(file.exists(fo_all)){
        file.remove(fo_all)
    }
    
    
    for(i in 1:length(disease_ls)){
        disease <- disease_ls[i]
        print(disease)
        disease_dir <- paste0(home_dir, '/',disease, '_mouse_model_project/')
        result_dir <- paste0(disease_dir, '/mixed_model/', mm_folder,'/', phase,'/')  ## mixed model only
        print(result_dir)
        
        # if(is.na(f_dir)){ ## if the mixed_model_results is not provided, search the given disease dir and phase
        ## get the most result result folder and sub folder with the key word
        (f_dir <- max(grep('mixed_model_results.tsv', list.files(result_dir, full.names = T), value = T)))
        #}
        if(length(f_dir) !=1){
            stop(paste0(result_dir, 'doesnt contain mixed_model_results.tsv'))
        }
        print(f_dir)
        
        
        df_input <- read.delim(f_dir, comment.char = '#')
        ## add gene name
        df_input <- noWarnings(left_join(df_input,df_gene))
        ## add pvalue adj (adj for raw, not up or down pval)
        df_input$pvalue_adj <- p.adjust(df_input$pvalue, method='BH')
        
        ## save a r data with mm result dataframe
        df <- noWarnings(left_join(df_input, df_gene_all))
        
        need_col <- c("geneSymbol", 'Name', 'pvalue', "pvalue_adj", "up_rank", "down_rank",
                      "Feature.Type", "Chr","Strand","Start","End","Human_Disease")
        
        (col_order <- c(need_col, setdiff(colnames(df), need_col)))
        df <- df[, col_order]
        assign(paste0(disease, "_", phase, "_mm_variable"), value = df)
        
        
        ##### get the top genes
        if(top_genes){
            print('RETURN TOP GENES')
            #### get the top genes
            
            ## upregulated
            top_col <- c("geneSymbol",'Name',"pvalue", "pvalue_adj","Estimate","Std..Error")
            ## get top genes
            df_up <- df_input[which(df_input$up_rank <= threshold), top_col]%>%droplevels()
            
            df_up <- orderCol(df_up,cols = 'pvalue')
            df_up$regulation= ''
            df_up$regulation[1]= paste0(phase, ' phase: up-regulated')
            colnames(df_up)[3:4] <- c('p', 'adj_p')
            
            ## down regulated
            ## get top genes
            df_down <- df_input[which(df_input$down_rank <= threshold), top_col]%>%droplevels()
            
            df_down <- orderCol(df_down,cols = 'pvalue')
            df_down$regulation= ''
            df_down$regulation[1]= paste0(phase, ' phase: down-regulated')
            colnames(df_down)[3:4] <- c('p', 'adj_p')
            
            df <- rbind(df_up, df_down)
            df$Disease <- ''
            df$Disease[1] <- disease
            
            
            ## round the numbers
            for(col_r in c("adj_p","Estimate","Std..Error")){
                df[, col_r] <- round(df[, col_r], digits = 4)
            }
            
            for(col_r in c("adj_p","p")){
                df[, col_r] <- format(df[, col_r], scientific=T, digits = 3)
            }
            
            ## update column names and order
            order_col <- c('Disease', 'regulation',"geneSymbol",'Name', "p", "adj_p","Estimate","Std..Error")
            df <- df[, order_col]
            new_col <- c('Disease', 'Regulation',"Gene",'Name', "P-value", "Adj. p","Slope estimate","Std. error")
            colnames(df) <- new_col
            ## write to table
            fo <- paste0(f_out, disease, '_', phase, '_top_', threshold, '_genes_', Sys.Date(), '.tsv')
            writeTable(df, fo_all, file_append =T)
            print(paste0('# ', Sys.Date(), ': ', f_dir))
            print(fo)
        }
    }## loop end for each disease
    
    save_ls <- grep('mm_variable', ls(), value = T)
    (cmd <- paste0('save(', paste0(save_ls, collapse = ", "), ", file= '", f_rdata_out,rdata_pre, phase, "_mm_results.Rdata')"))
    print(cmd)
    eval(parse(text = cmd))
}


#########################################################
## mixed model summary of top pathways
#########################################################
summarizeTopPathways <- function(disease_ls, f_out,model_ls,
                                 phase_ls=c('early', 'late'), 
                                 regulation_ls = c('up', 'down'),  
                                 threshold =5,
                                 mm_jack_folder){
    # e.g. mm_jack_folder = 'random_intercept_include_NA_low_exp_rm' # which ermineJ subfolder to grep files from 
    
    dir.create(f_out, recursive = T, showWarnings = F)
    

    
    for (model in model_ls){
        (fo_all <- paste0(f_out, model, '_top_', threshold, '_pathways_', Sys.Date(), '.tsv'))
        if(file.exists(fo_all)){
            file.remove(fo_all)
        }
        for(i in 1:length(disease_ls)){
            disease <- disease_ls[i]
            print(disease)
            disease_dir <- paste0(home_dir, disease, '_mouse_model_project/')
            result_dir <- paste0(disease_dir, 'ermineJ/')
            ## get the most recents
            (result_dir <- max(grep(paste0(model, '\\/', mm_jack_folder,'\\/.*analysis$'), 
                                list.dirs(result_dir), value = T)))
            print(result_dir)
            
            ## get the most recent result folder and sub folder with the key word
            (f_dir_ls <- grep('mixed_model.tsv$', list.files(result_dir, full.names = T), value = T))
            if(length(f_dir_ls) !=4){
                stop(paste0(result_dir, 'doesnt contain 4 mixed_model.tsv files'))
            }
            for(phase in phase_ls){
                print(phase)  
                df_all <- NULL
                
                for(regulation in regulation_ls){ ## loop reg
                    f_dir <- grep(paste0(phase, '.*', regulation), f_dir_ls,value = T)
                    print(f_dir)
                    df_input <- read.delim(f_dir, comment.char = '#')
                    
                    colnames(df_input)
                    
                    need_col <- c('ID','Name','Pval', "CorrectedPvalue",'top_genes')
                    
                    df <- df_input[1:threshold, need_col]%>%droplevels()
                    df$Disease <- ''
                    df$Disease[1] <- disease
                    df$regulation= ''
                    df$regulation[1]= paste0(phase, ' phase: ', regulation, '-regulated')
                    
                    ## reorder columns
                    df <- df[, c('Disease', 'regulation', need_col)]
                    colnames(df) <- c('Disease', 'Regulation', 'GO term','Description','P-value','Adj. p', 'Top genes')
                    
                    ## round numbers
                    for(col_r in c('P-value','Adj. p')){
                        df[, col_r] <- round(df[, col_r], digits = 4)
                    }
                    
                    if(is.null(df_all)){
                        df_all <- df
                        
                    }else{
                        df_all <- rbind(df_all,df)
                    }
                } ## loop reg
                
                
                ## write all results per phase
                noWarnings(writeTable(df_all, fo_all, file_append = T))
                
            }
        } ## end loop for disease
    }
 
    print(fo_all)
}





#########################################################
## jackkife mixed model summary of top genes in all diseases, and save an rdata of all the results
## must run aftersummarizeTopGenes()
#########################################################

summarizeTopJackGenes <- function(disease_ls, phase, f_out, f_rdata_out=f_out, threshold =10, top_genes =T,
                                  mm_jack_folder){
    #' top_genes: write a table with the top genes of all disease and phase
    #' f_out: outdir for topgenes results
    #' f_rdata_out outdir for mixed model results rdata
    #' 
    #
    dir.create(f_out, recursive = T, showWarnings = F)
    dir.create(f_rdata_out, recursive = T, showWarnings = F)
    print(phase)  
    (fo_all <- paste0(f_out,'jackknife_', phase,'_top_', threshold, '_genes_', Sys.Date(), '.tsv'))
    if(file.exists(fo_all)){
        file.remove(fo_all)
    }
    
    
    ## load the precalculated mixed model for all diseases of the phase
    (rdata <- grep(paste0(phase, "_mm_results.Rdata"), list.files(f_rdata_out, full.names = T) , value = T))
    (rdata <- grep('jackknife', rdata, value = T, invert = T)) ## grep not the jackknife result
    load(rdata)
    print(paste0('INPUT RDATA: ', rdata))
    
    
    for(i in 1:length(disease_ls)){
        disease <- disease_ls[i]
        print(disease)
        (disease_dir <- paste0(home_dir, disease, '_mouse_model_project/'))
        ## get the variable with disease and phase
        (r_obj <- grep(paste0(disease, '_', phase,'_mm_variable'), ls(), value = T))
        (cmd <- paste0('df <- ', r_obj))
        eval(parse(text = cmd))
        
        ## get the jackknife results of the phase and regulation
        disease_dir <- paste0(home_dir,'/' ,disease, '_mouse_model_project/')
        (jack_dir <- paste0(disease_dir, 'mixed_model_jackknife/',mm_jack_folder,'/'))
        print(jack_dir)
        
        ### for each regulation:
        df_all <- NULL
        for (regulation in c('up','down')){
            (f_jack <- grep(paste0(phase,'_', regulation, '.*mixed_model_results.tsv' ),list.files(jack_dir, full.names = T), value = T))
            
            df_jack <- read.delim(f_jack, comment.char = '#')
            need_col <- which(colnames(df_jack) %in% c('max_jackknife_rank', 'final_jackknife_rank'))
            colnames(df_jack)[need_col] <- paste0(regulation, '_',colnames(df_jack)[need_col])
            df_jack <- df_jack[, c(1, need_col)]
            
            ## join
            if(is.null(df_all)){
                df_all <- noWarnings(left_join(df_jack, df))
            }else{
                df_all <- noWarnings(left_join(df_jack, df_all))
            }
        }
        
        ## reassign the variable value
        (cmd <- paste0(r_obj, '<- df_all'))
        eval(parse(text = cmd))
        
        ## get the top genes for up and down
        if(top_genes){
            print('RETURN TOP GENES')
            #### get the top genes
            
            ## upregulated
            top_col <- c("geneSymbol",'Name',"pvalue", "pvalue_adj","Estimate","Std..Error")
            ## get top genes
            df_up <- orderCol(df_all,cols = c('up_final_jackknife_rank', 'up_max_jackknife_rank') )
            df_up <- df_up[1:threshold, c('up_final_jackknife_rank', top_col)]%>%droplevels()
            df_up$regulation= ''
            df_up$regulation[1]= paste0(phase, ' phase: up-regulated')
            colnames(df_up)[c(1, 4:5)] <- c('rank', 'p', 'adj_p')
            
            ## down regulated
            ## get top genes
            df_down <- orderCol(df_all,cols = c('down_final_jackknife_rank', 'down_max_jackknife_rank') )
            df_down <- df_all[1:threshold, c('down_final_jackknife_rank',top_col)]%>%droplevels()
            df_down$regulation= ''
            df_down$regulation[1]= paste0(phase, ' phase: down-regulated')
            colnames(df_down)[c(1, 4:5)] <- c('rank','p', 'adj_p')
            
            df <- rbind(df_up, df_down)
            df$Disease <- ''
            df$Disease[1] <- disease
            
            
            ## round the numbers
            for(col_r in c("adj_p","Estimate","Std..Error")){
                df[, col_r] <- round(df[, col_r], digits = 4)
            }
            
            for(col_r in c("p")){
                df[, col_r] <- format(df[, col_r], scientific=T, digits = 3)
            }
            
            ## update column names and order
            order_col <- c('Disease', 'regulation','rank',"geneSymbol",'Name', "p", "adj_p","Estimate","Std..Error")
            df <- df[, order_col]
            new_col <- c('Disease', 'Regulation','Rank',"Gene",'Name', "P-value", "Adj. p","Slope estimate","Std. error")
            colnames(df) <- new_col
            ## write to table
            writeTable(df, fo_all, file_append =T)  ## write out for all disease
        }## for top genes
    } ## for each disease
    
    
    print(fo_all)
    
    save_ls <- grep('mm_variable', ls(), value = T)
    (cmd <- paste0('save(', paste0(save_ls, collapse = ", "), ", file= '", f_rdata_out,'jackknife_', phase, "_mm_results.Rdata')"))
    print(cmd)
    eval(parse(text = cmd))
}

