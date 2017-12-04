#' 2016-08-04
#' get the top 500 genes from meta, jackknife or mixed model results
#' and get the top [threshold], pathways of the correponding erminej file
#' look at the ranking to the geneMembers of these pathways
#' write the dataframe to the ermineJ folder
# see wrapper
# (jack_meta_folder <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/'))
# (ermineJ_folder <- max(list.dirs(paste0(disease_dir, 'ermineJ/mixed_model/', model,model_keyword,'/'), recursive = F)))
# source('result_explore/top_genes_for_ermineJ.R')



#---------------------------------------------------------------------------#
# PART 1: SET SYS DEFAULTS
#---------------------------------------------------------------------------#
# 
source("helper_functions.R")

## only look at the top genes
if(!('top_threshold' %in% ls())){
  top_threshold = 200
}
#---------------------------------------------------------------------------#
# PART 3: LOOPs
#---------------------------------------------------------------------------#

for(keyword in keyword_ls){ ## loop1 keyword (meta, jackknife, mixed_model)
    print(paste0('start ', keyword))
    ## get the meta or jack files or mm up and down regulations
    (f_gene_list_all <- list.files(jack_meta_folder, full.names=T,
                                   pattern=paste0('_regulation_',keyword, '.*.tsv')))
    #print(f_gene_list_all)
    # get the top up and down gene list
    
    (f_ermineJ <- list.files(ermineJ_folder, full.names = T,pattern=paste0('_regulation_',keyword, '.*.tsv')))
    
    ## for combined disease files (must specify which model and model keyword for mixed model combined results:
    if('mixed_model_only' %in% ls()){
        if(mixed_model_only){
            (f_gene_list_all <-grep(paste0(model, model_keyword), f_gene_list_all, value = T))
            (f_ermineJ <- grep(paste0(model, model_keyword), f_ermineJ, value = T))
        }
    }
    
    
    # grep not table
    (f_ermineJ <- grep('table.tsv', f_ermineJ, invert = T, value = T))
    
    for (phase in phase_ls){##loop 2 phase
        print(phase)
        for(regulation in regulation_ls){ ## loop 3 regulation
            print(regulation)
            (f_gene_list <-grep(regulation, f_gene_list_all, value = T))  ## get the right regulation
            (f_gene_list <-grep(phase, f_gene_list, value = T)) ## get the right disease phase
            ## get the corresponding ermineJ files
            (f_ej_list <-grep(regulation, f_ermineJ, value = T))  ## get the right regulation
            (f_ej_list <-grep(phase, f_ej_list, value = T)) ## get the right disease phase
            #***************************
            #' loop for each disease phase and correcponding studies
            ## f_gene_list: a meta or jack file
            ## threshold: top # of genes in the f_gene_list
            ## info_df: a df with File(the limma object files), Phase, Extra_phase
            #***************************
            print(f_gene_list)
            print(f_ej_list)
            
            ## get the genes from jack or meta or mixed model file
            df <- read.delim(f_gene_list,comment.char = "#")
            
            ## reorder by the score column, and give rank, rank_avg_percent is the combined disease rank
            (score_col <- intersect(c('adj_combined_max_p','Fisher', 'up_pval', 'down_pval','rank_avg_percent','avg_rank_percent'), colnames(df)))
            df <- df[order(df[, score_col]), ]
            df$rank <- rank(df[, score_col])
            df$output <- paste0(df$geneSymbol, '(', df$rank, ')')
            
            ## only look at the top genes 
            df <- df[1:top_threshold, ]

            
            ## read ermineJ file, (skip the first column)make a tmp file and then read the tmp file
            df_line <- noWarnings(readLines(f_ej_list))
            (index <- max(grep('#!', df_line[1:100])-1))
            df_line <- df_line[-(1:index)]
            df_line <- gsub('#!\t','', df_line)
            df_line <- gsub('!\t', '', df_line)
            df_line <- gsub('GeneMembers', 'GeneMembers\t', df_line)
            
            (f_out <- paste0(ermineJ_folder,'/', phase, '_', regulation, '_regulation_',keyword, '_table.tsv'))
            sink(file=f_out, type = 'output')
            writeLines(df_line)
            sink()
            
            
            ermine_df <- read.delim(f_out, stringsAsFactors = F) 
            ## reorder by corrected p value
            ermine_df <- orderCol(ermine_df, c("CorrectedPvalue","Pval"))


            (index <- which(colnames(ermine_df) !='X'))
            ermine_df <- ermine_df[,index ]
            ## ad rewrite the table
            writeTable(ermine_df, f_out = f_out)
            
            
            ## get the top erminej rows, max(top pathways or top pathways <=FDR)
            if(!is.na(threshold)){
                fdr <- which(ermine_df$CorrectedPvalue <= fdr_threshold)
                if(length(fdr) > threshold){
                    threshold_index <-fdr 
                }else{
                    threshold_index <-1:threshold
                }
                
                print(paste0('only top ', length(threshold_index), ' pathways are analyized'))
                ermine_df <- ermine_df[threshold_index, ]
            }else{
                print('keep all pathways')
            }

            ermine_df$top_genes <- ""
            
            df_gene_list <- NULL
            for (i in 1:nrow(ermine_df)){
                (members <- data.frame(geneSymbol =unlist(strsplit(ermine_df[i, 'GeneMembers'], split = '\\|'))))
                df_member <- noWarnings(left_join(members, df[, c('geneSymbol', 'output', 'rank')]))
                ## remove NA and order by rank
                df_member <-na.omit(df_member)
                df_member <- orderCol(df_member, 'rank')
                ermine_df$top_genes[i] <- paste0(df_member$output, collapse = ', ')
                df_member$name <- paste0(ermine_df$Name[i], '- FDR:',ermine_df$CorrectedPvalue[i], '[', ermine_df$NumGenes[i], '](', ermine_df$ID[i] ,')')
                df_member$NumGenes <- ermine_df$NumGenes[i]
                df_member$CorrectedPvalue <- ermine_df$CorrectedPvalue[i]
                if(is.null(df_gene_list)){
                    df_gene_list <- df_member
                }else{
                    df_gene_list <- rbind(df_gene_list,df_member)
                }
            }
            
            
            ## save a table with the top genes involved and which pathways are involved 
            df_gene_list <- orderCol(df_gene_list, c('geneSymbol', 'NumGenes'))
            df_tmp1 <- as.data.frame(table(df_gene_list$geneSymbol))
            colnames(df_tmp1) <- c('geneSymbol', 'Count')
            df_tmp2 <- aggregate(name~geneSymbol, data = df_gene_list, function(x) paste0(x, collapse = '; '))
            df_tmp2 <- noWarnings(left_join(df_tmp1, df_tmp2))
            df_tmp2 <- noWarnings(left_join(rmDup(df_gene_list[,c('geneSymbol', 'rank')]), df_tmp2))
            df_tmp2 <- orderCol(df_tmp2, 'rank')
            
            ## which ones are under FDR threshold
            df_gene_list2 <- df_gene_list[which(df_gene_list$CorrectedPvalue <= fdr_threshold), ]%>%droplevels()
            if(nrow(df_gene_list2) >0){
                df_tmp3 <- aggregate(name~geneSymbol, data = df_gene_list2, function(x) paste0(x, collapse = '; '))
                colnames(df_tmp3)[2] <- 'sig_paths'
                
                df_sig_gene <- noWarnings(left_join(df_tmp2, df_tmp3))
            }else{
                df_sig_gene <- df_tmp2
                df_sig_gene$sig_paths <- NA
            }

            
            ## write
            f_out <- paste0(ermineJ_folder, '/analysis/')
            dir.create(f_out, showWarnings = F)
            
            if(all('mixed_model_only' %in% ls(), mixed_model_only)){
                (f_o <- paste0(f_out, model, model_keyword, '_', phase, '_', regulation, '_regulation_top_genes_',keyword, '.tsv'))
            }else{
                (f_o <- paste0(f_out, phase, '_', regulation, '_regulation_top_genes_',keyword, '.tsv'))
            }
            
            
            msg <- paste0('# ', Sys.Date(), '\n# Input score file: ', f_gene_list, 
                          '\n# Input erminej file: ', f_ej_list,
                          '\n# only top ', length(threshold_index), ' pathways',
                          '\n# only top ', top_threshold, ' genes from the score file are shown ranking in ()',
                          '\n# FDR: corrected p value, []: gene member size of the pathway')
            writeTable(df = df_sig_gene, f_out = f_o, msg = msg)
            print(f_out)
            
            ## reorder columns
            tmp <- c("Name","ID","CorrectedPvalue","top_genes","NumGenes","MFPvalue")
            ermine_df <- ermine_df[,c(tmp, setdiff(colnames(ermine_df), tmp))]
            
            
            

            
            
            
            if(all('mixed_model_only' %in% ls(), mixed_model_only)){
                (f_out <- paste0(f_out, model, model_keyword, '_', phase, '_', regulation, '_regulation_',keyword, '.tsv'))
            }else{
                (f_out <- paste0(f_out, phase, '_', regulation, '_regulation_',keyword, '.tsv'))
            }
            

            msg <- paste0('# ', Sys.Date(), '\n# Input score file: ', f_gene_list, 
                          '\n# Input erminej file: ', f_ej_list,
                          '\n# only top ', threshold, ' pathways',
                          '\n# only top ', top_threshold, ' genes from the score file are shown ranking in ()')
            writeTable(df = ermine_df, f_out = f_out, msg = msg)
            print(f_out)
        } # loop3
    }## loop2
}