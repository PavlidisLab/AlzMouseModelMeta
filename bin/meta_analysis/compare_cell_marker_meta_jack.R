# 2016-05-19
# updated 2016-05-24
#' this script is to look at the cell type markers and use Fisher's exact test 
#' to see if there's is over representation of cell type markers on the top genes
#' of meta and jackknife files
#' usage example:
#' 
# rm(list=setdiff(ls(),'home_dir'))
#  
# source('config/config_wrappers.R')
# source('meta_analysis/compare_cell_marker_meta_jack.R')
# 
# out_dir <-paste0(disease_dir, '/meta_analysis/cell_marker_comparison/')
# dir.create(out_dir,showWarnings=F, recursive=T)
# 
# (f_marker_genes <- cell_pop_marker_genes_config)
# (jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))
# ## meta files
# (f_mj_ls <- list.files(jack_meta_folder, pattern = 'meta.*tsv', full.names = T))
# threshold_method_ls <- c('top', 'percent', 'q_value')
# threshold_ls <- c(100, 0.1, 0.1)
# mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir)
# ## jack files
# (f_mj_ls <- list.files(jack_meta_folder, pattern = 'jackknife.*tsv', full.names = T))
# threshold_method_ls <- c('top', 'percent', 'q_value')
# threshold_ls <- c(100, 0.1, 0.2)
# mainCellTypeTest(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir)

source('helper_functions.R')

#---------------------------------------------------------------------------#
# Functions
#---------------------------------------------------------------------------#

#' input a folder with celltype and their marker genes (from Ogan)
#' output a df with geneSymbol(marker genes) and the cell type
markerGenes <- function(f_marker_genes){
    celltypes <- list.files(f_marker_genes, recursive = T)
    f_markers <- list.files(f_marker_genes, recursive = T, full.names = T)
    df_all <- NULL
    for(i in 1:length(celltypes)){
        cell <- celltypes[i]
        f_m <- f_markers[i]
        df <- read.delim(f_m, comment.char = '#', header = F)
        colnames(df) <- 'geneSymbol'
        df$cell_marker <- cell
        df_all <- rbind(df_all, df)
    }
    return(df_all)
}


#'For each input jack or meta, do the fisher's test for each cell type (for 1 threshold method)
#'NA in meta and jackknife are not removed
#'@param f_mj a meta or jack file dir
#'@param df_marker: a df with geneSymbol and cell_marker from markerGenes()
#'@param f_test_out: outdir for the test results
#'@param save_f, f_out: to save a jack or meta file with cell markers annotated, and the save dir
#'@param rm_na: rm the records with 'Fisher_BH_adjusted', 'adj_combined_max_p' = NA
#'@param f_result_out: results of the fisher test
#'
#'
cellTypeFisherTest <- function(f_mj, df_marker, threshold, f_test_out,
                               threshold_method = c('top', 'percent', 'q_value'), 
                               save_f =F, f_out = "", rm_na = T, f_result_out = 'celltype_fisher_test_results.tsv'){
    df_all <-NULL ## store all test results
    ## read meta or jack files
    df_mj <- read.delim(f_mj, comment.char = '#')
    (rank_col <- intersect(c('Fisher_BH_adjusted', 'adj_combined_max_p'), colnames(df_mj)))
    if(rm_na){
        index <- which(is.na(df_mj[, rank_col]))
        if(length(index >0)){
            df_mj <- df_mj[-index, ]
        }
        msg_na <- length(index)
    }else{msg_na <- 'not removed'}
    

    
    ## get the columns
    df_mj<- df_mj[, setdiff(colnames(df_mj), c("GOTerms","GemmaIDs","NCBIids", "platform"))] %>%droplevels()
    df <- noWarnings(left_join(df_mj, df_marker))
    df$cell_marker <- as.factor(df$cell_marker)
    df_other <- noWarnings(left_join( df_marker,df_mj))

    ## sort by the rank col
    df <- df[order(df[, "Rank"]), ]

    if(save_f){
        ## if save the jack or meta file with the cell type annotated
        writeTable(df, f_out = f_out)
    }
    
    #return()
    ## set the threshold (count of numbers under threshold)
    threshold_t <- switch(threshold_method,
                          top = threshold,
                          percent = round(threshold * nrow(df)),
                          q_value = length(which(df[, rank_col] <= threshold))
    )

    if(threshold_t ==0){
        print('No genes under the threshold')
        df_result <- data.frame(input_file = f_mj, 
                                NA_genes_removed = msg_na,
                                total_genes = nrow(df),
                                threshold_method = threshold_method, threshold= threshold,
                         genes_under_threshold = threshold_t, 
                         cell_type = NA,
                         marker_under_threshold = NA,
                         non_marker_under_threshold = NA,
                         marker_over_threshold = NA,
                         non_marker_over_threshold = NA,
                         Fisher_exact_test = NA,
                         p_value = NA,
                         odds_ratio =  NA)
        df_all <- rbind(df_all, df_result)
    }else{
        ## do the fisher's test for each cell type markers
        for (cell_t in levels(df$cell_marker)){

            ## mk the 2x2 table
            (total_marker <- length(which(df$cell_marker == cell_t)))
            (threshold_marker <- length(which(df$cell_marker[1:threshold_t] == cell_t)))
            (after_marker <- total_marker - threshold_marker)
            (threshold_non_marker <- threshold_t - threshold_marker)
            (after_non_marker <- nrow(df) - threshold_t - after_marker)
            
            fisher_tbl <- data.frame(x1 = c(threshold_marker, after_marker),
                                     x2 = c(threshold_non_marker, after_non_marker))
            row.names(fisher_tbl) <- c('before threshold', 'after threshold')
            colnames(fisher_tbl) <- c(cell_t, paste0('Non-', cell_t))
            
            for (alt in c('greater', 'less')){
                x <- fisher.test(fisher_tbl, alternative = alt)

                df_result <- data.frame(input_file = f_mj, 
                                        NA_genes_removed = msg_na,
                                        total_genes = nrow(df),
                                        threshold_method = threshold_method, threshold= threshold,
                                 genes_under_threshold = threshold_t, 
                                 cell_type = cell_t,
                                 marker_under_threshold = threshold_marker,
                                 non_marker_under_threshold = threshold_non_marker,
                                 marker_over_threshold = after_marker,
                                 non_marker_over_threshold = after_non_marker,
                                 Fisher_exact_test = alt,
                                 p_value = x$p.value,
                                 odds_ratio =  x$estimate)
                df_all <- rbind(df_all, df_result)
                #print(paste0(cell_t, ' loop2 ', alt))
            }
        }
    }
    return(df_all)
}


#---------------------------------------------------------------------------#
# Main function
#---------------------------------------------------------------------------#

mainCellTypeTest <- function(f_marker_genes, f_mj_ls, threshold_method_ls, threshold_ls, out_dir, 
                             f_result_out = 'celltype_fisher_test_results.tsv', rm_na = T){
    df_marker <- markerGenes(f_marker_genes)
    df_all <-NULL  # a df to store all test restuls for each file and each method
    for(i in 1: length(f_mj_ls)){
        ## loop for each meta or jack file
        (f_mj <- f_mj_ls[i])
        (prefix <-unlist(strsplit(f_mj, split = "/|\\.")))
        (prefix <- prefix[length(prefix)-1])
        
        for(j in 1:length(threshold_method_ls)){
            #### loop 2 for each threshold method
            (threshold_method <- threshold_method_ls[j])
            threshold <- threshold_ls[j]
            print(paste0('testing threshold method, ', threshold_method, ', threshold: ', threshold))
            
            (f_test_out <-paste0(out_dir, prefix, '_celltype_fisher_test_threshold_', threshold_method,'.txt')) 
            (f_out <-paste0(out_dir, prefix, '_celltype.tsv')) # out f for meta/jack file with celltype marker annotated
            
            if(j==1){ # save the annotated file once

                df_result <- cellTypeFisherTest(f_mj, df_marker, threshold = threshold, f_test_out=f_test_out,
                                   threshold_method = threshold_method, 
                                   save_f =T, f_out = f_out, rm_na = rm_na)
                df_all <- rbind(df_all, df_result)
            }else{
                
                df_result <- cellTypeFisherTest(f_mj, df_marker, threshold = threshold, f_test_out=f_test_out,
                                   threshold_method = threshold_method, 
                                   save_f =F, rm_na = rm_na, f_out = f_out)
                df_all <- rbind(df_all, df_result)
            }
        }
    }
    
    # write the result df
    file_out <-paste0(out_dir, f_result_out)
    print("WRITE TABLE")
    writeTable(df_all[order(df_all$p_value), ], file_out, msg = paste0("# ", Sys.Date())) ## order by pvalue
}





