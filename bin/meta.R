#' 2016-01-01
#' update 2016-06-13
#' ################################################################################################
#'
#' 
#' input: up or down 1-sided p values of genes (not adjusted!)
#' output: meta analysis combined p values by specified method
#' 
#' Notes: from limma.R
#'            - probes not mapped to a gene are removed
#'            - convert gene symbol as id. If multiple probes mapped to a gene, 
#'             get the min p value of one of the probes. 
#'            - If a probe mapped to more than 1 gene, which means the probe is not specific 
#'                and will be removed from the analysis (in the processInput())
#'             - probes are aggregated into genes
#' 
#' ----------------------------------------------------- #
#' USAGE: 
# source('meta.R')
# f_meta <- '/home/bzhuang/gemma_data/explore_data/meta_analysis_pval/'
# keyword <- "down.*3XTG"
# output <- '/home/bzhuang/gemma_data/explore_data/meta_analysis_result/meta_analysis_down_3XTG.tsv'
# metaPval(f_meta, keyword, output, method = 'Fisher', output_df =F)
#' ----------------------------------------------------- # 
#' 
#' Done: change to 1-sided p values, not padj -done
#' Done: probes mapped to 2 genes are removed, but insulin may got mapped Ins1|Ins2, may want to keep these probes -done
#' Done: (in limma_DEA.R)In the case where a gene had more than one probeset assigned to it, the p-values for 
#'             the probesets are Bonferroni-corrected and the lowest corrected p-value was used to represent the gene for that dataset
#' Done: add the jacknife method
#' 
#' UPDates:
#' 2016-06-13: add jack_rm_n in the jackknife procesdure, choose to remove n files for each run (default remove 1 at a time), add better msg


#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
#install.packages('MetaDE')
library('MetaDE')


#---------------------------------------------------------------------------#
# PART 1: SET SYS DEFAULTS
#---------------------------------------------------------------------------#
# load packages
package_ls <- c('plyr', 'dplyr','ggplot2', 'MetaDE')
for(i in package_ls){
    cmd <- paste0("if('package:", i, "' %in% search() == FALSE){library(", i, ")}")
    eval(parse(text =cmd))
}

 
library("HelperFunctions")

## set defaults
## image size (pch)
width = 1200
height = 800

rm(package_ls)
rm(i)
rm(cmd)

if(!exists("adj_method")){adj_method <- "BH"}

#---------------------------------------------------------------------------#
# PART 2: SET INPUTS
#---------------------------------------------------------------------------#
#f_meta <- '/home/bzhuang/gemma_data/explore_data/meta_analysis_pval/'
#f_to_process <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/2016-02-10/AD_late_6_15months/GSE1556_gene_down_regulation_1_sided_pval_GenotypeTg2576-vs-baseline_WT.tsv'
#---------------------------------------------------------------------------#
# PART 3: FUNCTIONS
#---------------------------------------------------------------------------#
processInput <- function(f_to_process, plot = F, verbose =T){
        #' f_to_process (str) a file with 2 columns: pvalue and genesymbols
        #' return a df with 1 column and genesymbols as row names
        df <- read.delim(f_to_process ,check.names = F, comment.char = "#")
        
        rownames(df) <- df$GeneSymbols
        col_name <- grep("GeneSymbols", colnames(df), invert=T, value =T)
        df2 <- data.frame(df[, col_name])
        colnames(df2) <- col_name
        rownames(df2) <- df$GeneSymbols
        return(df2)
}


#' Grep files from a given folder by keyword
#' @param f_meta folder dir for 1-sided pval from each study (gene pvalues) 
#' @param keyword (string regex) to select the pval files, 
#'  e.g. "up_regulation", "down_regulation", "up.*Tg2576" etc.
getFiles <- function(f_meta, keyword){
    (file_list <- grep(keyword, ignore.case = T,
                       list.files(path = f_meta, pattern = '.tsv', full.names=T),
                       value = T))
    return(file_list)
}


#' Combine P values for all input files
#' 
#' Input a list of files contain 1-sided pval from each study (gene pvalues) (from Limma_DEA.R)
#' 
#' @param file_list (str list) list of files containing 1-sided pval(gene-level) from eacj study
#' @param gene_annotation (file dir) a file contain two columns, geneSymbol and geneNames
#' @param output (str): output dir for new table, default = "./"
#' @param method (str): method to combine pvalues (from MetaDE.pvalue())
#' @param output_df (bool): whether to output the dataframe as a variable, default =F
#' @param NA_rm: whether to remove probes with NA values in some datasets before combining p values, default = F
#' @param jackknife (bool) if use jackknife mode, a table of adj pvalues is returned as an object. default = F
#' @return A table with combined pvalues, adjusted combined p values,  and pvalues from each study (write to table)
#' @note
#' require MetaDE::MetaDE.pvalue() to calculate combined p values
#' The combined p-values are adjusted for multiple testing using the selected method (default is Benjamini-Hochberg) method. 
#' The genes that meet the threshold of FDR<0.05 were considered to be meta-signature genes.
metaPval <- function(file_list, gene_annotation, output = "./", method = 'Fisher', 
                     adj_method = "BH", output_df =T,NA_rm =F, write_df = T, jackknife =F){
    df_all <- T
    for(i in (1: length(file_list))){
        f_to_process <- file_list[i]
        df <- processInput(f_to_process) 
        if(!is.data.frame(df_all)){
            df_all <- df
        } else {
            df_all <- mapBindAllColumns(df_all, df)
            if(NA_rm){
                df_all <- na.omit(df_all)
            }
        }
    }
    
    ## get the combined p value by MetaDE::MetaDE.pvalue()
    x <- list(p=df_all)    #
    pval <- MetaDE.pvalue(x, meta.method = method)
    combined_p_df <- as.data.frame(pval$meta.analysis[2]) ## get the meta pvalue

    ## get the adjusted p value
    combined_p <- as.numeric(unlist(pval$meta.analysis[2]))
    adj_p_df <- as.data.frame(p.adjust(combined_p,adj_method))
    colnames(adj_p_df) <- paste0(method, "_", adj_method, "_adjusted")
    
    ## jackknife mode, only return adj combined p values
    if(jackknife){
        row.names(adj_p_df) <- row.names(df_all)
        return(adj_p_df)
    }
    
    
    ## get the gene annotation
    df_gene <- read.delim(gene_annotation, comment.char="#")
    df_tmp <- data.frame(geneSymbol = row.names(df_all))
    df_tmp <- noWarnings(left_join(df_tmp, df_gene))
    
    ## combine all columns
    pval_df <- cbind(df_tmp, combined_p_df, adj_p_df, df_all)    ## add gene name
    (cmd <- paste0("pval_df <- pval_df[order(pval_df[ , '",method,"']), ] "))
    eval(parse(text = cmd))   ## order by meta p
    
    #plot(density(pval_df[,2]))
    #hist((pval_df[,2]))
    
    ## add the rank
    (rank_col <- intersect(c(method, 'adj_combined_max_p'), colnames(pval_df)))
    pval_df$Rank <- rank(pval_df[, rank_col])
    
    
    ## write to file
    if(write_df){
        tmp_msg <- paste(file_list, collapse="\n# ")
        tmp_msg <- paste0("#", Sys.Date(), "\n# combine Pvalue method: ", method,
                          "\n# Input files for meta-analysis\n#", 
                          tmp_msg, "\n#")
        sink(output, type="output")
        writeLines(tmp_msg)
        sink()
        
        noWarnings(write.table(pval_df, file = output, row.names = F,
                               sep ='\t', quote = F, append = T))
        print(paste0(output," saved"))
    }

    if(output_df){
        return(pval_df)    
    }
}

# 
# 164 To obtain core signature genes we employed a jackknife procedure, which performs n sub-meta-analyses, 
# 
# 165 where n is the number of data sets considered, removing sequentially one data set at a time and then 
# 
# 166 finally combining the results of all n runs. We combined the results by intersecting n resulting meta-
#     
#     167 signatures at FDR<0.1.

#' @param jack_rm_n remove n number of files from all files
#' @param file_list (str list) list of files containing 1-sided pval(gene-level) from eacj study
#' @param gene_annotation (file dir) a file contain two columns, geneSymbol and geneNames
#' @param output (str): output dir for new table, default = "./"
#' @param method (str): method to combine pvalues (from MetaDE.pvalue())
#' @param output_df (bool): whether to output the dataframe as a variable, default =F
#' @param NA_rm: whether to remove probes with NA values in some datasets before combining p values, default = F
#' @param jackknife (bool) if use jackknife mode, a table of adj pvalues is returned as an object. default = F
#' @return A table with combined pvalues, adjusted combined p values,  and pvalues from each study (write to table)
jackknifeProcedure <- function(file_list, gene_annotation, jack_rm_n = 1,
                               threshold = 0.1, output="./", method = 'Fisher', adj_method = "BH", 
                               output_df =F, NA_rm =F, write_df = T){
    
    print("Jackkinfe procedure")
    ## recalculate adj combined p for each subset of files
    jack_all <- T
    (m <- length(file_list) - jack_rm_n)  # the number files in each combination (rm some files each run)
    (all_combo <- combn(file_list, m))
    msg_jack <- ''
    
    for (i in 1:dim(all_combo)[2]){
        # get the files in the combn
        (jk_file_list <- all_combo[, i])

        jack <- metaPval(jk_file_list, method = method, adj_method = adj_method, NA_rm =NA_rm, jackknife =T)
        
        colnames(jack) <- paste0("jackknife_", i)
        if(!is.data.frame(jack_all)){
            jack_all <- jack
        } else {
            jack_all <- mapBindAllColumns(jack_all, jack)
        }
        print(jk_file_list)
        
        ## msg for each jack
        msg_jack <- paste0(msg_jack, '\n# ',"jackknife_", i, ': removed files: \n# ', 
                           paste0(setdiff(file_list, jk_file_list), collapse = '\n# '))
    }
    
    
    
#     for (i in 1:length(file_list)){
#         ## remove the ith file from the list
#         (jk_file_list <- file_list[-i])
#         jack <- metaPval(jk_file_list, method = method, adj_method = adj_method, NA_rm =NA_rm, jackknife =T)
#         colnames(jack) <- paste0("jackknife_", i)
#         if(!is.data.frame(jack_all)){
#             jack_all <- jack
#         } else {
#             jack_all <- mapBindAllColumns(jack_all, jack)
#         }
#     }
    

    tmpF <- function(x){
        ## get the max value, return NA if all elements are NA
        ## if number of NAs more than have of the jackknife procedure, return NA
        if(length(which(is.na(x) == T)) > length(x)/2){
            return(NA)
        }else{
            return(max(x, na.rm = T))
            }
    }
    
    ## get the genes with all p < threshold
    df <- data.frame(adj_combined_p = apply(jack_all, 1, function(x) tmpF(x)))  # get the max adj p for each row
    
    index <- which(df$adj_combined_p <= threshold)
    jack_df <- cbind(row.names(df)[index], as.data.frame(df[index, ]))
    
    
    jack_df <- cbind(row.names(df), as.data.frame(df), jack_all)
    colnames(jack_df)[1:2] <- c("geneSymbol", "adj_combined_max_p")
    
    ## reorder by p
    jack_df <- jack_df[order(jack_df$adj_combined_max_p), ]
    
    ## get the gene annotation
    df_gene <- read.delim(gene_annotation, comment.char = "#")
    jack_df <-noWarnings(left_join(jack_df, df_gene)) 
    
    if(NA_rm){
        jack_df <- na.omit(jack_df)
    }
    ## write to file
    if(write_df){
        (rank_col <- intersect(c(method, 'adj_combined_max_p'), colnames(jack_df)))
        jack_df$Rank <- rank(jack_df[, rank_col])
        
        
        tmp_msg <- paste(file_list, collapse="\n# ")
        tmp_msg <- paste0("# Jackknife Procedure: ", Sys.Date(), "\n# combine Pvalue method: ", method,
                          "\n# adjusted combined p value threshold: ", threshold,
                          "\n# Input files for meta-analysis\n#", 
                          tmp_msg, "\n#")
        if(NA_rm){
            tmp_msg <- paste0(tmp_msg, "\n# NA removed.")
        }
        tmp_msg <- paste0(tmp_msg,msg_jack) ## which file(s) removed each run
        
        sink(output, type="output")
        writeLines(tmp_msg)
        sink()
        
        noWarnings(write.table(jack_df, file = output, row.names = F,
                               sep ='\t', quote = F, append = T))
        print(paste0(output," saved"))
    }
    
    ## return all the jackknife values (before filtering)
    if(output_df){
        return(jack_all)    
    }
    
}
        

    
    
