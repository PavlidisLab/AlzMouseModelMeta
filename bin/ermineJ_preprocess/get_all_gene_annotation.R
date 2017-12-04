# 2016-02-09
# last update: 2016-03-16
# use wrapper_filter_out_low_enrichment.R

## make a table of geneSymbol and full gene names of specified platforms
## make a table with genesymbols and GO terms for ermineJ


#' @note the final output
#'  the input is a list of platforms, and generate a summary file
#'  with all gene symbols (ignore probe names from different platforms) and their correponding gene name and GO terms
#'  if the same gene symbol has different GO terms, or some with empty GO entry, only the first one with the GO term
#'  is retained.
#' @note used in scripts filter_out_low_exp_for_enrichment.R and other scripts 
#'  that want to match geneSymbol with names and GO terms

 
source("helper_functions.R")



#########################################################
## helper function
#########################################################
aggregateTerms <- function(geneSymbol){
    #' return a df with geneSYmbol and GOterms. GOterms are aggregated among multiple entries and duplicated GO terms are removed
    #' @param geneSymbol: a dataframe with geneSymbol and GO terms for aggregation
    
    aggregateGO <- function(a_string, delim = "\\|"){
        (x <- sort(unique(unlist(strsplit(as.character(a_string), split = delim)))))
        return(paste0(x, collapse="|"))
    }
    
    df <- geneSymbol
    if(!all(c("GOTerms", "geneSymbol") %in% colnames(df))){
        stop(paste0("input doesn't has both columns GOTerms and geneSymbol. Columns:", paste0(colnames(df), collapse =", " )))
    }
    new_df <- aggregate(GOTerms ~ geneSymbol, df, function(x) aggregateGO(x))
    return(new_df)
}
    

#########################################################
## main function
#########################################################

#' @param platform_ls: list of files of gene annotations (in gemma format) to combine
#' @param platform_folder: folder for the platforms
#' @param ermineJ_format: in ermineJ format (with GO terms)
#' @param f_out: output dir
#' @return write a table with genesymbols, genenames, GO id etc
geneAnno <- function(platform_ls, platform_folder, f_out,  df_return = F,ermineJ_format = F){
    df_all <-T
    msg_platform <- ""
    for (platform in platform_ls){
        msg_platform <- paste0(msg_platform, "\n# ", platform)
        x <- processPlatformAnnotation(platform, platform_folder)
        annotation <- x[[1]]
        annotation$platform <- platform
        # change the probe names to gene names
        annotation$ProbeName <- annotation$GeneSymbols
        colnames(annotation)[1:3]<-c("geneID","geneSymbol", "geneNames")
        
        if(class(df_all) == "logical"){
            df_all <- annotation
        }else{
            df_all <- rbind(df_all, annotation)
        }
    } # end of loop with a prefiltered df of all probes
    
    #************
    # clean the entries (rm no genesymbol, duplicated entries)
    #************
    # rm rows without a genesymbol
    df_all <- excludeMatch(df_all, "geneSymbol", "")
    # rm duplicates
    df_all <- rmDup(df_all)
    
    #************
    # get GO terms
    #************
    ## make sure only 1 GO entry associate with 1 gene
    df_go <- aggregateTerms(df_all)
    
    tempFunction <- function(x){
        # select a column (x) in df_all, and return a df with aggregated x column by geneNames
        df <- rmDup(df_all[, c("geneSymbol", x)])
        (cmd <- paste0("df <- aggregate(", x, " ~ geneSymbol, df, function(x) paste0(x, collapse='|'))"))
        eval(parse(text = cmd))
        return(df)
    }
    
    
    #************
    # get Gene names, and other ids
    #************
    df_genenames <- tempFunction("geneNames")
    df_gemmaid <- tempFunction("GemmaIDs")
    df_ncbi <- tempFunction("NCBIids")
    df_platform <- tempFunction("platform")

    if(!all(c(nrow(df_gemmaid), nrow(df_genenames), nrow(df_ncbi)) == nrow(df_go))){
        stop(paste0("not all aggregated GemmaID, Genenames, and NCBI id match with number of entries for GO terms"))
    }
    
    #************
    # combine all info to a table and write table
    #************
    df <- noWarnings(Reduce(left_join, list(df_genenames, df_go, df_gemmaid, df_ncbi, df_platform)))
    if(ermineJ_format){
        df <- cbind(data.frame(geneID = df$geneSymbol), df)
    }
    
    # write to table
    msg <- paste0("#", Sys.Date(), "\n# Included platforms: ", msg_platform)
    sink(f_out, type="output")
    writeLines(msg)
    sink()
    write.table(df, file =f_out, quote = F, sep ='\t',
                row.names =F, append =T)
    print(paste0("output file ", f_out))
    
    if(df_return){
     return(df)   
    }
}


