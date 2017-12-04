#' 2016-03-15
#' ## input the jackknife/or meta files and get the GOterms for these genes (for ermineJ background of prefiltered probes)
#' remove genes that are not in the annotation files (e.g. if annotation is BP only, some probes are lost (but why?))
#' use wrapper
# 'low_exp_rm_analysis/wrapper_filter_low.R' # for prefiltered genes with gender gene threshold
# 'ermineJ_preprocess/wrapper_filter_out_low_enrichment.R' # for prefiltered genes with hard threshold

getGOterms <- function(f_anno, f_up, f_down, f_out, notes = ""){
    #' f_anno, a file with all the genes and GOterms
    df_anno<- read.delim(f_anno, comment.char="#")
    df_up <- read.delim(f_up, comment.char = "#")
    df_down<- read.delim(f_down, comment.char = "#")
    
    # get the union of genes from meta/jackknife files
    df_all <- rbind(df_up[, c('geneSymbol', 'geneNames')], df_down[, c('geneSymbol', 'geneNames')])
    # rm duplicates
    df_all <- df_all[!duplicated(df_all), ]%>% droplevels()
    
    genes <- union(df_down$geneSymbol, df_up$geneSymbol)
    genes_anno <- intersect(df_anno$geneSymbol, genes)
    
    df <- df_all
    # get the go terms
    df <- noWarnings(left_join(df, df_anno[, c('geneSymbol', 'GOTerms')]))
    
    ## save
    msg <- paste0("#", Sys.Date(), "\n# background genes for ermineJ with GO terms",
                  "\n# Gene counts from files: ", length(genes),
                  "\n# Genes without annotation entry: ", length(genes) - length(genes_anno),
                  "\n# non-expressed probes filtered and NA filtered", 
                  "\n# background annotation: ", f_anno,
                  "\n# input files: ", f_up, "\n# ", f_down,
                  "\n# ", notes)
    sink(f_out, type="output")
    writeLines(msg)
    sink()
    noWarnings(write.table(df, file =f_out, quote = F, sep ='\t',
                row.names =F, append =T))
    print(paste0("output file ", f_out))
}
