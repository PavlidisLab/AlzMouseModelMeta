#' 2017-01-16
#' summary of all samples and genes counts for input data
#' genes are after filter
rm(list=setdiff(ls(),'home_dir'))
 

library(HelperFunctions)
library(dplyr)

summaryDf <- function(array_design, array_dat, phase){
    df <- data.frame(disease = disease, phase=phase, n_studies = length(levels(array_design$Study)),
                     control= table(array_design$Disease_stage)['WT'],
                     disease_samples = 
                         table(array_design$Disease_stage)['Disease'],
                     total = nrow(array_design),
                     total_genes = nrow(array_dat))
    return(df)
}





summaryDisease <- function(disease){
    phase ='early'
    
    robj <- paste0(home_dir, disease, "_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/", phase,"/mixed_model_results_exp_corrected.Rdata")
    
    load(robj)
    df <- summaryDf(array_design, array_dat,phase)
    array_design_e <- array_design
    array_dat_e <- array_dat
    
    phase ='late'
    robj <- paste0(home_dir, disease, "_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/", phase,"/mixed_model_results_exp_corrected.Rdata")
    
    load(robj)
    df_l <- summaryDf(array_design, array_dat,phase)
    array_design_l <- array_design
    array_dat_l <- array_dat
    
    length(union(array_design_l$Sample, array_design_e$Sample))
    
    df <- rbind(df, df_l)
    
    (adj = length(intersect(array_design_l$Sample, array_design_e$Sample))/2)
    
    df$control <- df$control-adj
    df$total <- df$total-adj
    
    
    tmp <- rmDup(rbind(array_design_l[, c('Study', 'Sample', 'Disease_stage')],array_design_e[, c('Study', 'Sample', 'Disease_stage')]))
    print(table(tmp[, c('Study', 'Disease_stage')]))
    print(table(tmp[, c('Study')]))
    
    
    
    ## return
    (study_count <- length(union(array_design_l$Study, array_design_e$Study)))
    gene_count <- union(rownames(array_dat_e), rownames(array_dat_l))
    return(list(df, study_count, gene_count))
}


disease='AD'
x <- summaryDisease(disease)
df <- x[[1]]
study <- x[[2]]
gene_count <- x[[3]]

disease='HD'
x <- summaryDisease(disease)
df2 <- x[[1]]
study2 <- x[[2]]
gene_count2 <- x[[3]]

df <- rbind(df, df2)


df2 <- data.frame(disease ='Total', phase ='', 
           n_studies = study+study2, control = sum(df$control),
           disease_samples = sum(df$disease_samples),
           total = sum(df$total),
           total_genes= length(union(gene_count, gene_count2)))
df <- rbind(df, df2)

# save file
f_out <- paste0('../../results/ND_results/input_data_summary_', Sys.Date(), '.tsv')
writeTable(df=df, f_out= f_out)
print(f_out)
