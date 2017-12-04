#' 2016-09-25
#' take the top up and down genes from mm jack results
#' and make a new expression.Rdata that only contains these genes
#' in mixed_model/random_intercept_include_NA_low_exp_rm_top_genes/early(late)
#' redo mixed model later


library(HelperFunctions)
library(dplyr)




for(disease in disease_ls){
  source('config/config_wrappers.R')
  (jack <- paste0(disease_dir, 'mixed_model_jackknife/random_intercept_include_NA_low_exp_rm/'))
  jack_ls <- list.files(jack, pattern = 'mixed_model_results.tsv', full.names = T)

  

  for(phase in phase_ls){
    ## get genes from up and down regulation mm jack results
    genes <- NULL
    for(regulation in regulation_ls){ 
      (f <- grep(paste0(phase, "_", regulation), jack_ls, value = T))
      print(f)
      df <- read.delim(f, comment.char = '#')
      df <- orderCol(df, cols = 'final_jackknife_rank')  # order by jackknife rank
      genes <- unique(c(genes, as.character(df$geneSymbol[1:threshold])))
    }
    
    ## make new expression.Rdata for the top genes (redo the mixed model to get all estimates)
    (mm_dir <- paste0(disease_dir, 'mixed_model/random_intercept_include_NA_low_exp_rm/', phase, '/'))
    (out_dir <- paste0(disease_dir, 'mixed_model/random_intercept_include_NA_low_exp_rm_top_genes/', phase, '/'))
    dir.create(out_dir, recursive = T, showWarnings = F)
    (expression <- grep('expression.Rdata', list.files(mm_dir, full.names = T), value = T)) 
    load(expression)
    
    ## get the expression for the genes only
    array_dat <- array_dat[genes, ]%>%droplevels()
    array_dat_not_qn <- array_dat_not_qn[genes, ]%>%droplevels()
    plot_dir <- out_dir
    
    msg <- paste0(Sys.Date(), ': Top ', threshold, 'genes from up and down regulation. \n', msg)
    
    save(array_dat_not_qn, array_dat, array_design, phase, msg, plot_dir,
         file = paste0(plot_dir, 'expression.Rdata'))
    print(paste0(plot_dir, 'expression.Rdata'))
    }
}
