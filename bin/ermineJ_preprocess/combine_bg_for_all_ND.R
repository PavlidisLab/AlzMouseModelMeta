#' 2016-08-11 make ermineJ background for all diseases
#' updated 2016-09-23: use the background from mixed model jackknife results
#' combine all the lowly expressed removed background of all 3 disesaes in 1, or any given diseases

#' process_keyword '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories



dir.create(out_dir,recursive = T, showWarnings = F)
for(meta_jack in meta_jack_ls){
    for (phase in phase_ls){
        for(regulation in regulation_ls){
            f=vector()
            for (i in 1:length(disease_ls)){
                ## get all the background files
                (f[i] <- paste0(home_dir, '/', disease_ls[i], "_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/ermineJ_background", process_keyword,"/",
                               'bg_jackknife_', phase,'_', regulation, '_regulation_',meta_jack,'_results.tsv'))
            }
            for(i in 1:length(f)){
                df <- read.delim(f[i], comment.char = '#')
                if (i ==1){
                    df_all <- df
                }else{
                    df_all <- noWarnings(full_join(df_all, df))
                }
            }
            
            ## remove duplicated entry:
            index <- which(!duplicated(df_all$geneSymbol))
            df_all <- df_all[index,] %>%droplevels()
            
            
            ## save the background file
            (f_out = paste0(out_dir, 'bg_low_expr_rm_', phase,'_', regulation, '_regulation_',meta_jack, '.tsv'))
            msg <- paste0('# ', Sys.Date(), ': ermineJ background for ', phase, ', ', regulation, ' regulated',
                          '\n# Input diseases : ', paste0(disease_ls, collapse = ', '), 
                          '\n# ', paste0(f, collapse = '\n# '))
            
            writeTable(df_all, f_out = f_out, msg = msg)
            print(f_out)
        }
    }
    
}





