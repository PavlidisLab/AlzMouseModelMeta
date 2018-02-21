# 2016-05-26
# update 2016-06-21 # add mixed model format

#' @param input_folder : folder for meta, jack files and mixed model files
#' @param bg_folder: folder for backgrounds
#' @param erminej_dir: dir for output sh and result dir for sh
#' @param scorecol, define scorecol (which column is the score column), if NULL, then auto select for meta, jack and mixed model for each disease

mkErminejSH <- function(disease,input_folder, bg_folder, erminej_dir,
                        scorecol = NULL,
                        xml,
                        maxsize = 100, iteration = '200000'){
    cmd <- paste0('#!/bin/bash',
                  '\n# ErmineJ script for ', disease, ', ', Sys.Date(),
                  '\n# Input folder: ', input_folder,
                  '\n# Background folder: ', bg_folder,
                  '\n# Result out dir: ', erminej_dir)
    
    
    cmd <- paste0(cmd, '
                  \n#############
                  \n# settings 
                  \n#############
                  \niteration=', iteration, 
                  '\nORAthreshold=0.1
                  \ngsrmethod=PRECISIONRECALL   # set method for GSR
                  \nmaxsize=', maxsize,
                  '\nminsize=10	# set the min class size
                  \nxml=', xml,
                  '\ntestmethod=GSR  #choose one of GSR (need to use -m for method), ORA(need to use -t for threshold, default =0.001), ROC, CORR')
    dir.create(erminej_dir, recursive = T, showWarnings = F)
    
    (mj_f_ls <- list.files(input_folder, pattern = 'mixed_model|jackknife|meta_genes.*tsv'))
    (bg_f_ls <- list.files(bg_folder, pattern = '.tsv', full.names = T))
    for (i in 1:length(mj_f_ls)){
        mj_f <- mj_f_ls[i]
        print(mj_f)
        (infile <- paste0(input_folder, '/',mj_f))
        (background <- grep(mj_f, bg_f_ls, value = T)) ## for meta and jack style, and jack mm
        if(length(background) ==0){
            (keyword <- gsub('mixed_model_jackknife_random_intercept_include_NA_low_exp_rm_jackknife_', '', mj_f)) 
            (keyword <- gsub('_mixed_model.tsv$|_mixed_model_results.tsv$', '', keyword)) ## for jack mixed model
            (keyword <- gsub('jackknife_', '', keyword)) ## for jackknife mixed model results
            (background <- grep(keyword, bg_f_ls, value = T))
        }
        
        if(length(background) ==0){ # this is for combined all disesase output
            (keyword <- gsub('random_intercept_include_NA_|random_slope_include_NA_', '', mj_f)) ## use jackknife for mixed model file (slightly less genes than meta background)
            (keyword <- gsub('low_exp_rm_', '', keyword)) ## use jackknife for mixed model file (slightly less genes than meta background)
            (keyword <- gsub('_mixed_model.tsv', '.*jackknife', keyword)) ## use jackknife for mixed model file (slightly less genes than meta background)
            (background <- grep(keyword, bg_f_ls, value = T))
        }
        
        if(length(background) ==0){
            for(keyword in c('early_up', 'early_down', 'late_up', 'late_down')){
                if(grepl(keyword, mj_f)){
                    (background <- grep(keyword, bg_f_ls, value = T))
                }
            }
        }
        
        print(background)
        if(length(background) ==0){
            stop(paste0('no background is found in ', bg_f_ls))}
        
        (outfile <- paste0(erminej_dir, 'ermineJ_', mj_f))
        if(is.null(scorecol)){
            if(grepl('jackknife', mj_f)){
                scorecol=2   # jackknife score (adj_combined_max_p)
            }else if(grepl('meta', mj_f)){# meta score (Fisher)
                scorecol = 7
            }else{scorecol = 8} # mixed model (up or down p value): colname: up_pval, or down_pval
        }

        
        cmd <- paste0(cmd, '\n\n\n############\n# ', mj_f,
                      '\n############\n\n',
                      '\nbackground=', shQuote(background),
                      '\ninfile=', shQuote(infile),
                      '\noutfile=', shQuote(outfile),
                      '\nscorecol=', scorecol, '\n\n',
                      "\nsh \"$ERMINEJ_HOME/bin/ermineJ.sh\" -a \"$background\" -c $xml -n \"$testmethod\" -o \"$outfile\" -s \"$infile\" -e \"$scorecol\" -g BEST -x $maxsize -y $minsize -j -i $iteration -m $gsrmethod -t $ORAthreshold -l")
        
    }
    file_name <- paste0(erminej_dir, '/', disease, '_erminej_script.sh')
    sink(file_name, type="output")
    writeLines(cmd)
    sink()
    print(paste0('Saved sh: ', file_name))
    return(file_name)
}


