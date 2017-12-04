# 2016-09-01

## corrected for gender effect, and remove study intercepts



## intercept(b0) + b1*disease-stage + b2*gender + error
#' intercept is different by study
#' b1: the slope of disease-stage (ie the FE of disease stage), the same for all study in randome intercept model(WT is the baseline)
#' b2: the slope of gender (i.e. FE_gender), male is the baseline


#' note that the original array_dat has more rows that the mixed model result(df_all), because mixed model remove gene when <2/3 of 
#' the study contain the gene expression
#' the final corrected matrix will have nrow of df_all genes and ncol of samples

# make a matrix of rows of gene(NA gene removed) * columns of samples for the correction value that remove intercept and remove gender value
# then add the original array dat to the correction matrix to get the corrected matrix
repCol<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

correctInterceptGender <- function(mm_dir,f_original_genotype,correct_gender=F){
    
    (rdata <- paste0(mm_dir, 'mixed_model_results.Rdata'))
    if(!file.exists(rdata)){
        stop(paste0(rdata, ' does not exist'))
    }
    result_dir = mm_dir
    load(rdata)
    
    if(nrow(array_dat) != nrow(df_all_expression)){
        
        if(all(rownames(df_all_expression) %in% rownames(array_dat))){
            print('get only result expression, continue')
            array_dat <- array_dat[as.character(rownames(df_all_expression)), ]%>%droplevels()
        }else{
            stop(print('array dat and result genes dont match, STOP'))
        }
    }
    
    ##########
    ## get the intercept for all studies, and make a correction matrix
    ##########
    (intercept_ls <- grep('intercept', colnames(df_all), value =T))
    study_ls <- gsub('_intercept', '', intercept_ls)
    
    intercept_m <- NULL
    print('Correct Intercept')
    for(i in 1:length(study_ls)){
        (study = study_ls[i])
        print(study)
        (intercept = intercept_ls[i])
        ## get the samples of the study
        (sample_df <- filterContain(array_design, 'Study', study))
        ## get the intercept (the same for all samples in the study) and
        ## make the matrix
        df <- repCol(df_all[, intercept], n = nrow(sample_df))
        colnames(df) <- sample_df$Sample
        rownames(df) <- rownames(df_all)
        df <- as.data.frame(df)
        if(is.null(intercept_m)){
            intercept_m <- df
        }else{
            intercept_m <- cbind(intercept_m, df)
        }
    }
    
    
    ##########
    ## get gender correction matrix (male is the baseline =0), only female need adjustment
    ## correction (i.e. slope for female) is the same for all studies and only apply to female
    ##########
    ##########
    ## correct matrix
    ##########
    ## make sure all the samples are in the same order
    array_dat_corrected <- array_dat[rownames(df_all), colnames(intercept_m)]
    
    if(correct_gender){
        ## skip if gender slope is NA
        if(!any(is.na(df_all$Gender))){
            print('Correct Gender')
            ## get default
            sample_m <- as.data.frame(repCol(df_all[, 'FE_Gender'], n = nrow(array_design)))
            colnames(sample_m) <- array_design$Sample
            rownames(sample_m) <- rownames(df_all)
            
            ## replace all male to 0
            for(i in 1:nrow(array_design)){
                if(array_design$Gender[i] == 'M'){
                    sample_m[, i] <- 0
                }
            }
            ## remove intercept and add gender effect for female
            df <- as.data.frame(as.matrix(array_dat_corrected) - as.matrix(intercept_m)+ as.matrix(sample_m))
        }
    } else {
            print('No Gender effect to correct')
            sample_m <- as.data.frame(repCol(rep(0, nrow(df_all)), n = nrow(array_design)))
            ## remove intercept and add gender effect for female
            df <- as.data.frame(as.matrix(array_dat_corrected) - as.matrix(intercept_m))
            sample_m=NULL
    }
 

    

    
    
    ## reorder by the sample order in design
    array_design <- orderCol(array_design, cols = c('Study', 'Disease_stage'))
    (x <- sort(levels(array_design$Study)))
    
    array_design_corrected <- array_design %>%
        mutate(Study =  factor(Study, levels = x)) %>%
        arrange(Study) 
    
    array_design_corrected$Disease_stage <- as.factor(array_design_corrected$Disease_stage)
    array_design_corrected <- array_design_corrected %>%
        mutate(Disease_stage =  factor(Disease_stage, levels = c("WT", "Disease"))) %>%
        arrange(Study, Disease_stage)
    
    array_dat_corrected  <- df[, array_design_corrected$Sample]%>%droplevels()
    
    array_dat <- array_dat_corrected
    #array_design <- array_design_corrected[, c("Study","Sample","Disease_stage")]%>%droplevels()
    array_design <- array_design_corrected
    
    ####
    ## add the original genotype to the array design
    ####
    df <- read.delim(f_original_genotype, comment.char = '#')
    df <- rmDup(df)
    
    ## check if the sampels names matched, if not use the simplified sample names
    if(!(all(df$Sample %in% array_design$Sample))){
        print('match without month')
        df <- df[, c('Sample_simplified', 'Original_genotype')]
        colnames(df) <- c('Sample', 'Original_genotype')
    }
    
    array_design <- noWarnings(left_join(array_design, df[, c('Sample', 'Original_genotype')]))
    if(any(is.na(array_design$Original_genotype))){## if still na left,then use genotype)
        array_design$Original_genotype <- array_design$Original_genotype
    }
    
    array_design$Original_genotype <- as.factor(array_design$Original_genotype)
    
    msg <- paste0(msg, '\n ', Sys.Date(), '\n gender is corrected and intercept is subtracted')
    f=paste0(result_dir, "mixed_model_results_exp_corrected.Rdata")
    if('df_all_expression' %in% ls()){
        save(array_dat,array_design, df_all, df_all_expression, msg, intercept_m, sample_m, file=f)
    }else{
        save(array_dat,array_design, df_all, msg, intercept_m, sample_m, file=f)
    }

    #return(list(array_design, df ))
}


