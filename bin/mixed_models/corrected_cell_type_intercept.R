# 2016-10-027

## corrected for cell_types, and remove study intercepts

#'
#'


## intercept(b0) + b1*disease-stage + b2*microglial + b3* neuros +.. + error
#' intercept is different by study
#' b1: the slope of disease-stage (ie the FE of disease stage), the same for all study in randome intercept model(WT is the baseline)
#' b2,3,4...: the slope of allcelltypes


#' note that the original array_dat has more rows that the mixed model result(df_all), because mixed model remove gene when <2/3 of 
#' the study contain the gene expression
#' the final corrected matrix will have nrow of df_all genes and ncol of samples

# make a matrix of rows of gene(NA gene removed) * columns of samples for the correction value that remove intercept and remove gender value
# then add the original array dat to the correction matrix to get the corrected matrix
repCol<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

correctInterceptCellType <- function(mm_dir,f_original_genotype,correct_gender=F, cell_markers_f, cell_types){
    
    (rdata <- paste0(mm_dir, 'mixed_model_results.Rdata'))
    if(!file.exists(rdata)){
        stop(paste0(rdata, ' does not exist'))
    }
    result_dir = mm_dir
    load(rdata)
    
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
    
    #-------------------------
    ## adj for cell types
    #-------------------------
    print(paste0('Correct for cell types: ', paste0(cell_types, collapse = ', ') ))
    
    cell_markers <-read.delim(cell_markers_f, comment.char = '#')
    
    ## double check sample names of cell type and expression matrix
    if(length(unique(as.character(cell_markers$Sample[which(cell_markers$Sample %in% array_design$Sample)]))) != length(unique(as.character(array_design$Sample)))){
        stop(print('samples names do not match'))
    }else{
        print('samples names match, continue')
    }
    
    ## get only the samples
    cell_markers <- cell_markers[which(cell_markers$Sample %in% array_design$Sample), ] %>%droplevels()
    
    ### get the cells
    
    
    cell_markers <- cell_markers[which(cell_markers$cell_type %in% cell_types), ]%>%droplevels()
    
    #get the PC1 scaled for each cell type (from long to wide table)
    df_markers <- cell_markers[, c('Sample', 'cell_type', 'PC1_scaled')]%>%droplevels()
    df_markers <- rmDup(df_markers)
    #df_markers <- tidyr::spread(df_markers, cell_type, PC1_scaled)
    
    df_markers <- noWarnings(reshape2::dcast(df_markers,Sample~cell_type))
    rownames(df_markers) <- df_markers$Sample
    if(!all(colnames(array_dat_corrected) %in%c(rownames(df_markers)))){
        return(df_markers)
        stop(print(paste0(cell_markers_f, ' does not contain all samples in the array')))
    }else{
        ## reorder cell proportions by sample names
        df_markers <- df_markers[colnames(array_dat_corrected), ]
    }
    
    adj_m <- NULL
    ## get the beta of each marker
    for(marker in intersect(cell_types, colnames(df_all))){
        print(paste0('Adjusting ', marker))
        proportion_m <- repCol(df_markers[, marker], n= nrow(array_dat_corrected))
        beta_m <- repCol(df_all[, marker], n= nrow(df_markers))
        
        if(is.null(adj_m)){
            adj_m <- t(proportion_m) * beta_m
        }else{
            adj_m <- adj_m + t(proportion_m) * beta_m
        }
        
    }
    
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
            df <- as.data.frame(as.matrix(array_dat_corrected) - as.matrix(intercept_m)- as.matrix(sample_m) - as.matrix(adj_m))
        }
    } else {
        print('No Gender effect to correct')
        sample_m <- as.data.frame(repCol(rep(0, nrow(df_all)), n = nrow(array_design)))
        ## remove intercept and add gender effect for female
        df <- as.data.frame(as.matrix(array_dat_corrected) - as.matrix(intercept_m)- as.matrix(adj_m))
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
        save(array_dat,array_design, df_all, df_all_expression, msg, intercept_m, sample_m,df_markers, adj_m, file=f)
    }else{
        save(array_dat,array_design, df_all, msg, intercept_m, sample_m, df_markers, adj_m,file=f)
    }
    print(f)
    
    #return(list(array_design, df ))
}


