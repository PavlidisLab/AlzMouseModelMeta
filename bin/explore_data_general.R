#' 2015-11-16
#' last update 2016-05-16
#' 
#' ################################################################################################
#' a general function to explore datasets
#' explore a new dataset, quality check, visualize data, to spot outlier, batch effect etc
#' 
#' INPUT: 
#'    - output_folder: default is'/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/'
#'    - f_dataset_ls: A file contain info about the dataset and whether to subset or remove outliers
#'                  default is "/home/bzhuang/AD_mouse_model_project/config_files/mouse_datasets_explore.tsv"
#'      - required columns:
#'      'eeName': GSE_id, e.g. GSE1234
#'      'platform_id': e.g. GPL1261
#'      'subset': 'yes' indicate that samples needed to be subsetted, 'no' otherwise
#'      'subset_by': if subset is 'yes', fill in this column of which factor needs to be subsetted 
#'                    and separated by '; '(with a space) if multiple factors are required to be subsetted. e.g. 'OrganismPart; Timepoint'
#'      'keep_subset': indicate the levels of a factor or multiple factors to keep for analysis. order is the same as order in 'subset_by'
#'                    column. e.g. "c('Hippocampus'); c('6_months', '15_months')"
#'    - f_design: an output file from the combine_experimental_design.R. 
#'                default is'/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv'
#'    - rm_outlier: if T, remove identified outlier by sample name. default is F
#'    - to_subset: if T, subset samples. default is T
#'    - explore_method: choose 'explore_mode'(default) or 'analysis_mode'
#' default not to subset the samples(by give factors and levels to keep in df_process:  to_subset <- T
#' 
#' ----------------------------------------------------- #
#' USAGE: to run by range of rows, specified rows, or datasets
# rm(list=setdiff(ls(),'home_dir'))
# rm_outlier <- T
# to_subset <- T
# f_dataset_ls <- "/home/bzhuang/AD_mouse_model_project/config_files/mouse_datasets_explore.tsv"
# output_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/'
# #unfiltered_exp_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/unfiltered_exp_data/'
# unfiltered_exp_folder <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/'
# explore_method <- 'explore_mode'
# start_row <- 1
# #to_do_ls <- 4
# #dataset_todo <-'GSE65067'
# source('explore_data_general.R')
#' ----------------------------------------------------- # 
#' Notes: 
#' 1. for expression matrix, column names(i.e. sample names) cannot contain'#', must change to '.' and 
#' the name should match with the processed desigen file too
#' 2. experimental design is read from '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv'
#' 3. none-log transformed data will automatically log2 transformed
#' 4. must contain at least 4 samples
#' 
#' depedencies
#' 'helper_functions.R'
#' 
#' processed data are saved as R object (for limma DEA)
#' ################################################################################################
#' 


#'*************************************************************************#
#' PART 1: SET SYS DEFAULTS
#' - Load packages, source and define defaults for img size, and set wd
#' PART 2: SET INPUTS
#' - Dirs for platforms, data matrix
#' - file: - a file contain the list of dataset GSE ids for processing (experiment list)
#'        - a gene list for sanity check (e.g. gender genes)
#' PART 3: DEFAULTS FOR THE MAIN FUNCTION
#' - default to whether remove outliers (default false)
#' - set index to run which dataset (default: all dataset in the input to process file)
#' PART 4: MAIN FUNCTION
#' - loop
#**************************************************************************#


#---------------------------------------------------------------------------#
# PART 1: SET SYS DEFAULTS
#---------------------------------------------------------------------------#
 
source('helper_functions.R')

## set defaults
## image size (pch)
width = 1600
height = 1000

#---------------------------------------------------------------------------#
# PART 2a: SET INPUTS
#---------------------------------------------------------------------------#

## input folders
# platform_folder <- paste0(disease_dir, '/data_and_QC/design_and_annotation/platforms/')
# design_folder <- paste0(disease_dir, '/data_and_QC/design_and_annotation/experimental_design_details/')

# ## output folder (each dataset has its own folder under the output folder)
# if(!exists('output_folder')){
#   output_folder <- paste0('/home/bzhuang/', disease, '_mouse_model_project/data_and_QC/gemma_data/explore_data/')
# }
# ## load the table with downloaded dataset and platform info and for loop to process
# if(!exists('f_dataset_ls')){
#   f_dataset_ls <- paste0('/home/bzhuang/', disease, '_mouse_model_project/config_files/mouse_datasets_explore.tsv')
# }


## design files for all datasets
# if(!exists('f_design')){
#     f_design <- paste0(disease_dir, '/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv')
# }

## a gene list for sanity heatmap check
f_gene_list <- f_sanity_check_genes


## set the explore_method
if(!exists("explore_method")){
    explore_method <- 'explore_mode'
}


#---------------------------------------------------------------------------#
# PART 2b: read df and create folders
#---------------------------------------------------------------------------#
dir.create(output_folder, showWarnings=F, recursive=T)
df_process <- read.delim(f_dataset_ls, comment.char="#")  # must add comment.char="#"

print(paste0("disease: ", disease, ". Datasets to process: ", f_dataset_ls))



#---------------------------------------------------------------------------#
# PART 3: DEFAULTS FOR THE MAIN FUNCTION
#---------------------------------------------------------------------------#
## default not to remove outliers
if(!exists("rm_outlier")){rm_outlier <- T}

## default not to subset the samples
if(!exists("to_subset")){to_subset <- T}



#---------------------------------------------------------------------------#
# PART 4: DEFINE HELPER FUNCTIONS
# 
# loopSampleOutlier(): find sample outliers by sample to sample correlation
# processGemmaDataset(): process genmma data matrix
#---------------------------------------------------------------------------#


############################
## detect sample outliers by sample to sample correlation (after remove already identified outliers,and after subsetting the samples)
############################
#' detect samples that are 2sd sample-sample correlation away from mean
#' @return a msg and df of the outlier samples
#' @param prefix a note to put in the notes column
sampleOutlier <- function(array_dat, array_design, notes=""){
    msg <- ""
    sample_cor <- cor(array_dat)
    ## replace diag cells with NA value
    diag(sample_cor) <- NA
    (lower_r <- mean(sample_cor,na.rm =T) - sd(sample_cor,na.rm =T)*2)
    (avg_cor <- apply(sample_cor, 1, function(x) mean(x, na.rm = T)))
    (index <- which(avg_cor < lower_r))
    if(length(index) > 0){
        df <- data.frame(sample_avg_correlation = avg_cor[index], 
                         mean_cor = mean(sample_cor,na.rm =T), sd = sd(sample_cor,na.rm =T),
                         outlier_notes = notes)
        df <- cbind(array_design[names(avg_cor[index]), ], df) # combine with sample infomation
        msg <- paste0(msg, 'Outliers found: ', paste0(names(avg_cor[index]), collapse = ", "))
        array_dat <- array_dat[, -index]%>% droplevels()
    }else{
        df <- NULL
        msg <- paste0(msg, 'No outliers detected')
    }
    returnlist <- list(msg, df, index, array_dat)
    return(returnlist)
}

#' loop iteratively to get the sample outliers
#' repeat the process each time by removing the outlier, until not outliers are detected
loopSampleOutlier <- function(array_dat, array_design){
    array_dat_tmp <- array_dat
    index = 1
    outlier_round = 0
    df_outlier <- NULL
    msg_outlier <- "\n#-------------------------------#\n## Sample Outliers\n#-------------------------------#\nSample outlier defined by 2 sd away from mean sample to sample correlation. \n"
    while(length(index) > 0){
        outlier_round = outlier_round +1
        x <- sampleOutlier(array_dat_tmp, array_design, notes = paste0('outlier_round_', outlier_round))
        (msg_x <- x[[1]])
        (df_x <- x[[2]])
        (index <- x[[3]])
        array_dat_tmp <- x[[4]]
        df_outlier <- rbind(df_outlier, df_x)
        msg_outlier <- paste0(msg_outlier, msg_x)
    }
    returnlist <- list(msg_outlier, df_outlier)
    return(returnlist)
}

############################
## Process a dataset expression
############################

#' Process a dataset expression
#' @description
#' for a dataset (e.g 'GSE1556')
#' reformat data matrix, log transform(if needed) and return an R object with the processed exp matrix, and array design
#' plot density and sample corr heatmaps 
#' 
#' @param dataset (str) the GSE id, e.g 'GSE1556'
#' @param unfiltered_exp_folder(str): folder for expression matrix downloaded from Gemma, or noramlized CEL exp files
#' @param annotation (df) processed platform annotation
#' @param plot_pre (str) folder and prefix of the plots
#' @param result_pre (str) folder and prefix of the results
#' @param model_test (bool) whether use gt() function to test the model (e.g. id batch effect etc), default = T
#' @return
#' A list containing 1. processed array_dat (df), 2. modified array_design(df), 3. 
#' save plots
processGemmaDataset <- function(dataset, array_design, unfiltered_exp_folder, annotation, plot_pre, result_pre, model_test = F, to_plot=T, to_log =T){  
    ###################  check inputs  ################### 
    if(!dir.exists(unfiltered_exp_folder)){
        stop(paste0("Error: Dir for expression files ", unfiltered_exp_folder, " is not found."))
    }
    if(class(annotation) != "data.frame"){
        stop(paste0("Input Annotation is a ", class(annotation), ". A data.frame is required."))
    }
    if(class(dataset) != "character"){
        stop(paste0("Input dataset is a ", class(dataset), ". A string is required, e.g. 'GSE1234'."))
    }
    if(class(plot_pre) != "character"){
        stop(paste0("Input plot_pre is a ", class(plot_pre), ". A string is required, e.g. '/home/folder/plot/prefix'."))
    }
    if(class(result_pre) != "character"){
        stop(paste0("Input result_pre is a ", class(result_pre), ". A string is required, e.g. '/home/folder/results/prefix'."))
    }
    ###################  load data  ###################             
    ## grep the exp data file dir (gemma or CEL RMA normalized files)
    (x <- grep(paste0(dataset,"_"), 
               list.files(unfiltered_exp_folder), value=T))
    if(length(x) == 0){
        stop(paste0("Error: data file of ", dataset, " is not found."))
    }
    if(length(x) != 1){
        stop(paste0("Error: more than 1 data files found: ", paste0(x, collapse = ", ")))
    }
    (data_file <- paste0(unfiltered_exp_folder, "/", x))
    print(paste0("Data file: ", data_file))
    
    ## turn the first column(probenames) as the row names, 
    ## must use check.names = F to avoid R changing special char to dots
    data <- read.table(data_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1,fill=T,quote="",
                       check.names = F, comment.char = "#")
    
    ## Merge annotation and data dataframes by probe name. turn the first column(probenames) as the row names
    data2 <- merge(annotation, data, by = "row.names", all=FALSE, sort=TRUE)
    if(nrow(data2) == 0){
        stop(paste0("maybe use the wrong annotation file? no match of probe names"))
    }
    row.names(data2) <- data2[,1]
    data2<-data2[,-1]
    
    ## the format for gemma expression file and normalized exp from raw CEL are different. Gemma sample names are
    ## descriptive samples names after the column "NCBIid"; CEL file names are using GSM id, after the column "NCBIids" (with the "s")
    ## both cases, should contain a column called "NCBIid" (from the annotation file)
    
    ## gemma style or GEO style:
    if("NCBIid" %in% colnames(data2)){
        ## get only expression data matrix
        ## expression column range (all columns after "NCBIid" in gemma, or "NCBIids" in CEL),
        (col_range <- (which(colnames(data2) == "NCBIid") +1): ncol(data2))
    }else{
        if(!("NCBIids" %in% colnames(data2))){
            stop("Error: NCBIid column is not in expression data table")
        } ## CEL style
        (col_range <- (which(colnames(data2) == "NCBIids") +1): ncol(data2))
    }
    colnames(data2)[col_range]
    array_dat<-data2[, col_range]
    
    ## simplfy col names for data matrix (to match with array design) and reorder columns
    x <- colnames(array_dat)
    
    if(all(grepl("^GSM", x))){      # CEL file contains GSM id(must start with)
        new_col <- plyr::mapvalues(x, from = as.character(array_design$ExternalID), 
                                   to = as.character(array_design$Sample), warn_missing=F)
    }else {  # other are by externalID
        new_col <- plyr::mapvalues(x, from = as.character(array_design$ExternalID), 
                                   to = as.character(array_design$Sample), warn_missing=F)
    }
    
    ## if not updated then it's gemma style
    if(all(x == new_col)){ # gemma style
        print('sample names updated')
    new_col <- plyr::mapvalues(x, from = as.character(array_design$Bioassay), 
                               to = as.character(array_design$Sample), warn_missing=T)}
    
    colnames(array_dat) <- new_col
    
    ## get col names from data design (Sample)
    (sample_order <- as.character(array_design$Sample))
    ## check if colnames match
    if(length(setdiff(colnames(array_dat), sample_order)) != 0){
        stop(paste0("Error: sample names in expression data and design matrix do not match ", 
                    paste0(setdiff(colnames(array_dat), sample_order), collapse=", ")))
    }
    
    array_dat <- array_dat[, sample_order]
    
    ###############################
    ### remove NaN columns
    ###############################
    ## remove NaN, NA columns of expression matrix (if the full column is NA or NaN)
    ## when NaN column is removed, the design file should also be updated (otherwise design file has more exp than data matrix)
    df <- array_dat
    df <- df[, colSums(is.na(df)) != nrow(df)] %>% droplevels()
    removed_col <- setdiff(colnames(array_dat), colnames(df))
    msg_col_rm <- ""
    if(length(removed_col) !=0){
        f_rm <-paste0(result_pre, "removed_samples.tsv")
        array_dat <- df
        ## record removed experiments
        array_design_removed <- array_design[which(array_design$Sample %in% removed_col), ] %>% droplevels()
        write.table(array_design_removed, file =f_rm, quote = F, sep ='\t',
                    row.names =F)
        msg_col_rm <- paste0(msg_col_rm, length(removed_col), " sample(s) removed ",  
                             paste(removed_col, collapse =", "),  " (no value recorded) \n see ", f_rm)
        ## remove the column in array design too
        array_design <- array_design[-which(array_design$Sample %in% removed_col), ] %>% droplevels()
        print(msg_col_rm)
        print(paste0("new sample size: ", nrow(array_design)))
    }
    
    ###############################
    ######### log2 transformation if not done
    ###############################
    # log transformation must be done before removing missing values. (if max(value >0) than consider not transformed)
    # for illumina data, values are usually not log transformed and background 
    # correction will create negative values (which can't be log transformed)
    # these negative values are treated as missing value (but not in the case of two color agilent, where the exp is ratio)
    msg_plot <-""
    msg_log <-""
    if(max(na.omit(array_dat))>100){
        # density plot before transformation
        main_title <- paste0('Sample Density: ', dataset, msg_col_rm,
                             '\n', ncol(array_dat), ' samples, ', nrow(array_dat), ' probes')
        png(filename= paste0(plot_pre,'density_before_transformation.png'), width = width, height = height)
        plotDensities(array_dat, legend=F, main = main_title)
        dev.off()
        
        ## transform data
        (msg_plot <- "Log2 tranformation performed")
        msg_log <- paste0(msg_plot, "\nBefore transformation range: ", paste(range(na.omit(array_dat)),collapse =" to "))
        array_dat <- log2(array_dat)
        (msg_log <- paste0(msg_log, "\nAfter log2 transformation: ", paste(range(na.omit(array_dat)),collapse =" to ")))
    }
    
    ###############################
    ### remove NaN rows
    ###############################
    ## inspect NAN and remove them by rows
    temp <- which(is.na(array_dat) == TRUE)
    old_rows <- nrow(array_dat)
    array_dat <- data.frame(na.omit(array_dat), check.names=F) # must turn off the check names, keep special char in col names
    if(length(temp) > 0){
        rm_rows <- old_rows - nrow(array_dat)
        rm_p <- round(rm_rows/old_rows*100, 2)
        msg1 <- paste0("Expression data contain NAs, ", rm_rows, " out of ", old_rows, " removed (", rm_p,"%).")
        print(msg1)
    }else{
        msg1 <- "No probes removed due to missing value."
    }
    
    
    
    
    ###############################
    ### remove rows with no variation
    ###############################
    # input an array_dat, and remove rows with no variation (ie all readings are the same)
    # if not removed, there will be problem with PCA
    row_variance <- apply(array_dat, 1, var)
    index <- which(row_variance ==0)
    if (length(index) >0){
        array_dat <- array_dat[-index, ]
        msg_rm_pb <- paste0("\n", length(index), ' probes/probesets are removed due to 0 variation')
    }else{
        msg_rm_pb <- ""
    }
    
    
    
    ###############################
    ######### plots: optional
    ###############################
    # density plots of each sample
    main_title <- paste0('Sample Density: ', dataset, msg_col_rm,
                         '\n', ncol(array_dat), ' samples, ', nrow(array_dat), ' probes',
                         '\n', msg1,
                         '; ', msg_plot
    )
    png(filename= paste0(plot_pre,'density_exp.png'), width = width, height = height)
    plotDensities(array_dat, legend=F, main = main_title)
    dev.off()
    
    if(to_plot){
        ## plot heatmap and PCA plots
        plotFunctions(array_dat, plot_pre, dataset, array_design)
    }
    
    
    ###############################
    ######### test the factor and expression matrix
    ###############################
    (response_ls <- intersect(x_col, colnames(array_design)))
    (response_ls <- response_ls[-length(response_ls)]) # remove the last element, Sample
    if(model_test){msg_model <- modelTest(array_dat, array_design, response_ls)} else {msg_model <- ""}
    
    ###############################
    ######### output log: optional
    ###############################
    if(to_log){
        (msg2 <-printSampleCountsByFactor(array_design))
        ## save the process messages
        msg3 <- paste0("Expression data file: ", data_file, 
                       "\n", msg1,
                       msg_rm_pb,
                       "\n", msg_col_rm, 
                       "\nnumber of rows: ", nrow(array_dat),
                       "\nnumber of samples: ", ncol(array_dat),
                       "\n", msg2,
                       "\n", msg_log, msg_model)
        
    }
    ###############################
    ######### return dataframe
    ###############################
    returnlist <- list(array_dat, array_design, msg3)
    return(returnlist)
}


#---------------------------------------------------------------------------#
# PART 5: MAIN FUNCTION
#---------------------------------------------------------------------------#

exploreGSE <- function(df_process, index, method, model_test=F){
    ## process a dataset and output processed R objects (for data explore, spot problems, or for limma DEA)
    ## INPUT:
    ##   - df_process: A df contain info about the dataset and whether to subset or remove outliers
    ##   - index: which line of the df_process (to process 1 dataset at a time)
    ##   - method: choose one of 'explore_mode', 'analysis_mode'
    ## Notes:
    ##   - 
    
    if(class(method) != 'character'){
        stop(paste0('method is a ', class(method), "; must a character defined as one of 'explore_mode', 'analysis_mode'"))
    }
    if(!(method %in% c('explore_mode', 'analysis_mode'))){
        stop(paste0('method is ', method, "; must be one of 'explore_mode', 'analysis_mode'"))
    }
    
    ## check if df_process has all the required columns
    need_col <- c("eeName", "platform_id", "subset", "subset_by", "keep_subset", "outlier", "outlier_to_rm", "rm_notes")
    diff_col <- setdiff(need_col, colnames(df_process))
    if(length(diff_col) >0 ){
        stop(paste0("Require column(s): ", paste(diff_col, collapse = ', ')))
    }
    
    ######################################################################################
    ###################  Get dataset and platform info
    ######################################################################################
    model_info <- df_process[index, ]
    platform <- as.character(model_info$platform_id)
    dataset <- as.character(model_info$eeName)
    print(index)
    
    print(paste0("processing ", dataset, " and platform ", platform))
    
    ## if mode is analysis, only one entry of the dataset name is accepted, 
    ## explore_mode can allow multiple matching of the same dataset
    ## change to_do to to_do_ls
    if(method == 'analysis_mode' & !exists("to_do_ls")){
        tmp <- which(df_process$eeName == dataset)
        if(length(tmp) != 1){
            stop(paste0("There are ", length(tmp), " records of ", dataset, " in the input. Only 1 is allowed for ", method))
        }
    }
    if(method == 'analysis_mode' & exists("to_do_ls")){
        tmp <- as.character(df_process[to_do_ls, "eeName"])
        if(any(duplicated(tmp))){
            stop(paste0("There are duplicated records of ", dataset, " in the input. Only 1 is allowed for ", method))
        }
    }
    
    
    ######################################################################################
    ###################  set folder dir for each dataset
    # folder names are updated based on subsets or no subsets
    ######################################################################################
    ## output folders
    (dataset_folder <- paste0(output_folder, method, '/', dataset, '/'))  # subfolder of 'explore_mode' or 'analysis_mode'
    
    ## if it's in explore mode, folder are based on subsets or no_subset
    if(method == 'explore_mode'){
        if(tolower(model_info$subset) == 'yes'){
            (subset_factor_ls <- unlist(strsplit(as.character(model_info$subset_by), split ="; ")))
            (subset_keep_factor_ls <- unlist(strsplit(as.character(model_info$keep_subset), split ="; ")))
            ## reformat subfolder name
            tmp <- paste(paste(subset_factor_ls,subset_keep_factor_ls, sep = ""), collapse="_")
            tmp <- gsub("\'", "", tmp)
            tmp <- gsub("c\\(", "\\_", tmp)
            tmp <- gsub("\\)", "\\", tmp)
            tmp <- gsub("\\ ", "", tmp)
            (dataset_folder <- paste0(dataset_folder, 'subset/',tmp, '/'))
        } else {
            dataset_folder <- paste0(dataset_folder, 'no_subset/')
        }
    }
    
    print(paste0("output folder is ", dataset_folder))
    
    plot_folder <- paste0(dataset_folder, 'plots/')
    result_folder <- paste0(dataset_folder, 'results/')
    
    dir.create(plot_folder, showWarnings = FALSE, recursive=T)
    dir.create(result_folder, showWarnings = FALSE, recursive=T)
    
    plot_pre <- paste0(plot_folder, dataset, '_')
    result_pre <- paste0(result_folder, dataset, '_')
    
    ######################################################################################
    ###################  main: process the dataset and explore
    ######################################################################################
    array_design <- processDesign(dataset, s_file = f_design)
    ## return annotation and process msg
    x <- processPlatformAnnotation(platform, platform_folder)
    annotation <- x[[1]]
    msg_annotation <- as.character(x[[2]])
    
    ## return data matrix
    x <- processGemmaDataset(dataset, array_design, unfiltered_exp_folder, annotation, plot_pre, result_pre, model_test,
                             to_plot=T, to_log =T)
    array_dat <- x[[1]]
    array_design <- x[[2]] ## processed array design(with NaN samples removed)
    msg_data_exp <- as.character(x[[3]])
    
    ## save the correlation matrix
    sample_cor <- cor(array_dat)
    ## replace diag cells with NA value
    diag(sample_cor) <- NA
    sample_cor <- cbind(data.frame(Sample = row.names(sample_cor)), sample_cor)
    f_cor <- paste0(result_pre, "sample_correlations.tsv")
    write.table(sample_cor, file =f_cor, quote = F, sep ='\t',
                row.names =F)
    
    
    ## get the sanity check heatmap
    sanityCheckHeatmap(array_dat, f_gene_list, array_design, annotation, plot_pre, height, width, dataset = dataset, top_200 = T)
    
    ## remove outliers (optional), and save a df with outliers removed
    msg_rm_samples <- ""
    if(rm_outlier){
        if(tolower(as.character(model_info$outlier)) == "yes"){
            print("Removing outliers")
            to_rm_ls <- as.character(model_info$outlier_to_rm)
            rm_notes <- as.character(model_info$rm_notes)
            x <- rmSamples(to_rm_ls, array_dat, array_design, to_log =T, rm_notes=rm_notes)
            array_dat <- x[[1]]
            array_design <- x[[2]]
            msg_rm_samples <- as.character(x[[3]])
            outlier_rm_df<- x[[4]]
            if(!is.null(outlier_rm_df)){
                write.table(outlier_rm_df, file =paste0(result_pre, "outliers_removed.tsv"), quote = F, sep ='\t',
                            row.names =F)
                print(paste0('removed samples saved: ', paste0(result_pre, "outliers_removed.tsv")))
            }
            plot_pre <- paste0(plot_pre, "outlier_rm_")
            prefix <- " (outlier removed)"
            ## plot after outlier removed
            plotFunctions(array_dat, plot_pre, dataset, array_design, prefix)
            print(msg_rm_samples)
        }
    }
    
    
    ############################
    ## subset the data (optional)
    ############################
    msg_subset_samples <-""
    if(to_subset){
        if(tolower(as.character(model_info$subset)) == "yes"){
            print("Subset the samples")
            print(model_info)
            x <- subsetSamples(model_info, dataset, array_dat, array_design, to_log =T, model_test=model_test)
            array_dat <- x[[1]]
            array_design <- x[[2]]
            msg_subset_samples <- as.character(x[[3]])
            subset_rm_df<- x[[4]]
            if(!is.null(subset_rm_df)){
                write.table(subset_rm_df, file =paste0(result_pre, "subset_removed.tsv"), quote = F, sep ='\t',
                            row.names =F)
                #print(paste0('removed samples saved: ', subset_rm_df))
            }
            
            
            plot_pre <- paste0(plot_pre, "subset_")
            prefix <- " (subset)"
            ## plot after subset
            plotFunctions(array_dat, plot_pre, dataset, array_design, prefix)
            print(msg_subset_samples)
        }
    }
    
    
    x <- loopSampleOutlier(array_dat, array_design)
    msg_outlier <- x[[1]]
    df_outlier <- x[[2]]
    if(!is.null(df_outlier)){
        f_outlier <- paste0(result_pre, "outliers.tsv")
        write.table(df_outlier, file =f_outlier, quote = F, sep ='\t',
                    row.names =F)
    }
    
    ############################
    ## save process messages
    ############################
    msg <-paste0(Sys.Date(), '\n', msg_annotation, '\n', msg_data_exp, msg_rm_samples,msg_subset_samples, msg_outlier)
    sink(paste0(result_pre, "exp_process_log.txt"), type="output")
    writeLines(msg)
    sink()
    
    ############################
    ## set WT as the base level
    ############################
    array_design$Genotype <- factor(array_design$Genotype, levels = c('WT', setdiff(levels(array_design$Genotype), "WT")))
    
    ############################
    ## save r objects
    ############################
    save(array_design, annotation, array_dat, dataset, platform, plot_pre, result_pre, msg,
         file = paste0(result_pre, "objects.Rdata"))
}



#---------------------------------------------------------------------------#
# PART 5: LOOP (main function)
#---------------------------------------------------------------------------#
# ## set the loop index


exploreExpData <- function(){
    # 1. set a start_row to loop the experiment list, default start from first row
    if(!exists("start_row")){start_row <- 1}
    # 2. can select a list of experiments by row index (over-write default to_do) 
    if(!exists("to_do_ls")){to_do <- c(start_row:nrow(df_process))}
    # 3. pick a specific dataset by GSE id (over-write default to_do) 
    if(!exists("to_do_ls") & exists("dataset_todo")){
        to_do <- which(df_process$eeName %in% dataset_todo)}
    if(exists("to_do_ls")){
        to_do <- to_do_ls}
    for (i in to_do){
        exploreGSE(df_process, index=i, method=explore_method, model_test)
    }
}
