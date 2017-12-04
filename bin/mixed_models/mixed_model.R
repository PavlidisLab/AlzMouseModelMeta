#'created 2016-06-16
#'updated 2016-09-19 # updated for the new test model after model selection
#'
#'@see LME4 tutorial www.bodowinter.com/tutorial/bw_LME_tutorial.pdf
#'
#'aggregate expression (early or late phasese) based on dataset info_mixed model
#'go to limma R data to load expression and design (batch corrected, outlier removed)
#'match genotype and timeline to get each sample
#'
#'for each R.data, aggregate expression by gene (max median)
#'remove "" genes, rm multiple genes(see limma p value for genes)
#'for each design, change samples names in design and expression to dataset_sample names ( to combine with other datasets)
#'
#'loop aggregate studies from all examples from the phase (experssion), aggregate design
#'rm genes with > 1/3 NA or remove all NA (default)
#'disease_stage: WT, Disease
#'quantile normalize by sample
#'genes can be removed by affy_threshold and filter_method

# **********************************************************************************************#
# **********************************************************************************************#
#' ###################################
#' mixed model: random intercept model
#' ###################################
#' 
#' assumptions: fixed effects: Disease_stage and Gender are all the same for all studies (slope is the same for all studies).
#' only account for baseline differences in expression of a gene
#' for Gender, male is the baseline
#' 
#' model: (Expr ~ Disease_stage + Gender + (1|Study),  data=df, REML=FALSE)
#' fixed effect disease_stage and gender
#' random intercept: study (each study has it's own intercept)
#' 
#'compare to models (with or without the effect) ## must use the fixed effect after model selection
#' ~ 1+ gender + (1|study)
#' ~ 1+ gender + disease_stage + (1|study)
#' 
#' Interpretation for each gene
# disease change expression of [gene] by [estimate] +- [standard error], (model test p=[] ... chisq(df =[]) = .. ) in model?
# gender change ....

#' ###################################
#' mixed model: random slope model ## do not use this
#' ###################################
#' the effects of disease_stage is different for different studies. 
#' (gender is assumed the same for all study, but different intercept and different slope for each study). 
#' ie some mouse models
#' change the gene expression more( or less)) than other models in a different study.

# full_model = lmer(Expr ~ Disease_stage + Gender + (1+Disease_stage|Study), data=df, REML=FALSE)
# null_model = lmer(Expr ~ Gender + (1+Disease_stage|Study),  
#                   data=df, REML=FALSE)

# **********************************************************************************************#
# **********************************************************************************************#

#' check the residual, fitted value, and other plots to make sure the model fits well -not done
#' 
#' get the coefficients (per gene) and calculated CI and t, if CI includes 0 - no effect
#' get the p value by comparing to the null model
#' 
#' input:
#' a list of Robjects, dataset_info_mixed_model.tsv
#'@examples_1: from r files to MM result
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease ='AD'
# source('config/config_wrappers.R')
# source('mixed_models/mixed_model.R')
# md_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv')  #dataset label, phase label, order of the datasets etc.
# datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
# (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
# phase='early'
# out_dir <- paste0(disease_dir, 'mixed_model/')
# NA_filter='all'
# width = 1200
# height = 1200
# 
# ### just get the combined expression data
# x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
#                        width = width, height = height)
# 
# 
# ### get the expression data and MM results
# mixedModelAll(rfile_ls= rfile_ls, md_info=md_info, phase =phase, out_dir =out_dir, 
#               NA_filter=NA_filter, width = width, height = height, model='random_intercept')


#' @examples_2: from expression data to MM result
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease ='AD'
# source('config/config_wrappers.R')
# source('mixed_models/mixed_model.R')
# 
# phase='early'
# out_dir <- paste0(disease_dir, 'mixed_model/')
# (exprdata <- paste0(out_dir, phase, '/expression.Rdata'))
# 
# ### get the expression data and MM results
# mixedModelAll(phase =phase, out_dir =out_dir, 
#                           exprdata = exprdata)

#' @examples_3: from result data(partially done) to MM result
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease ='AD'
# source('config/config_wrappers.R')
# source('mixed_models/mixed_model.R')
# 
# phase='early'
# out_dir <- paste0(disease_dir, 'mixed_model/')
# (previous_result<- paste0(out_dir, phase, '/mixed_model_results.Rdata'))
# 
# ### get the expression data and MM results
# mixedModelAll(phase =phase, out_dir =out_dir, start_row = 123,
#               rm_genes = "",
#               previous_result = previous_result)



#----------------------------------------------------#
## pre req
#----------------------------------------------------#
library(lme4)
library(HelperFunctions)
library(preprocessCore)  # quantile normalization
source("helper_functions.R")

#----------------------------------------------------#
## functions
#----------------------------------------------------#

aggregateGeneExp <- function(array_dat){
    # input a array dat with expression data and probenames and genesymbols
    # return a array dat with genesymbol match to 1 most expressed probe (max median)
    # return a df with genesymbol and probes that represent the gene
    df <- array_dat
    rownames(df) <- array_dat$ProbeName
    
    ###############################
    #  get defaults and sanity check
    ###############################
    ## check the format of the input table
    if(!all(c("ProbeName", "GeneSymbols") %in% colnames(df))){
        stop(paste0("Input array matrix need to have columns 'ProbeName', 'GeneSymbols'"))}
    
    
    ###############################
    #  remove probes that don't map to any gene symbols
    ###############################
    df <- excludeMatch(df, "GeneSymbols", "")
    (msg <- paste0('Mapped probesets: ', nrow(df), "; "))
    
    
    ###############################
    #  rm probes that mapped to multiple genes
    ###############################
    ## rm probes that mapped to multiple genes (these probes are not specific to 1 gene, thus not useful)
    ## get index for the insulin genes (keep these probes even they mapped to insulin 1 and 2)
    (index_ins <- grep('Ins.\\|', df$GeneSymbols))  # index for insulin
    index <- setdiff(grep('\\|', df$GeneSymbols), index_ins)
    (msg <- paste0(msg, 'probes mapped to multiple genes to be removed (except for insulin): ', length(index), "; "))
    df <- df[-index, ] %>%droplevels()
    
    # change the "Ins2|Ins1" to "Ins1|Ins2" to make it consistant across datasets
    (index <- which(df$GeneSymbols == "Ins2|Ins1"))
    if(length(index >0 )){df[index, "GeneSymbols"] <- "Ins1|Ins2"}
    
    
    ###############################
    ## aggregate expression value
    ## by ....maxMedian
    ###############################
    
    df$median <- apply(df[, setdiff(colnames(df), c('ProbeName', 'GeneSymbols'))], 1, function(x)median(x, na.rm = T))
    
    tmp <- aggregate(median ~ GeneSymbols, data = df, max) ## get the max median
    probe_df <- noWarnings(left_join(tmp, df[, c('ProbeName', 'GeneSymbols','median')])) ## get the matching probenames
    
    ## get the new matrix based on the probes
    df <- df[probe_df$ProbeName, ]%>% droplevels()
    
    if(nrow(tmp) != nrow(df)){
        df2 <- df
        df2$avg <- apply(df2[, setdiff(colnames(df), c('ProbeName', 'GeneSymbols', 'median'))], 1, 
                         function(x) mean(x, na.rm = T))
        tmp <- aggregate(avg ~ GeneSymbols, data = df2, max) ## get the max avg
        probe_df <- noWarnings(left_join(tmp, df2[, c('ProbeName', 'GeneSymbols','avg')])) ## get the matching probenames
        ## get the new matrix based on the probes
        df <- df[probe_df$ProbeName, ]%>% droplevels()
    }
    
    ## make sure no probes mapped to the same gene has the same median
    
    
    ## save the probenames and genesymbols as another df
    df_genes <- df[,  c('ProbeName', 'GeneSymbols')] %>% droplevels()
    rownames(df_genes) <- NULL
    
    ## remove median and probeName symbols, use genesybols as rownames
    rownames(df) <- df$GeneSymbols
    df <- df[, setdiff(colnames(df), c('ProbeName', 'GeneSymbols','median'))]
    
    returnlist <- list(df, df_genes)
    return(returnlist)
}


processExpDesign <- function(rfile, disease_genotype_ls){
    ## input an robject file and a list of disease genotypes
    ## output processed array_dat and array_design
    
    #-----------------#
    ## read r object, get the samples (exp and design)
    #-----------------#
    
    ## load r data
    load(rfile)
    
    ## filter design and expression for the needed samples only, assign samples disease or not
    array_design$Disease_stage <- 'WT'  # categorized mouse model (disease); control (WT)
    array_design[which(array_design[, 'Genotype'] %in% disease_genotype_ls), 'Disease_stage'] <- 'Disease'
    
    
    ## define study
    (array_design$Study <- as.factor(gsub("\\.[0-9]", "", array_design$Dataset)))  ## replace the dataset name if there's a split, e.g. GSE1234.1 = GSE1234
    
    ###########
    ## define age as numerical (unit: m), convert months and weeks to numerical months
    ##########
    array_design$Age <- as.character(array_design$Timepoint)
    (index <- grep('_weeks', array_design$Age))
    if(length(index)>0){
        for(i in index){
            (tmp <- unlist(strsplit(array_design$Age[i], split = '-|_')))
            if(length(tmp)==3){  ## if the format is 3-10_weeks
                weeks <- mean(c(as.numeric(tmp[1]),as.numeric(tmp[2])))
            }else{weeks <- as.numeric(tmp[1])}
            (array_design$Age[i] <- weeks/4)  ## from week to month
        }}
    (index <- grep('_months', array_design$Age))
    if(length(index)>0){
        for(i in index){
            (tmp <- unlist(strsplit(array_design$Age[i], split = '-|_')))
            if(length(tmp)==3){  ## if the format is 3-10_weeks
                array_design$Age[i] <- mean(c(as.numeric(tmp[1]),as.numeric(tmp[2])))
            }else{array_design$Age[i] <- as.numeric(tmp[1])}
        }}
    array_design$Age <- as.numeric(array_design$Age)
    
    ###########
    ## or dataset is subset by batch and different genotypes compared to different WT
    #(array_design$Study <- as.factor(paste0(array_design$Dataset,'_', array_design$Timepoint)))  ## do not seperate subsets of dataset, because WT are different for comparison, 
    ##########
    
    need_col <- c('Study','Sample','Disease_stage','Genotype','Timepoint', 'Gender','Dataset', 'Age')
    array_design <- array_design[, need_col]%>% droplevels()
    
    ## get the needed samples by genotype
    array_design <- array_design[which(array_design$Genotype %in% c(disease_genotype_ls, 'WT')), ]%>% droplevels()
    
    ## get the correponding exp matrix
    array_dat <- array_dat[, as.character(array_design$Sample)]
    
    ## add "ProbeName", "GeneSymbols" to matrix
    array_dat$ProbeName <- rownames(array_dat)
    array_dat <- noWarnings(left_join(array_dat, annotation[, c('ProbeName', 'GeneSymbols')]))
    
    ## save the matrix before aggregation
    ## update sample names
    # mk this as a list
    array_dat_bf_aggregation <- array_dat
    rownames(array_dat_bf_aggregation) <- array_dat_bf_aggregation[, 'ProbeName']
    colnames(array_dat_bf_aggregation) <- paste0(levels(array_design$Study), '_', colnames(array_dat_bf_aggregation))
    list_name <- paste0(levels(array_design$Study), '_', gsub('onths|ays|eeks', '', levels(array_design$Timepoint)))
    
    ## aggregate the probes to gene expression, gene symbol as rownames
    x <- aggregateGeneExp(array_dat)
    array_dat <- x[[1]]
    df_gene_match <- x[[2]] ## this is the info of which probes are selected for the gene
    df_gene_match$Dataset <- levels(array_design$Dataset)
    
    ## update sample names 
    colnames(array_dat) <- paste0(array_design$Study, '_', colnames(array_dat))
    array_design$Sample <- paste0(array_design$Study, '_', array_design$Sample)
    rownames(array_design) <- array_design$Sample
    
    
    returnlist <- list(array_dat, array_design, array_dat_bf_aggregation, df_gene_match, list_name)
    return(returnlist)
}



combineAllStudies <- function(rfile_ls, md_info, phase, out_dir, NA_filter=c('all', '0.3', 'none')){
    #' input a list of r files and phase to select the correponding r files
    #' output a array dat with all samples and expression are genes (represented by the max median probes), and NA removed by NA_filter method
    #' output a array design correponding to all samples
    #' for the mixed model 
    #' output: save a r object
    #' low_exp_rm: remove genes with expression lower than affy_threshold, by the filter_method
    
    
    #' rfile_ls: list of limma R objects
    #' md_info: mixed model dataset info (similar to dataset info)
    #' phase: str to specify which phase (fto select correponding samples)
    #' out_dir: save r objects
    #' NA_filter: all: filter all NAs, if a probe contains at least 1 NA sample, removed (default)
    #'              0.3: if a probe has > 1/3 NAs in all samples, filter out
    #'              none: no NA filtered
    
    dir.create(out_dir, recursive = T, showWarnings = F)
    
    ## get corresponding r objects
    info_df <- getRObjectLabels(rfile_ls, md_info)
    
    ## get the phase
    info_df <- filterContain(info_df, "Phase", phase)
    
    genotype_ls <- info_df$Genotype
    file_ls <- unique(info_df$File)  # no duplicate R files
    
    
    #--------------#
    ## Loop for each r object and get all the expression gene value for all samples from all studies of the phase
    #--------------#
    all_array_dat <- NULL
    all_array_design <- NULL
    all_array_dat_bf_aggregation <- vector("list")  ## all probes before aggregation
    all_df_gene_match <- NULL  ## which gene match to which probe
    for( i in 1:length(file_ls)){
        (rfile <- file_ls[i])
        print(i)
        print(rfile)
        (indi_info <- filterContain(info_df, 'File', rfile))  # may contain multiple entries of the same file, but different disease genotypes share the same WT
        (disease_genotype_ls <- as.character(indi_info$Genotype))  # for disease model only (no WT)
        x <- processExpDesign(rfile, disease_genotype_ls)
        array_dat <- x[[1]]
        array_design <- x[[2]]
        array_dat_bf_aggregation <- x[[3]] ## return as a dataframe probes before aggregation
        df_gene_match <- x[[4]]
        list_name <- x[[5]]  ## name for the array_dat_bf_aggregation
        
        
        ## for aggregated genes
        if(is.null(all_array_dat)){
            all_array_dat <- array_dat
        }else{
            all_array_dat <- mapBindAllColumns(all_array_dat, array_dat)
        }
        
        ## for all probes
        all_array_dat_bf_aggregation[[i]] <- array_dat_bf_aggregation
        names(all_array_dat_bf_aggregation)[[i]] <- list_name
        
        ## which gene to which probe
        if(is.null(all_df_gene_match)){
            all_df_gene_match <- df_gene_match
        }else{
            all_df_gene_match <- rbind(all_df_gene_match, df_gene_match)
        }
        
        ## combine designs
        if(is.null(all_array_design)){
            all_array_design <- array_design
        }else{
            all_array_design <- rbindAllColumns(all_array_design, array_design)
        }
        
    }
    
    ##setdiff(colnames(all_array_dat), all_array_design$Sample)
    
    msg <- paste0(Sys.Date(), ': combined study expressions. \nPHASE: ', phase, 
                  '\nINPUT: ', length(file_ls),  ' Rdata files; \n    ',
                  paste0(file_ls, collapse = '\n    '),
                  '\n\nTOTAL SAMPLES: ', ncol(all_array_dat), 
                  '; from ', length(levels(all_array_design$Study)), ' studies: ', paste0(levels(all_array_design$Study), collapse=','))
    
    
    #--------------#
    ## Filter NAs...
    # no filter
    # filter all NA
    # filter when > 1/3 are NAs
    #--------------#
    NA_filter <- NA_filter[1]
    if(!NA_filter %in% c('all', '0.3', 'none')){
        stop(paste0('NA_filter is ', NA_filter, "; must select one of 'all', '0.3', 'none'"))
    }
    
    if(NA_filter=='all'){
        array_dat <- na.omit(all_array_dat)  ## filter all
        msg <- paste0(msg, 
                      '\nSTAT: \n NA filter method: remove all NA',
                      '\n   Original genes: ', nrow(all_array_dat),
                      '\n   NA genes removed: ', nrow(all_array_dat) - nrow(array_dat))
    }
    
    #     ## if filter when > 1/3 are NAs
    #     if(NA_filter== '0.3'){
    #         NA_count <- apply(is.na(all_array_dat),1,sum)
    #         keep_index <- which(NA_count <= ncol(all_array_dat)/3)
    #         array_dat <- all_array_dat[keep_index, ]%>% droplevels()
    #         msg <- paste0(msg, 
    #                       '\nSTAT: \n NA filter method: remove genes with 1/3 NAs',
    #                       '\n   Original genes: ', nrow(all_array_dat),
    #                       '\n   NA genes removed: ', nrow(all_array_dat) - nrow(array_dat),
    #                       '\n   Final gene count: ', nrow(array_dat))
    #     }
    
    ## if filter when > 1/3 the studies don't have the gene expression
    if(NA_filter== '0.3'){
        df_all <- NULL
        total_studies <- unique(as.character(all_array_design$Study))
        ## loop to find which study has NA for the gene
        for(study in total_studies){
            print(study)
            study_s <- filterContain(all_array_design, 'Study', study)$Sample
            NA_study <- apply(is.na(all_array_dat[, study_s]),1,max)  # see which gene has NA
            df <- as.data.frame(NA_study)
            colnames(df) <- study
            
            if(is.null(df_all)){
                df_all <- df
            }else{
                df_all <- mapBindAllColumns(df_all, df)}}
        # end of loop
        
        NA_count <- apply(df_all,1,sum)
        keep_index <- which(NA_count < length(total_studies)/3)
        array_dat <- all_array_dat[keep_index, ]%>% droplevels()
        msg <- paste0(msg, 
                      '\nSTAT: \n NA filter method: remove genes presented in less than 2/3 of the studies',
                      '\n   Original genes: ', nrow(all_array_dat),
                      '\n   NA genes removed: ', nrow(all_array_dat) - nrow(array_dat))
        
    }
    
    if(NA_filter== 'none'){
        array_dat <- all_array_dat
        msg <- paste0(msg, 
                      '\nSTAT: \n NA filter method: No NA filtered')
    }
    
    
    
    ## returns
    returnlist <- list(array_dat, all_array_design, msg, all_array_dat_bf_aggregation, all_df_gene_match)
    return(returnlist)
}


#*********
# Mixed model (random intercept) for each gene
#*********

mixedModelStudy <- function(array_dat, array_design, i, 
                            model=c('random_intercept', 'random_slope'),
 #                           fixed_effect, 
 #                           random_effect = c('Study'),
                            to_plot=F,
                            plot_dir ='',
                            full_report = F,
                            REML=FALSE,
                            estimateCI= T,
                            disease_stage_only =T,
                            study_color =NULL){
    #' array_dat (full array, with gene symbol as row names)
    #' array_design(full)
    #' i: index of a gene name from array_dat, or a gene name
    #' output a df with the model stats and CI, and pvalue(compare null model and model with disease stage, anova test)
    #' model: (Expr ~ Disease_stage + Gender + (1|Study),  data=df, REML=REML)
    #' fixed effect disease_stage and gender
    #' @param full_report: report all the coefs for disease and gender
    #' @param to_plot, plot the residuals, histogram, box plots etc 
    #' model: choose random_intercept (default) or random_slope (see assumptions)
    #' random intercept: study (each study has it's own intercept)
    
    ## plot name:
    plot_dir <- paste0(plot_dir, '/diagnostic_plots/')
    plot_pre <- paste0(plot_dir, i, '_')
    
    ##pick 1 model
    #model <- model[1] ## random_intercept is the default
    
    ## preprocess to get the expression df for lmer()
    #print(i)
    df <- as.data.frame(t(array_dat[i, ])) # expression of the genes from all samples
    
    i <- colnames(df)

    #print(i)
    colnames(df) <- 'Expr'
    df$Sample <- row.names(df)
    df <- noWarnings(left_join(df, array_design))
    row.names(df) <- df$Sample
    ## update level order
    df$Disease_stage <- factor(df$Disease_stage, levels = c('WT', 'Disease'))
    if('Gender' %in% colnames(df)){ ## set gender M as baseline
        if(length(levels(df$Gender)) ==2){
            df$Gender <- factor(df$Gender, levels = c('M', 'F'))
            gender_col <- 'GenderF'
        }else{gender_col <- NA}
    }
    
    
    ## keep non NA expressions
    df <- df[which(is.na(df$Expr)==F),] %>% droplevels()
    
    ##******************************##
    ## random intercept model (slope do not change for random variables)
    ## if only 1 gender in the design
    ##******************************##
    
#     #########
#     ## check for confound fixed effects(to study)
#     #########
#     ## if each study only contains 1 gender 
#     (tmp <- rmDup(df[,c('Gender','Study')]))
#     if(!any(as.numeric(table(tmp[,'Study']))>1)){## remove gender if T (F when at least 1 study contains both gender), also rm when gender =NA
#         fixed_effect <-setdiff(fixed_effect,'Gender')
#     }   
#     
#     ## if age group only belongs 1 study (i.e. each study has distinct 1 age group)
#     (tmp <- rmDup(df[,c('Age','Study')]))
#     if(!any(as.numeric(table(tmp[,'Age']))>1)){## rm age if T (F when at least 1 age belongs to multiple studies)
#         fixed_effect <-setdiff(fixed_effect,'Age')
#     }  
#     
#     
#     
#     
#     if('Gender' %in% fixed_effect){  ## add gender as fixed
#         model_formular <-  paste0('Gender + ', model_formular)
#         assign(paste0('model_formular_', counter),model_formular)
#         test_model_effects <- c(test_model_effects, model_formular)
#         counter = counter +1
#     }
#     if('Age' %in% fixed_effect){  ## add age as fixed
#         model_formular <-  paste0('Age + ', model_formular)
#         assign(paste0('model_formular_', counter),model_formular)
#         test_model_effects <- c(test_model_effects, model_formular)
#         counter = counter +1
#     }
#     
#     ## add disease
#     model_formular <-  paste0('Disease_stage + ', model_formular)
#     assign(paste0('model_formular_', counter),model_formular)
#     test_model_effects <- c(test_model_effects, model_formular)
#     counter = counter +1
#     
#     ## add disease age interaction
#     if('Age' %in% fixed_effect){  ## add age as fixed
#         model_formular <-  paste0('Age * Disease_stage + ', model_formular)
#         assign(paste0('model_formular_', counter),model_formular)
#         test_model_effects <- c(test_model_effects, model_formular)
#         counter = counter +1
#         effect_ls <- c(effect_ls, 'disease and age interaction')
#     }
#     
#     
    
    
    
    
    
    if(model == 'random_intercept'){
#         if(any(length(levels(df$Gender))==1, length(levels(df$Gender))==F)){  ## if only 1 gender or gender is NA
#             full_model = lmer(Expr ~ Disease_stage + Age + (1|Study),  
#                               data=df, REML=REML)
#             null_model = lmer(Expr ~ Age + (1|Study),  
#                               data=df, REML=REML)
#         }else{
#             full_model = lmer(Expr ~ Disease_stage + Gender + (1|Study),  
#                               data=df, REML=REML)
#             null_model = lmer(Expr ~ Gender + (1|Study),  
#                               data=df, REML=REML)
#         }
        if(disease_stage_only){
            full_model = lmer(Expr ~ Disease_stage + (1|Study),  
                              data=df, REML=REML)
            null_model = lmer(Expr ~ (1|Study),  
                              data=df, REML=REML)
        }else{
            full_model = lmer(Expr ~ Disease_stage * Age + (1|Study),  
                              data=df, REML=REML)
            null_model = lmer(Expr ~ Age + (1|Study),  
                              data=df, REML=REML)}
        }

#     } else{
#         ##******************************##
#         ## random slope model (slope change for random variables)
#         ## if only 1 gender in the subset df
#         ##******************************##
#         if(any(length(levels(df$Gender))==1, length(levels(df$Gender))==F)){
#             full_model = lmer(Expr ~ Disease_stage + (1+Disease_stage|Study),  
#                               data=df, REML=REML)
#             null_model = lmer(Expr ~ (1+Disease_stage|Study),  
#                               data=df, REML=REML)
#         }else{
#             full_model = lmer(Expr ~ Disease_stage + Gender + (1+Disease_stage|Study), data=df, REML=REML)
#             null_model = lmer(Expr ~ Gender + (1+Disease_stage|Study),  
#                               data=df, REML=REML)
#         }
#     }
    
    ## get the pvalue by comparing to the null model
    (anova_r <- anova(full_model,null_model))
    (pvalue <- anova_r[, "Pr(>Chisq)"][2])
    (chisq <- anova_r[, "Chisq"][2])
    (chi_df <- anova_r[, "Chi Df"][2])
    
    ## extract information:
    ## get estimate, se, t of the disease_stage
    (df_t <- data.frame(geneSymbol =i, t(coef(summary(full_model))['Disease_stageDisease',])))
    
    
    
    ##################
    ## estimate 95% CI
    ##################
    if(estimateCI){
        ## get the CI, and whether CI overlap with 0
        (tmp <- try(noWarnings(confint.merMod(full_model, level = 0.95)), silent = T))
        if(inherits(tmp, "try-error")){        
            #error handling code, Error in zetafun(np, ns) : profiling detected new, lower deviance
            ## cant get CI, then assign NA
            ci <- data.frame(CI_2.5=NA, CI_97.5=NA)
        }else{
            (ci <- as.data.frame(tmp)['Disease_stageDisease', ]) ## for slope estimate CI
            colnames(ci) <- c('CI_2.5', 'CI_97.5')
        }
    }else{
        ci <- data.frame(CI_2.5=NA, CI_97.5=NA)
    }
    
    
    if(full_report){ # report all coefficients
        ## get coefficients for each study
        x <- as.data.frame(t(coefficients(full_model)$Study))
        # get intercepts
        (all_intercepts <- as.data.frame(x['(Intercept)', ]))
        colnames(all_intercepts) <- paste0(colnames(all_intercepts), '_intercept')
        # get fixed effects: (fixed intercept model)
        df_effect <- as.data.frame(x["Disease_stageDisease", ])
        colnames(df_effect) <- paste0(colnames(df_effect), '_FE_Disease')
        if(is.null(gender_col)){ ## if gender is not in the model
            df_effect <-cbind(df_effect, data.frame(FE_Gender = NA, Gender = levels(factor(df$Gender))))
        }else{ ## if gender is in the model
            df_effect <-cbind(df_effect, data.frame(FE_Gender = x["GenderF", 1], Gender = gsub('Gender', '', gender_col)))
        }
        ## report all pvalues and coefficients
        df_t <- cbind(df_t, ci, pvalue, chisq, chi_df,all_intercepts, df_effect, data_frame(ME_model = model))
    }else{
        ## only report p values
        df_t <- cbind(df_t, ci, pvalue, chisq, chi_df, data_frame(ME_model = model))
    }
    
    ## if including plots (residual, box, histograms etc)
    if(to_plot){
      ## plot name:
        dir.create(plot_dir,showWarnings = F, recursive = T)
        ##################
        ## output model information of the gene
        ##################
        dir.create(plot_dir, showWarnings = F)
        f_out = paste0(plot_pre, 'model_summary.txt')
        sink(f_out, type="output")
        cat(paste0(i, '\n######### full model ',Sys.time(),'\n'))
        print(full_model)
        cat(paste0('\n\n######### full model ','\n'))
        print(summary(full_model))
        cat(paste0('\n\n######### NULL model ','\n'))
        print(summary(null_model))
        cat(paste0('\n\n######### ANOVA test ','\n'))
        print(anova_r)
        sink()
            
      ######
      ## diagnostic plots
      ######
        plot_name_dia <- paste0(plot_pre, i, '_diagnostic_plots.png')
      
        png(filename = plot_name_dia, width =800, height =800)
        #graphics.off()
        par(mfrow=c(2,2))
        #1 QQ plot
        qqnorm(residuals(full_model), main =paste0(i, ': Normal Q-Q plot'))
        qqline(residuals(full_model), col = 2,lwd=2,lty=2)
        #2 residual vs fitted
        plot(fitted(full_model),residuals(full_model), main = 'Residuals vs. Fitted',xlab = 'Fitted', ylab = 'Residuals')
        abline(h=0,col = 2,lwd=2,lty=2)
        #3 check normality of residuals
        hist(residuals(full_model), main =paste0('Residuals'), xlab = 'Residuals')
        #4 sample expression
        boxplot(Expr ~ Disease_stage, data=df, ylab = 'Expression', main ='Gene expression')
        dev.off()
        
        ########
        ## plot the predicted and the original expr
        ########
        plot_predicted <- paste0(plot_pre, i, '_predicted_intercept_slope.png')
        
        df_plt <- mutate(df, predicted = predict(full_model))
        p <- ggplot(df_plt, aes(Disease_stage, predicted, group=Study, color = Study)) + 
            geom_line() +
            theme_bw() +
            geom_point(aes(Disease_stage, Expr, group = Study)) +
            ggtitle(paste0(i, ": Predicted value vs. expression"))
        
        
        if(!is.null(study_color)){
            p <- p +
            scale_colour_manual(values = study_color) +  # customize color
            scale_fill_manual(values = study_color)}
        
        ggsave(filename = plot_predicted, plot=p, width = 6, height =8, units = "in")
        
        ####
        ## corrected for the intercepts
        ####
        plot_cor <- paste0(plot_pre, '_corrected_expression_', i, '.png')
        
        x <- as.data.frame(t(coefficients(full_model)$Study))
        # get intercepts and subtract from expression
        (cor_intercept <- as.numeric(x['(Intercept)', ]))
        names(cor_intercept) <- colnames(x)
        cor_value <- mapValue(as.character(df_plt$Study), cor_intercept)
        df_plt$Expr_corrected <- df_plt$Expr-cor_value
        
        ## the line is the mean
        p <- ggplot(df_plt, aes(Disease_stage, Expr_corrected, group = Study, color= Study)) + 
          theme_bw() +
          geom_point(alpha=0.6) +
          stat_summary(fun.y = mean, geom="line")+
          ggtitle(paste0(i))+
            ylab("")+
            xlab("")
        if(!is.null(study_color)){
            p <- p +
                scale_colour_manual(values = study_color) +  # customize color
                scale_fill_manual(values = study_color)}
        ggsave(filename = plot_cor, plot=p, width = 6, height =8, units = "in")
        
    }
    
    
    rownames(df_t) <- i
    
    ## if need to check residual etc TO DO
    
    
    return(df_t)
}


#*********
# combine all samples, and process design
#*********
prepareAllSamples <- function(rfile_ls, md_info, phase, out_dir, NA_filter=c('all', '0.3', 'none'),
                              width = 1200, height = 1200,
                              affy_threshold = 6, 
                              filter_method =c("median", "max", "min",  "mean"),
                              low_exp_rm=F){
    ## input a list of r files, only select the phase rfiles
    ## output: combine all samples expression, and a design file
    ## plot before quantile normalization and after normalization
    ## output Rdata with non-QN expression matrix, QN expression matrix, and design
    
    (plot_dir <- paste0(out_dir, '/',phase,'/'))
    dir.create(plot_dir, recursive = T, showWarnings = F)
    (plot_pre <- paste0(plot_dir, 'plot_', Sys.Date(), "_"))
    
    ## get the combined array dat and array design
    cat('COMBINE ALL STUDIES')
    x <- combineAllStudies(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter)
    array_dat_not_qn <- x[[1]] # array not quantile normalized
    array_design <- x[[2]]
    msg <- x[[3]]
    all_array_dat_bf_aggregation <- x[[4]]
    all_df_gene_match <- x[[5]]
    
    
    ## density plots before normalization:
    main_title <- paste0('Sample Density before normalization: ', 
                         '\n', ncol(array_dat_not_qn), ' samples, ', nrow(array_dat_not_qn), ' genes')
    png(filename= paste0(plot_pre,'density_before_normalization.png'), width = width, height = height)
    main_title=''
    plotDensities(na.omit(array_dat_not_qn), legend=F, main = main_title)
    dev.off()
    #return(array_dat_not_qn)
    
    cat(toupper('\n#############quantile normalization#############\n'))
    ## quantile normalization:
    (samples <- colnames(array_dat_not_qn)) ## need to get the sample names first
    genes <- rownames(array_dat_not_qn)
    array_dat <- as.data.frame(normalize.quantiles(as.matrix(array_dat_not_qn))) # col names are removed ## must double check this
    colnames(array_dat) <- samples
    rownames(array_dat) <- genes
    
    
    ## if remove lowly expressed genes from the array expression
    if(low_exp_rm){
        cat(paste0('\nremove lowly expressed genes lower than ', affy_threshold))
        filter_method <- filter_method[1]
        assign("index_exp", eval(parse(text = paste0("apply(as.matrix(array_dat),1,",filter_method, ")"))))
        
        index_exp <- apply(as.matrix(array_dat),1,function(x) median(x, na.rm = T))
        ## remove low exp in array dat
        array_dat <- array_dat[which(index_exp > affy_threshold), ]   #(for Affy)
        
        msg <- paste0(msg, 
                      '\nREMOVE lowly expressed genes, threshold = ', affy_threshold, '; filter method = ', filter_method,
                      '\n   Genes removed: ', length(which(index_exp <= affy_threshold)))
    }
    
    
    cat(toupper('\n#############making plots#############\n'))
    ## QC plots before scale
    plotFunctions(na.omit(array_dat), plot_pre, dataset="", array_design=array_design, prefix="", simple_heatmap = F)
    
    
    ## QC plots after row mean subtract
    ## subtract row mean from each row, don't use scale function, it's by column mean
    #    array_dat_mean <- data.frame(t(apply(array_dat, 1, function(y) y - mean(y))))
    ## same as array_dat_mean <- data.frame(t(scale(t(array_dat), scale=F)))
    #    plot_pre_scale <- paste0(plot_pre, 'scaled_')
    #     plotFunctions(array_dat= array_dat_mean,plot_pre_scale , dataset="", 
    #                   array_design=array_design, prefix="scaled", simple_heatmap = F)
    
    save(array_dat_not_qn, array_dat, array_design, phase, msg, plot_dir,all_array_dat_bf_aggregation, all_df_gene_match,
         file = paste0(plot_dir, 'expression.Rdata'))
    
    
    msg <- paste0(msg,
                  '\n   Final gene count: ', nrow(array_dat))
    cat(paste0('\n',msg))
    
    cat(paste0('\nSaved expression Rdata: ', plot_dir, 'expression.Rdata'))
    ## return
    returnlist <- list(array_dat_not_qn, array_dat, array_design, phase, msg, plot_dir)
    return(returnlist)
}

#----------------------------------------------------#
## main functions
#----------------------------------------------------#

#'quantile normalize by sample
#' mixed model
#'model expression ~ genotype + gender + (1|study) + error
#'
#'compare to models (with or without the effect)

#' ~ 1+ gender + (1|study)
#' ~ 1+ gender + disease_stage + (1|study)
#' 
#' check the residual, fitted value, and other plots to make sure the model fits well
#' 
#' get the coefficients (per gene) and calculated CI and t, if CI includes 0 - no effect
#'     
#' 
#' check the residual, fitted value, and other plots to make sure the model fits well -not done yet
#' 
#' get the coefficients (per gene) and calculated CI and t, if CI includes 0 - no effect
#' 
#' @param exprdata : the expression.Rdata that contains array dat, msg, and array design
#' @param previous_result: mixed_model_results.Rdata, continue from the previous run, specify start rows

mixedModelAll <- function(rfile_ls=NULL, md_info=NULL, phase, out_dir, NA_filter=c('all', '0.3', 'none'),
                          width = 1200, height = 1200, 
                          rm_genes = "",
                          start_row =1, 
                          exprdata = NULL,
                          previous_result=NULL,
                          model=c('random_intercept', 'random_slope'),
                          full_report = F,
                          to_plot= F,
                          tmp_result_output=T,
                          REML=FALSE,
                          affy_threshold = 6, 
                          filter_method =c("median", "max", "min",  "mean"),
                          low_exp_rm=F,
                          estimateCI= T,
                          disease_stage_only =T){
    print(Sys.Date())
    
    ## dir to save rdata
    (result_dir <- paste0(out_dir, '/', phase,'/'))
    dir.create(result_dir, recursive = T, showWarnings = F)
    #--------------#
    ## prep expression data if no rdata as input
    #--------------#
    if(!is.null(rfile_ls)){
        x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
                               width = width, height = height,
                               affy_threshold = affy_threshold, 
                               filter_method =filter_method,
                               low_exp_rm=low_exp_rm)
        array_dat <- x[[2]]
        array_design <- x[[3]]
        msg <- x[[4]]
    }else if(!is.null(exprdata)){
        load(exprdata)
        #start_row = 1
    }
    if(!is.null(previous_result)){load(previous_result)}
    save(array_dat,array_design, msg, file=paste0(result_dir, "mixed_model_before_results.Rdata"))
    
    if(tmp_result_output){
        tmp_result_output_f <- paste0(result_dir, 'tmp_mm_results.tsv')
    }
    
    
    ## define color for each study for ggplots
    study_color <- levels(array_design$Study)
    tmp <- ggColorHue(length(study_color))
    (names(tmp) <- study_color)
    study_color <- tmp
    
    
    #--------------#
    ## mixed model per gene
    #--------------#
    ## for loop here for each gene for the mixed model
    if(!('df_all' %in% ls())){df_all <- NULL} ## if a previous result is loaded
    print(paste0('startrow: ', start_row))
    genes_index <- start_row:nrow(array_dat)
    genes_index <- setdiff(genes_index, which(row.names(array_dat) %in% rm_genes))
    
    print(paste0('total genes : ', length(genes_index)))
    #####
    
    for (i in genes_index){
        df_t <- mixedModelStudy(array_dat, array_design, i, model=model, 
                                plot_dir = plot_dir,
                                full_report = full_report,
                                to_plot= to_plot,
                                REML=REML,
                                estimateCI= estimateCI,
                                disease_stage_only =disease_stage_only,
                                study_color = study_color)
        
        if(tmp_result_output){
            if(i==1){ # get the column names
                msg_t <- gsub('\n', '\n# ', msg)
                writeTable(df=df_t,f_out = tmp_result_output_f, msg = paste0('# ', Sys.Date(), '\n#', msg))
            }else{ ## else append the results 
                writeTable(df=df_t, f_out=tmp_result_output_f, append=T, col.names =F)
            }
        }
        
        if(is.null(df_all)){
            df_all <- df_t
        }else{
            df_all <- rbindAllColumns(df_all, df_t)  ## can bind NAs
        }
    }
    
    ## remove duplicates
    df_all <- rmDup(df_all)
    
    ## log
    msg <- paste0(Sys.Date(), ": Model testing: ", model, '\n', msg) 
    if(all(rm_genes != "")){
        msg <- paste0(msg, '\nRemove ', length(rm_genes), ' genes in the mixed model due to ERROR, profiling detected new, lower deviance: \n   ',
                      paste0(rm_genes, collapse = ', '))
    }else{
        msg <- paste0(msg, '\n Remove 0 genes in the mixed model due to ERROR, profiling detected new, lower deviance')
    }
    
    
    # attach expression value
    df_all_expression <- mapBindAllColumns(df_all, array_dat)
    
    save(array_dat,array_design, df_all, msg, df_all_expression,file=paste0(result_dir, "mixed_model_results.Rdata"))
    writeTable(df=NULL, msg = msg, f_out = paste0(result_dir, 'mixed_model_results_log_', Sys.Date(), '.txt'))
    
    ##2016-09-20 temp
    model_formular = "Expr ~ Disease_stage * Age + (1|Study)"
    model_formular = "Expr ~ Disease_stage + (1|Study)"
    
    returnlist=list(array_dat,array_design, df_all, msg, model_formular)
    print(paste0(result_dir, "mixed_model_results.Rdata"))
    return(returnlist)
}