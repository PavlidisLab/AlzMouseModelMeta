#'created 2017-02-15

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
#' Fixed model using lm()
#' ###################################
#' 


#----------------------------------------------------#
## pre req
#----------------------------------------------------#

library(HelperFunctions)
source("helper_functions.R")
source('mixed_models/archive/mixed_model.R') ## to prepare matrix

#----------------------------------------------------#
## functions
#----------------------------------------------------#


#*********
# fixed model (random intercept) for each gene
#*********

fixedModelStudy <- function(array_dat, array_design, i, 
                            model=c('fixed_model'),
 #                           fixed_effect, 
 #                           random_effect = c('Study'),
                            to_plot=F,
                            plot_dir ='',
                            full_report = F,
                            # REML=FALSE,
                            #estimateCI= T,
                            # disease_stage_only =T,
                            study_color=NULL){
    #' array_dat (full array, with gene symbol as row names)
    #' array_design(full)
    #' i: index of a gene name from array_dat, or a gene name
    #' output a df with the model stats and CI, and pvalue(compare null model and model with disease stage, anova test)
    #' model: (Expr ~ Disease_stage + Study,  data=df)
    #' fixed effect disease_stage and study
    #' @param full_report: report all the coefs for disease and gender
    #' @param to_plot, plot the residuals, histogram, box plots etc 
    #' model: fixed model

    
    ## plot name:
    plot_dir <- paste0(plot_dir, '/diagnostic_plots/')
    plot_pre <- paste0(plot_dir, i, '_')
    

    
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
    
    #print(str(df))
    full_model = lm(Expr ~ Disease_stage + Study, data=df)
    null_model = lm(Expr ~ Study,  data=df)
    
    
    ## get the pvalue by comparing to the null model
    (anova_r <- anova(full_model,null_model))
    (pvalue <- anova_r[, "Pr(>F)"][2])
    (chisq <- NA)
    (chi_df <- NA)
    
    ## extract information:
    ## get estimate, se, t of the disease_stage
    (df_t <- data.frame(geneSymbol =i, t(coef(summary(full_model))['Disease_stageDisease',])))
    
    
    
    if(full_report){ # report all coefficients
        ## get coefficients for each study (as fixed effects)

        (x <- as.data.frame(t(coefficients(full_model))))
        #print(i)
        ### get fixed effects of each study (label as intercept just for consistency with fixed models column names)
        ## one of the study is the baseline (=0)
        (base_study <- setdiff(levels(df$Study), gsub('Study', '', grep('Study', colnames(x), value = T))))
        baseline_study <- data.frame(x=0)
        colnames(baseline_study) <- paste0('Study', base_study)
        
        study_index <- grep('Study', colnames(x))
        (all_intercepts <- as.data.frame(x[,study_index ]))
        # cbind the baseline
        all_intercepts <- cbind(all_intercepts, baseline_study)
        colnames(all_intercepts) <- gsub('Study', '', colnames(all_intercepts))
        colnames(all_intercepts) <- paste0(colnames(all_intercepts), '_intercept')
        # get fixed effects of disease stage, FE of disease stage is the same for all studies
        df_effect <- all_intercepts
        df_effect[1, ] <- x[, "Disease_stageDisease"]
        colnames(df_effect) <- gsub('_intercept', '_FE_Disease',colnames(df_effect))
        
        ## confidence interval are not estimated
        ci <- data.frame(CI_2.5=NA, CI_97.5=NA)
        
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
#' fixed model : disease stage + study
#' @param exprdata : the expression.Rdata that contains array dat, msg, and array design
#' @param previous_result: mixed_model_results.Rdata, continue from the previous run, specify start rows

fixedModelAll <- function(rfile_ls=NULL, md_info=NULL, phase, out_dir, NA_filter=c('all', '0.3', 'none'),
                          width = 1200, height = 1200, 
                          rm_genes = "",
                          start_row =1, 
                          exprdata = NULL,
                          previous_result=NULL,
                          model=c('fixed_model'),
                          full_report = F,
                          to_plot= F,
                          tmp_result_output=T,
                          affy_threshold = 6, 
                          filter_method =c("median", "max", "min",  "mean"),
                          low_exp_rm=F){
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
    ## fixed model per gene
    #--------------#
    ## for loop here for each gene for the fixed model
    if(!('df_all' %in% ls())){df_all <- NULL} ## if a previous result is loaded
    print(paste0('startrow: ', start_row))
    genes_index <- start_row:nrow(array_dat)
    genes_index <- setdiff(genes_index, which(row.names(array_dat) %in% rm_genes))
    
    print(paste0('total genes : ', length(genes_index)))
    #####
    
    for (i in genes_index){
        df_t <- fixedModelStudy(array_dat, array_design, i, model=model, 
                                plot_dir = plot_dir,
                                full_report = full_report,
                                to_plot= to_plot,
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
        msg <- paste0(msg, '\nRemove ', length(rm_genes), ' genes in the fixed model due to ERROR, profiling detected new, lower deviance: \n   ',
                      paste0(rm_genes, collapse = ', '))
    }else{
        msg <- paste0(msg, '\n Remove 0 genes in the fixed model due to ERROR, profiling detected new, lower deviance')
    }
    
    
    # attach expression value
    df_all_expression <- mapBindAllColumns(df_all, array_dat)
    
    save(array_dat,array_design, df_all, msg, df_all_expression,file=paste0(result_dir, "mixed_model_results.Rdata"))
    writeTable(df=NULL, msg = msg, f_out = paste0(result_dir, 'mixed_model_results_log_', Sys.Date(), '.txt'))
    
    ##2016-09-20 temp
    model_formular = "Expr ~ Disease_stage + Study"
    
    returnlist=list(array_dat,array_design, df_all, msg, model_formular)
    print(paste0(result_dir, "mixed_model_results.Rdata"))
    return(returnlist)
}