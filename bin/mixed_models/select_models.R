#' 2016-09-20
#' this is for select models
#' input a random list of genes and get p values from the model test (from null to add more fixed effects) and 
#' see if adding an effect will affect the model

## plot the models for each gene (multiple models)
plotModelTest <- function(df, plot_name, test_model, i, plot_title){
    #i is the gene
    png(filename = plot_name, width =800, height =800)
    #graphics.off()
    par(mfrow=c(2,2))
    #1 QQ plot
    qqnorm(residuals(test_model), main =paste0(i, ': Normal Q-Q plot\n',plot_title))
    qqline(residuals(test_model), col = 2,lwd=2,lty=2)
    #2 residual vs fitted
    plot(fitted(test_model),residuals(test_model), main = 'Residuals vs. Fitted',xlab = 'Fitted', ylab = 'Residuals')
    abline(h=0,col = 2,lwd=2,lty=2)
    #3 check normality of residuals
    hist(residuals(test_model), main =paste0('Residuals'), xlab = 'Residuals')
    #4 sample expression
    boxplot(Expr ~ Disease_stage, data=df, ylab = 'Expression', main ='Gene expression')
    dev.off()
}




## get anova p values
## gender compared to null
#' age compared to null
#' if interaction: age + gender + disease_stage + age * disease  compared to age + gender + disease_stage + age * disease
modelTestPval <- function(df, fixed_effect, random_effect,significance_level,plot_dir,plt_prefix='',i){
    ## i is the gene name for plots
    (formular_rand <- paste0('(1|', random_effect, ')'))
    model0 <- lmer(Expr ~ (1|Study), data=df, REML=REML)
    effect_ls <- NULL
    
    test_model_effects <- NULL
    counter = 1
    if('Gender' %in% fixed_effect){  ## add gender as fixed
        (model_formular <-  paste0('Gender + ', formular_rand))
        assign(paste0('model_formular_', counter),model_formular)
        test_model_effects <- c(test_model_effects, model_formular)
        counter = counter +1
        effect_ls <- c(effect_ls, 'Gender')
    }
    if('Age' %in% fixed_effect){  ## add age as fixed
        model_formular <-  paste0('Age + ', formular_rand)
        assign(paste0('model_formular_', counter),model_formular)
        test_model_effects <- c(test_model_effects, model_formular)
        counter = counter +1
        effect_ls <- c(effect_ls, 'Age')
    }
    
    ############
    ## run lmer different effect against NULL
    ############
    (formular_ls <- grep('model_formular_',ls(), value = T))
    
    for(j in 1:length(formular_ls)){
        (formular <- formular_ls[j])
        (cmd <- paste0('lmer(Expr ~ ', eval(as.symbol(formular)), ', data=df, REML=REML)'))
        test_model <- eval(parse(text = cmd))
        assign(paste0('model_lmer_', j), test_model)
        
        ## plot the model test
        (plot_name <- paste0(plot_dir, '/plots/'))
        dir.create(plot_name, recursive = T, showWarnings = F)
        (plot_name <- paste0(plot_name, i, '_model_test_plots_model_',j, plt_prefix, '.png'))
        (plot_title <- eval(as.symbol(formular)))
        plotModelTest(df, plot_name, test_model, i, plot_title)

    }
    (modellmer_ls <- grep('model_lmer_',ls(), value = T))
    
    
    ### test for each effect against NULL
    df_test_ls <- NULL
    for(j in 1:length(modellmer_ls)){
        (cmd <- paste0('anova_r <- anova(model0, ', modellmer_ls[j], ')'))
        eval(parse(text = cmd))
        ## get the pvalues
        (pvalue <- anova_r[, "Pr(>Chisq)"][2])
        (chisq <- anova_r[, "Chisq"][2])
        (dfree <- anova_r[, "Df"][2])
        df_test <- data.frame(geneSymbol =i, anova_p = pvalue, chisq = chisq, degreefreedom = dfree, 
                              model = test_model_effects[j], effect = effect_ls[j])

        if(is.null(df_test_ls)){
            df_test_ls <- df_test
        }else{
            df_test_ls <- rbind(df_test_ls, df_test)
        }
    }

    


    
    #############
    ####### check if disease age interaction is important
    #############
        ## add disease age interaction
        if('Age' %in% fixed_effect){  ## add age as fixed
            if('Gender' %in% fixed_effect){
                (model_formular <-  'Age * Disease_stage + Gender +(1|Study)')
                null_formular <- 'Age + Disease_stage + Gender +(1|Study)'
            }else{
                (model_formular <-  'Age * Disease_stage +(1|Study)')
                null_formular <- 'Age + Disease_stage +(1|Study)'
            }
            (cmd <- paste0('lmer(Expr ~ ', eval(as.symbol(formular)), ', data=df, REML=REML)'))
            test_model <- eval(parse(text = cmd))
            
            
            assign('model_interaction',eval(parse(text = paste0('lmer(Expr ~ ', model_formular, ', data=df, REML=REML)'))))
            assign('model_null',eval(parse(text = paste0('lmer(Expr ~ ', null_formular, ', data=df, REML=REML)'))))
            anova_r <- anova(model_null, model_interaction)

            ## get the pvalues
            (pvalue <- anova_r[, "Pr(>Chisq)"][2])
            (chisq <- anova_r[, "Chisq"][2])
            (dfree <- anova_r[, "Df"][2])
            
            df_test <- data.frame(geneSymbol =i, anova_p = pvalue, chisq = chisq, degreefreedom = dfree, 
                                  model = model_formular, effect = 'disease and age interaction')
            df_test_ls <- rbind(df_test_ls, df_test)
        }
    
    
    df_test_ls$significant <- NA
    df_test_ls$significant[which(df_test_ls$anova_p<=significance_level)] <- 'yes'
    
    
    
    
    (index <- intersect(which(df_test_ls$effect %in% setdiff(fixed_effect, 'Disease_stage')), which(is.na(df_test_ls$significant))))
    if(length(index) >0){
        (not_sig_effect <- as.character(df_test_ls$effect[index]))
        (fixed_effect <- setdiff(fixed_effect, not_sig_effect))
        with_non_sig <- 'yes'
    }else{
        with_non_sig <- 'no'
    }
    returnlist <- list(df_test_ls,fixed_effect, with_non_sig)
    return(returnlist)
}

#------------------------------------------------------------------#
# for each gene, do the model selection
#------------------------------------------------------------------#
mixedModelSelection <- function(array_dat, array_design, i, 
                            model=c('random_intercept'),
                            fixed_effect = c('Gender', 'Age'), 
                            random_effect = c('Study'),
                            plot_dir,
                            REML=FALSE,
                            significance_level=0.05){
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

    
    ## preprocess to get the expression df for lmer()
    df <- as.data.frame(t(array_dat[i, ])) # expression of the genes from all samples
    
    i <- colnames(df)
    print(i)
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

    #########
    ## check for confound fixed effects(to study)
    #########
    ## if each study only contains 1 gender 
    
    if(any(length(levels(df$Gender))==1, length(levels(df$Gender))==F)){ ## if only 1 gender or gender is NA
        fixed_effect <-setdiff(fixed_effect,'Gender')
    }
#     (tmp <- rmDup(df[,c('Gender','Study')]))
#     if(!any(as.numeric(table(tmp[,'Study']))>1)){## remove gender if T (F when at least 1 study contains both gender), also rm when gender =NA
#         fixed_effect <-setdiff(fixed_effect,'Gender')
#     }   
    
    ## if age group only belongs 1 study (i.e. each study has distinct 1 age group)
    (tmp <- rmDup(df[,c('Age','Study')]))
    if(!any(as.numeric(table(tmp[,'Age']))>1)){## rm age if T (F when at least 1 age belongs to multiple studies)
        fixed_effect <-setdiff(fixed_effect,'Age')
    }  
    ############
    ##formulate the different models
    ############
    x <- modelTestPval(df, fixed_effect, random_effect,significance_level, plot_dir = plot_dir, plt_prefix = '_first',i)
    df_test <- x[[1]]
    fixed_effect <- x[[2]]
    with_non_sig <- x[[3]]
    df_test$test <- 'first'

    
    returnlist = list(df, i, df_test)
    return(returnlist)
}



#------------------------------------------------------------------#
# loop for a gene list
#------------------------------------------------------------------#
mixedModelSelectionLoop <- function(array_dat, array_design, genelist, 
                                model=c('random_intercept'),
                                fixed_effect = c('Gender', 'Age','Disease_stage'), 
                                random_effect = c('Study'),
                                plot_dir,
                                REML=FALSE,
                                significance_level=0.05){
    #' genelist: a list of genes names or index that corresponding to the array_dat
    
    ## output dir name:
    plot_dir <-  paste0(plot_dir, '/', 'model_selection/', Sys.Date(),'/')
    dir.create(plot_dir,showWarnings = F, recursive = T)
    
    df_test_all <- NULL
    for(i in genelist){
        print(i)
        x <- mixedModelSelection(array_dat, array_design, i, 
                                        model,
                                        fixed_effect, 
                                        random_effect,
                                        plot_dir,
                                        REML,
                                        significance_level)
        df_test <- x[[3]]  ## the pvalue test for each gene
        if(is.null(df_test_all)){
            df_test_all <- df_test
        }else{
            df_test_all <- rbind(df_test_all, df_test)
        }
    }
    ## write test pvalues
    out_dir <- paste0(plot_dir, '/model_selection_results/')
    dir.create(out_dir,showWarnings = F)
    writeTable(df_test_all, f_out = paste0(out_dir,'model_test_pvalues.tsv'))
    
    
    ## summary
    ## genes with any model with a significant
    (sig_genes<- unique(as.character(df_test_all[which(df_test_all$significant=='yes'), 'geneSymbol'])))
    df <- df_test_all[which(df_test_all$geneSymbol %in% sig_genes), ]
    df_t <- na.omit(df)%>%droplevels()
    df_summary <- data.frame(table(df_t$effect))
    colnames(df_summary) <- c('effect', 'freq')
    df_summary$total_sig_genes <- NA
    for(k in 1:nrow(df_summary)){
        eff <- as.character(df_summary$effect[k])
        df_summary$total_sig_genes[k] <- length(which(df$effect %in% eff))
    }
    
    df_summary$ratio<- df_summary$freq/df_summary$total_sig_genes
    writeTable(df_summary, f_out = paste0(out_dir,'model_test_pvalues_summary.tsv'))
    writeTable(data.frame(sig_genes = sig_genes), f_out = paste0(out_dir,'model_test_sig_genes.tsv'))
    return(df_test_all)
}



# ## get anova p values
# modelTestPval <- function(df, fixed_effect, random_effect,significance_level,plot_dir,plt_prefix='',i){
#     ## i is the gene name for plots
#     (formular_rand <- paste0('(1|', random_effect, ')'))
#     model_formular_0 <- formular_rand  ## null only with study(random)
#     model_formular <- formular_rand  ## null only with study(random)
#     test_model_effects <- model_formular
#     effect_ls <- c(random_effect, fixed_effect)
#     
#     counter = 1
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
#     ############
#     ## run lmer on the different models
#     ############
#     formular_ls <- grep('model_formular_',ls(), value = T)
#     
#     for(j in 1:length(formular_ls)){
#         formular <- formular_ls[j]
#         (cmd <- paste0('lmer(Expr ~ ', eval(as.symbol(formular)), ', data=df, REML=REML)'))
#         test_model <- eval(parse(text = cmd))
#         assign(paste0('model_lmer_', j), test_model)
#         
#         ## plot the model test
#         (plot_name <- paste0(plot_dir, '/', i, '_model_test_plots_model_',j, plt_prefix, '.png'))
#         (plot_title <- eval(as.symbol(formular)))
#         plotModelTest(df, plot_name, test_model, i, plot_title)
# 
#     }
#     (modellmer_ls <- grep('model_lmer_',ls(), value = T))
#     
#     (cmd <- paste0('anova_r <- anova(', paste0(modellmer_ls, collapse = ', '), ')'))
#     eval(parse(text = cmd))
#     
#     ## get the pvalues
#     (pvalue <- anova_r[, "Pr(>Chisq)"])
#     (chisq <- anova_r[, "Chisq"])
#     (dfree <- anova_r[, "Df"])
#     
#     #############
#     ## pvalue df
#     #############
#     df_test <- data.frame(geneSymbol =i, anova_p = pvalue, chisq = chisq, degreefreedom = dfree, model = test_model_effects, effect = effect_ls)
#     df_test$significant <- NA
#     df_test$significant[which(df_test$anova_p<=significance_level)] <- 'yes'
#     
#     #############
#     ####### check if a fixed effect if not significant (except for disease stage)
#     #############
#     (index <- intersect(which(df_test$effect %in% setdiff(fixed_effect, 'Disease_stage')), which(is.na(df_test$significant))))
#     if(length(index) >0){
#         (not_sig_effect <- as.character(df_test$effect[index]))
#         (fixed_effect <- setdiff(fixed_effect, not_sig_effect))
#         with_non_sig <- 'yes'
#     }else{
#         with_non_sig <- 'no'
#     }
#     returnlist <- list(df_test,fixed_effect, with_non_sig)
#     return(returnlist)
# }
