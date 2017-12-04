#' 2016-03-14
#' last update 2016-04-25
#' use wrapper
#' the script only finds non-expressed gender genes thesholds 
#' by selected method (one of method_t_ls <-c("max", "median", "75quantile", "95quantile") )
#' expressed data are not actually filtered. use 'filter_low_exp_and_save.R' to actually filter 
#' genes by gender gene expression
#' 
#' 
#' ################################################################################################

#' ----------------------------------------------------- #
#' USAGE:

# datadir <-'/home/bzhuang/AD_mouse_model_project/data_and_QC/all_data/explore_data//analysis_mode/'  ## dir for limma toptables
# df_info <- '/home/bzhuang/AD_mouse_model_project/config_files/dataset_info.tsv'  #dataset label, phase label, order of the datasets etc.
# plot_dir <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/all_data/explore_data/gender_expression_threshold/"
# 
# (r_obj_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
# method_t_ls <-c("max", "median", "75quantile", "95quantile") 
# 
# df <- loopExpThreshold(method_t_ls, r_obj_ls, plot_dir, width=1000, height=800, df_return =T)
#' ----------------------------------------------------- # 
#' Notes: 


library(dplyr)
library(tidyr)
#############################
## functions
#############################

probeExp <- function(r_obj, df_gene){
    #' r_obj dir for object (array_dat, annotation etc)
    #' @return df_exp_t: exp of non expressed gender genes
    #'          df_exp: exp of all gender genes
    #'          dataset,platform, array_dat, 
    ###-------------------------
    ## probe expression
    ###-------------------------
    
    ## get the probes
    load(r_obj)
    df_gene <- na.omit(noWarnings(left_join(df_gene, annotation)))
    (probe_ls <- as.character(df_gene$ProbeName))
    
    ## make sure the probe is in the array dat too
    (probe_ls <- intersect(probe_ls, row.names(array_dat)))
    
    ## get the probes and expression
    (df_exp <- array_dat[probe_ls, ] %>% droplevels)
    df_exp$ProbeName <- row.names(df_exp)
    df_exp <- gather(df_exp, ProbeName, Expression)
    colnames(df_exp) <- c("ProbeName", "Sample", "Expression")
    
    ## combine with the array_design, and gene names
    df_exp <- noWarnings(left_join(df_exp, array_design[, c("Sample", "Gender")]))
    df_exp <- noWarnings(left_join(df_exp, df_gene[, c("GenderGenes", "GeneSymbols", "ProbeName")]))
    df_exp$Sample <- as.factor(df_exp$Sample)
    
    ## sort df_exp by gender
    df_exp <- df_exp[order(df_exp$Gender), ]

    ###-------------------------
    ## get the gender of samples and threshold cut off for non expressed geneder gene
    ###-------------------------
    ## get the expression other gendergenes (i.e. probes that represent no expression)
    ## if the sample is mixed gender (multiple samples combined, then ignore these samples)
    index <-NULL
    for (i in 1: nrow(df_exp)){
        if (as.character(df_exp$Gender)[i] != 'mixed'){
            if(as.character(df_exp$Gender)[i] != as.character(df_exp$GenderGenes)[i]){
                index <- c(index, i)
            }
        }

    }
    df_exp_t <- df_exp[index, ]%>%droplevels()
    
    returnlist <- list(df_exp_t, dataset,platform, array_dat, df_exp, 
                       annotation, array_design) # return all to save for a new object for limma
    return(returnlist)
}





#### get the probe expression first and find threshold for each method and see the removal of probes
findThreshold <- function(method_t, df_exp_t, dataset,platform, array_dat, array_design, df_exp,
                          plot_dir="", width=1000, height=800, to_plot =T, msg ="#"){
    
    #df_exp : all gender probes expression
    #df_exp_t: gender probes of the opposite sex (i.e. none expressed probes)

    
    if(length(method_t) >1){method_t <- method_t[1]}
    if(method_t == "median"){
        (threshold <- median(df_exp_t$Expression))
    }else if (method_t == "max"){
        (threshold <- max(df_exp_t$Expression))
    }else if (method_t == "75quantile"){
        (threshold <- quantile(df_exp_t$Expression, probs = 0.75))
    }else if (method_t == "95quantile"){
        (threshold <- quantile(df_exp_t$Expression, probs = 0.95))
    }
    
    if('Timepoint' %in% colnames(array_design)){
        timepoint <- levels(array_design$Timepoint)
        if(length(timepoint) >1){ timepoint <- paste0(timepoint, collapse = '; ')}
    }else{
        timepoint <- ''
    }
    
    if(timepoint !=''){
        (plot_pre <- paste0(plot_dir, dataset, "_",timepoint, "_", method_t, '_'))
        dataset_p<- paste0(dataset, '(', timepoint,')')  ## add the timepoint to the plot title (not in the saved plotname)
    }else{
        dataset_p<- dataset
        (plot_pre <- paste0(plot_dir, dataset, "_", method_t, '_'))
    }
    
    if(to_plot){
        ## label the x labels by gender
        df_tmp <-df_exp[, c("Sample", "Gender")]
        df_tmp <- df_tmp[!duplicated(df_tmp), ] %>%droplevels()
        df_tmp <- noWarnings(left_join(data.frame(Sample = unique(df_exp$Sample)), df_tmp))
        df_tmp$Gender <- gsub("F", "female", df_tmp$Gender)
        df_tmp$Gender <- gsub("M", "male", df_tmp$Gender)
        (label_col <- mapValue(df_tmp$Gender, c(male="blue", female="red", mixed ="brown")))
        ## must reorder the factor levels same as the new sample order
        # (otherwise the label of gender colour will not match)
        df_exp$Sample <- factor(df_exp$Sample, levels = unique(df_exp$Sample))  
        
        
        ## save plots
        p<- plotXY(df_exp, y="Expression", x="Sample", group ="ProbeName",plt_alpha = 1, 
                   title = paste0(dataset_p, ": Expression of gender genes.\n Threshold: ", round(threshold, 3), 
                                  ", Method: ", method_t))
        p + geom_hline(yintercept = threshold, show.legend = T)+ theme(axis.text.x = element_text(colour=label_col))
        ggsave(filename= paste0(plot_pre,'scatterplot.png'),width = 6, height = 5)
        
        p<- plotXY(df_exp, y="Expression", x="Sample", group ="GenderGenes",plt_alpha = 1, 
                   title = paste0(dataset_p, ": Expression of gender genes.\n Threshold: ", round(threshold, 3), 
                                  ", Method: ", method_t))
        p + geom_hline(yintercept = threshold, show.legend = T)+ theme(axis.text.x = element_text(colour=label_col))
        ggsave(filename= paste0(plot_pre,'scatterplot_by_gender_genes.png'),width = 6, height = 5)
        
        png(filename= paste0(plot_pre,'dens_all.png'), width = width, height = height)
        plot(density(df_exp$Expression), main = paste0(dataset, ": all gender probes "))
        dev.off()
#         png(filename= paste0(plot_pre,'dens_non_exp.png'), width = width, height = height)
#         plot(density(df_exp_t$Expression), main = paste0(dataset, ": non-expressed gender probes "))
#         dev.off()
    }

    
    #**************************
    ## remove probes with median < threshold
    #**************************

    x <- Biobase::rowMedians(as.matrix(array_dat))
    index <- which(x > threshold)
    
    ## new array with non expressed probes removed
    array_new <- array_dat[index, ]%>% droplevels()
    
    (removed_ratio <- round((nrow(array_dat)-nrow(array_new))/nrow(array_dat), 4))
    
    ## df is recording the threshold
    df <- data.frame(dataset = dataset, Timepoint = timepoint, platform = platform, 
                     threshold = threshold, method = method_t,
                     removed_ratio = removed_ratio, prefilter_probes = nrow(array_dat), 
                     after_filter_probes =nrow(array_new))
    
    ## log
    msg_rm <- paste0("\n#-------------------------------#\n## FILTER LOW EXP PROBES\n#-------------------------------#",
                     "\n# ", Sys.Date(),
                     "\n# Non-expressed gender genes for filter: ", paste0(unique(df_exp_t$GeneSymbols), collapse=", "),
                     "\n# Expression threshold: ", threshold, "(", method_t, ")",
                     "\n# Before filter: ", nrow(array_dat), "; after filter: ", nrow(array_new), ". ", removed_ratio*100, "% under threshold.")
                     
    msg <- paste0(msg, msg_rm)
    returnlist <- list(df, array_new, msg)
    
    print(paste0(dataset, ": method: ", method_t, ": threshold: ", threshold, "; ",removed_ratio*100, "% below threshold."))
    return(returnlist)
}



#############################
## Main function
#############################
## check threshold for different methods for a dataset
expThreshold <-function(method_t_ls, r_obj, df_gene,
                            plot_dir, width=1000, height=800){
    
    x <- probeExp(r_obj, df_gene)
    df_exp_t <- x[[1]]
    dataset <- x[[2]]
    platform <- x[[3]]
    array_dat <- x[[4]]
    df_exp <- x[[5]]
    array_design <-x[[7]]
    
    df_all <- NULL
    for (method_t in method_t_ls){
        y <- findThreshold(method_t, df_exp_t, dataset,platform, array_dat, array_design, df_exp,
                      plot_dir, width, height, to_plot =T)
        df <- y[[1]]
        df_all <- rbind(df_all, df)
    }
    return(df_all)
}



## given a list of .Rdata objects (only from explore data Rdata (not after limma)) and calculate the thresholds
## and plot
loopExpThreshold <- function(method_t_ls, r_obj_ls,
                         plot_dir, width=1000, height=800, df_return =F){
    
    dir.create(plot_dir, showWarnings=F,recursive=T)
    
    ## define gender genes
    df_gene <- read.table(header=TRUE, text='
                         GenderGenes    GeneSymbols
                         F    Xist
                         M    Kdm5d
                         M    Rps4y1')
    
    ### loop for all datasets and all methods
    df_all <- NULL
    for (r_obj in r_obj_ls){
        df <- expThreshold (method_t_ls, r_obj, df_gene,
                            plot_dir, width=width, height=height)
        df_all <- rbind(df_all, df)
        write.table(df_all, file = paste0(plot_dir, "stat_thresholds.tsv"), quote = F, sep ='\t',
                    row.names =F)
    }
    print(paste0('figures are in ', plot_dir))
    if(df_return){return(df_all)}

}





# 
# 
# ## to use a method the filter out the probes and save as new object(specify 1 method)
# # # df_gene: the df for gender gene: e.g.
# # df_gene<- read.table(header=TRUE, text='
# #                      GenderGenes    GeneSymbols
# #                      F    Xist
# #                      M    Kdm5d
# #                      M    Rps4y1')
# mainExpThreshold <-function(method_t, r_obj, df_gene,
#                             plot_dir, width=1000, height=800){
#     
#     x <- probeExp(r_obj, df_gene)
#     df_exp_t <- x[[1]]
#     dataset <- x[[2]]
#     platform <- x[[3]]
#     array_dat <- x[[4]]
#     df_exp <- x[[5]]
#     
#     df_all <- NULL
#     for (method_t in method_t_ls){
#         y <- findThreshold(method_t, df_exp_t, dataset,platform, array_dat, df_exp,
#                             plot_dir, width, height, to_plot =T)
#         df <- y[[1]]
#         df_all <- rbind(df_all, df)
#     }
#     return(df_all)
# }





