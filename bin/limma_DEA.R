#' 2015-11-24
#' update 2016-05-18
#' ################################################################################################
#' Limma for exp data
#' 
#' INPUT: a file contain the model info'/home/bzhuang/AD_mouse_model_project/config_files/mouse_datasets_limma.tsv'
#'      - the table file must contain the following columns: 
#'      'eeName': GSE_id, e.g. GSE1234
#'      'platform_id': e.g. GPL1261
#'      'subset': 'yes' indicate that samples needed to be subsetted, 'no' otherwise
#'      'subset_by': if subset is 'yes', fill in this column of which factor needs to be subsetted 
#'                    and separated by '; '(with a space) if multiple factors are required to be subsetted. e.g. 'OrganismPart; Timepoint'
#'      'keep_subset': indicate the levels of a factor or multiple factors to keep for analysis. order is the same as order in 'subset_by'
#'                    column. e.g. "c('Hippocampus'); c('6_months', '15_months')"
#'      'batch_corrected': 'yes' to retrieve batch corrected expression matrix, 'no' otherwise
#'      'factors': factors to be used in the limma model, e.g. c('Genotype', 'Batch')
#'      'factor_levels': specify factor levels which update the baseline level in the limma model, sep by "; ", 
#'              e.g 'array_design$Genotype <- factor(array_design$Genotype, levels = c('WT', 'J20')); array_design$Treatment <- factor(array_design$Treatment, levels = c('reference', 'GFP'))'
#' 
#' OUTPUT: plots, up and down regulation 1 sided p values (default method 'BH') for meta analysis, and save processed toptable for each coef
#' NOTES:
#'    - 
#' ----------------------------------------------------- #
#' USAGE: to run by range of rows, specified rows, or datasets
# rm(list=setdiff(ls(),'home_dir'))
# q_threshold <- 0.05
# adj_method <- "BH"
# datadir <- '/home/bzhuang/gemma_data/explore_data/prioritized_limma/'
# object_dir <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/gemma_data/explore_data/'
# f_model <- '/home/bzhuang/AD_mouse_model_project/config_files/prioritized_mouse_datasets_limma.tsv'
# start_row <- 1 # or specify which rows: to_do <- c(9, 10) # or specify which dataset: dataset_todo <-c('GSE36237', 'GSE1556')
# source('limma_DEA.R')
#' ----------------------------------------------------- # 
#' 
#' ################################################################################################


#'*************************************************************************#
#' PART 1: SET SYS DEFAULTS
#'   - Load packages, source and define defaults for img size, and set wd
#' PART 2: CHECK INPUTS and THRESHOLDS
#'   - q_threshold is default 0.05
#'   - check if datadir(for output), object_dir(input), and f_model(input model info) are defined
#' PART 3: DEFINE HELPER FUNCTIONS
#'   - heatmapForTop200(): save heatmap for top 200 genes
#' PART 4: MAIN FUNCTION
#'   - mainLimmaDEA(): main DEA for one dataset and one model
#'   - output: up and down regulation 1 sided p values (default method 'BH') for meta analysis, and save processed toptable for each coef
#' Part 5: LOOP
#'   - loop for each predefined limma model for each dataset
#**************************************************************************#

#---------------------------------------------------------------------------#
# PART 1: SET SYS DEFAULTS
#---------------------------------------------------------------------------#
# load packages
package_ls <- c('plyr', 'dplyr','ggplot2', 'limma', 'RColorBrewer')
for(i in package_ls){
  cmd <- paste0("if('package:", i, "' %in% search() == FALSE){library(", i, ")}")
  eval(parse(text =cmd))
}

 
source('helper_functions.R')

## set defaults
## image size (pch)
width = 1200
height = 800



#---------------------------------------------------------------------------#
# PART 2: SET INPUTS and THRESHOLDS
#---------------------------------------------------------------------------#
## default p value adjustment for multiple testing
if(!exists("adj_method")){adj_method <- "BH"}


## FRD value threshold (default = 0.05)
if(!exists("q_threshold")){q_threshold <- 0.05}

## check if inputs are defined
if(!exists("datadir")){
  stop(paste0("need to define datadir for the result output"))}
if(!exists("object_dir")){
  stop(paste0("need to define object_dir for to retrieve R objects"))}
if(!exists("f_model")){
  stop(paste0("need to define f_model(dir) for the limma model input"))}

dir.create(datadir, showWarnings=F, recursive=T)


#---------------------------------------------------------------------------#
# PART 3.1: DEFINE HELPER FUNCTIONS-heatmap
#---------------------------------------------------------------------------#
heatmapForTop200 <- function(df, width, height,gene_annotation, array_design, target_coef, 
                             baseline, plot_pre, q_threshold,prefix = "all", timepoint =""){
  ## save heatmap for top 200 DE genes
  ## INPUT:
  ##   - df: the matrix with only the DE genes of all samples, samples are ordered
  ##   - prefix(string): specify input is all DE, only up or down regulated. prefix = one of "all", "up", "down"
  ##   - width and height: are the default img dimensions
  ##   - gene_annotation(df): with ProbeName and GeneSymbols
  ##   - target_coef, baseline are strings for titles
  ##   - q_threshold: for plot title and plot filename
  ##   - timepoint: a string define the mouse age (for plot name and figure title)    
  ## Notes:
  ##   - no plot output if df has 0 rows under threshold
  ##   - heatmap are clustered by rows
  ##   - scale_row = T: z_score scale by genes
  
  #' scale by row or column: 
  #' If we want to show differences between genes: Z-score by samples(columns) (force each sample to have zero mean and standard deviation=1). 
  #' If we want to show differences between samples (show DE genes), Z-score by genes(rows) (force each gene to have zero mean and standard deviation=1).
  
  if(nrow(df) ==0){
    print('No DE genes in the input for heatmap')
  }else{
    if(tolower(prefix) %in% c("all", "de")){
      prefix <- 'DE'
    }else if(tolower(prefix) == 'up'){
      prefix <- 'up_DE'
    }else{
      prefix <- 'down_DE'
    }
    if(timepoint ==""){
        dataset_p<- dataset
    }else{
        dataset_p<- paste0(dataset, '(', timepoint,')')  ## add the timepoint to the plot title (not in the saved plotname)
    }
    
    ## only plot the top 200 probes
    ## get the gene names
    gene_names_df <- data.frame(ProbeName = rownames(df))
    gene_names_df <- noWarnings(left_join(gene_names_df, gene_annotation))
    rownames(df) <- paste(gene_names_df$GeneSymbols, gene_names_df$ProbeName, sep = ':')
    ## truncate to the top 200 for heatmap
    if(nrow(df)>200){
      df_heatmap <- df[1:200, ] %>% droplevels()
      heatmap_title <- paste0(dataset_p, ' top 200 ', prefix, ' genes, q threshold: ', q_threshold, 
                              '\n', target_coef, '-vs-baseline_', baseline)
      f_heatmap <- paste0(plot_pre, prefix, '_top200_', target_coef, '-vs-baseline_', baseline, 
                          '_q_threshold', q_threshold, '.png')
    }else{
      df_heatmap <- df
      heatmap_title <- paste0(dataset_p, ' ', prefix, ' genes: ', nrow(df), ', q threshold: ', 
                              q_threshold, '\n', target_coef, '-vs-baseline_', baseline)
      f_heatmap <- paste0(plot_pre, prefix, '_', nrow(df), '_', target_coef, '-vs-baseline_', baseline, 
                          '_q_threshold', q_threshold, '.png')
    }
    
    ## set img dimensions
    heat_height <- nrow(df_heatmap)*15
    heat_width <- ncol(df_heatmap)*35
    if(heat_height < height){heat_height <- height}
    if(heat_width < width){heat_width <- width}
    
    ## plot heatmap for DE genes
    if(nrow(df_heatmap)==1){
    png(filename= f_heatmap, width = heat_width, height = heat_height)
    plotHeatmap(df_heatmap, title = heatmap_title, size_r=12, size_c =15, cluster_rows = F,array_design_df = array_design)
    dev.off()
    } else {
      png(filename= f_heatmap, width = heat_width, height = heat_height)
      plotHeatmap(df_heatmap, title = heatmap_title, size_r=12, size_c =15, cluster_rows = T,
                  array_design_df = array_design, scale_row = T)
      dev.off()
    }

    print(paste0("heatmap for ", prefix, " saved, ", f_heatmap))
  } 
}


#---------------------------------------------------------------------------#
# PART 3.2: DEFINE HELPER FUNCTIONS- process probesets into gene pvalues
#---------------------------------------------------------------------------#
#' Process a data.frame with 1-sided p values of probesets and return a df with 1-sided p values of genes for metaanalysis
#' @param df A dataframe with 33 columns: "ProbeName"   "GeneSymbols" "up"/"down" regulated p values
#' @param dataset (str) dataset id (for colname)
#' @param timepoint (prefixed with '_') for colname, e.g. _2_months
#' @param target_coef (str) coef lable e.g. 'GenotypeWT'
#' @param adj_method (str) A p value correction method for genes with mulitple probes mapped (e.g. 'BH')
#' @return
#' A list contain a df and msg
#' df: a processed dataframe with GeneSymbols and 1-sided pvalues (2 columns). If a gene has multiple probesets mapped, p values is adjusted and 
#'          the min adj p is selected for the gene. The returned df is for meta-analysis
#' msg: message
#' @note
#'  - probes not mapped to a gene are removed
#'  - probes mapped to multiple genes are removed
#'  - if a gene has multiple probes mapped, p values is adjusted and the min adj p is selected for the gene
#' @examples
#' x <- aggregateGenePvals(tmp_df_down, dataset, timepoint, target_coef, adj_method, plot = F, verbose =T)
#' df1<- x[[1]]
#' msg_meta <- x[[2]]
aggregateGenePvals <- function(df, dataset, timepoint, target_coef, adj_method, plot = F, verbose =T){
       ###############################
       #  get defaults and sanity check
       ###############################
       ## get the dataset id
       f_label <- dataset
       ## extract geneotype names for column names, otherwise use target coef
       if(grepl("Genotype", target_coef)){
           (f_label_1 <- unlist(strsplit(target_coef, split ="Genotype"))[2])
       }else {f_label_1 <- target_coef}
       ##add the timepoint
       (f_label <-paste0(f_label, timepoint, '_', f_label_1))
       
       ## check the format of the input table
       if(!all.equal(colnames(df)[1:2], c("ProbeName", "GeneSymbols"))){
           stop(paste0("Input file: ", f_to_process, " has columns:  ", 
                       paste(colnames(df), collapse=","),
                       ". Need to have 'ProbeName', 'GeneSymbols' as the first two columns."))}
       
       if(ncol(df) != 3){
           stop(paste0("Input file: ", f_to_process, " has columns:  ", 
                       paste(colnames(df), collapse=","),
                       ". Need to have 3 columns: 'ProbeName', 'GeneSymbols', 
                       and pvalue (up or down regualtion)."))}
       (msg <- paste0(f_label, ': original input has ', nrow(df), " probesets; "))
       ## update column names
       colnames(df)[3] <- 'pvalue'
       
       
       ###############################
       #  remove probes that don't map to any gene symbols
       ###############################
       df <- excludeMatch(df, "GeneSymbols", "")
       (msg <- paste0(msg, 'mapped probesets: ', nrow(df), "; "))
       
       
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
       #  #  adjust p values
       ###############################
       #' > In the case where a gene had more than one probeset assigned to it, the p-values for the probesets are 
       #' Bonferroni-corrected(default) and the lowest corrected p-value was used to represent the gene for that dataset
       pvalAdj <- function(x){
           min(p.adjust(x,adj_method))
       }
       
       #   how to deal with mulitple probes mapped to a gene:
       #     
       #     > In the case where a gene had more than one probeset assigned to it, the p-values for the probesets aere Bonferroni-corrected and the lowest corrected p-value was used to represent the gene for that dataset
       #   
       df <-  aggregate(pvalue ~ GeneSymbols, df, function(x) pvalAdj(x))
       ## sort by pval
       df <- df[order(df$pvalue),]
       
       
       #### sanity check
       if(plot){
           plot(density(df$pvalue), main = f_label)  
       }
       
       ## rename column
       df2 <- as.data.frame(df[, 'pvalue'])
       colnames(df2) <- f_label
       df2$GeneSymbols <- df$GeneSymbols

       (msg <- paste0(msg, 'genes with multiple probesets are ', adj_method, ' corrected. New data has ', nrow(df2), ' genes.'))
       if(verbose){print(msg)}
       
       returnlist <- list(df2, msg)
       return(returnlist)
}


#---------------------------------------------------------------------------#
# PART 4: DEFINE MAIN FUNCTIONS
#---------------------------------------------------------------------------#
mainLimmaDEA <- function(array_dat, array_design, annotation, model_info, q_threshold,
                         preprocess_msg,adj_method,plot_heatmap=T,
                         combined_genotype =F, write_table=T){
  ## for a given model info, process expression data and DE by limma (for one dataset and one model)
  ## INPUT:
  ##   - array_dat(df), array_design(df), annotation(df): loaded from preprocessed R objects 
  ##   - model_info(df): a df with one row of model info
  ##   - q_threshold(numeric), a threshold for FDR value for Limma results
  ##   - preprocess_msg(str): msg from the preprocessing, add as comments in the saved tables
  ##   - adj_method(str): method for p value correlation for multiple testing in function p.adjust(), e.g. 'BH', "bonferroni", 'fdr',
  ##      used for genes for multiple probesets
  ## OUTPUT:
  ##   - 
  ##   - 
  ##   - 
  ##   - 
  ## NOTES:
  ##   - the pvalue from limma is processed to up or down regualtion (1 sided pvalue)- will be used for meta analysis (don't use adj p)
  ##   - save toptable contains adj 1-sided pvalues (e.g. 'BH' for pvalue correction)
  
  if(nrow(model_info) != 1){
    stop("Error: model_info should only contain 1 row of record")
  }
  
  dataset <- as.character(model_info$eeName)
  print(paste0("DEA process of dataset ", dataset))
  dir.create(paste0(datadir, dataset), showWarnings=F, recursive=T)


  
  
  
  ###############################
  ###### subset data if required
  ###############################
  ### 1. get the list for subset factors
  ### 2. get the new data matrix and design
  ### 3. and assign folder name (subset and with or without batch correction)
  ###############################
  subset_folder <- paste0(datadir, dataset, '/')
  msg_batch <- ""
  if(tolower(as.character(model_info$batch_corrected)) == "yes"){
    subset_folder <- paste0(subset_folder, 'batch_corrected_')
    msg_batch <- "\nInput experssion matrix is batch corrected.\n"
  }
  if(tolower(model_info$subset) == 'yes'){
      subset_factor_ls <- unlist(strsplit(as.character(model_info$subset_by), split ="; "))
      subset_keep_factor_ls <- unlist(strsplit(as.character(model_info$keep_subset), split ="; "))
      x <- subsetSamples(model_info, dataset, array_dat, array_design, to_log =T, df_subset =F)
      array_dat <- x[[1]]
      array_design <- x[[2]]
      msg_subset <- paste0(msg_batch, as.character(x[[3]]))
      ## reformat subfolder name
      tmp <- paste(paste(subset_factor_ls,subset_keep_factor_ls, sep = ""), collapse="_")
      tmp <- gsub("\'", "", tmp)
      tmp <- gsub("c\\(", "\\_", tmp)
      tmp <- gsub("\\)", "\\", tmp)
      tmp <- gsub("\\ ", "", tmp)
      subset_folder <- paste0(subset_folder, 'subset/',tmp, '/')
  } else {
      subset_folder <- paste0(subset_folder, 'no_subset/')
      msg_subset <- msg_batch
  }

  

  ###############################
  ######### make dirs for output
  ###############################
  plot_dir <-paste0(subset_folder, 'plots/')
  result_dir <-paste0(subset_folder,'results/')
  dir.create(plot_dir, showWarnings= F, recursive =T)
  dir.create(result_dir, showWarnings= F, recursive =T)
  ## redefine plot and result dir (will add the timepoint after)
  plot_pre <- paste0(plot_dir, dataset, '_')
  result_pre <- paste0(result_dir, dataset, '_')
  
  ## get the timepoint from design (for dataset with only 1 timepoint, use for file names)
  if('Timepoint' %in% colnames(array_design)){
      timepoint <- levels(array_design$Timepoint)
      if(length(timepoint) >1){ timepoint <- ''}
  }else{
      timepoint <- ''
  }
  
  print(paste0('timepoint is ', timepoint))
  
  if(timepoint !=''){
      plot_pre <- paste0(plot_pre, timepoint, '_')
      result_pre <- paste0(result_pre, timepoint, '_')
      dataset_p<- paste0(dataset, '(', timepoint,')')  ## add the timepoint to the plot title (not in the saved plotname)
  }else{
      dataset_p<- dataset
  }
 
  ###############################
  ###### if combine all mouse models into 1 genotype
  ###############################
  if(combined_genotype){
      array_design$genotype_original = array_design$Genotype
      array_design$Genotype = as.character(array_design$Genotype)
      array_design$Genotype[which(array_design$Genotype != 'WT')] = disease ## assign all mouse models as the same genotype (disease)
      array_design$Genotype = factor(array_design$Genotype, levels = c('WT', disease))
      }

  ###############################
  ###### relevel the factors
  ###### and assign factors for DEA
  ###############################
  # relevel the factors as marked in the model file

  x <- as.character(model_info$factor_levels)
  if(!is.na(x)){
    print('Relevel factors')
    factor_relevels <- unlist(strsplit(x, split = '; '))
    for (i in factor_relevels){
      eval(parse(text =i))
    }
  }
  
  ## get the subset model list from file
  assign('factor_ls', eval(parse(text = as.character(model_info$factors))))
  
  ############################################
  ## model limma based on the given factors
  ############################################
  (model_formular <- paste0('model.matrix( ~ ', paste0(factor_ls, collapse=' + '), ', data = array_design)' ))
  
  msg_model <- paste0('\n\n#-------------------------------#\n MODEL INFO\n#-------------------------------#\n',
                      'Factors for limma analysis: ', paste0(factor_ls, collapse = ", "),
                      '\nModel: ', model_formular)
  
  cat(msg_model)
  
  assign('design_matrix', eval(parse(text = model_formular)))
  array_fit <- lmFit(array_dat, design_matrix)
  array_eb_fit <- eBayes(array_fit)
  
  (coef_ls <- colnames(array_eb_fit$coefficients))
  gene_annotation <- annotation[, c('ProbeName', 'GeneSymbols',  'GeneNames')]
  
  ########################################
  # loop for each contrast coef:
  #' 1. plot p value density, save top tables, save toptables after filter
  #######################################
  msg_coef_all <-"" 
  
  ## make a index for coef match to which factor (one factor can have multiple coefs)
  factor_index <- NULL
  for(i in 2:length(coef_ls)){
    y <- coef_ls[i]
    for(j in 1:length(factor_ls)){
      x <- factor_ls[j]
      (to_match <- substr(y, 1, nchar(x)))
      if(sum(grepl(x, to_match)) >0){
        factor_index <- c(factor_index, j)
      }    
    }
  }
  print("process top tables")
  for (i in 2:length(coef_ls)){
    #*****************************#
    # 1. get toptable of the coef and process
    #*****************************#
    ## get the baseline level
    target_factor <- factor_ls[factor_index[i-1]]
    print(paste0('target factor', target_factor))
    eval(parse(text = paste0('baseline <- ', 'as.character(levels(array_design$', target_factor, ')[1])')))

    
    ## reorder the array_dat samples by target_factor
    cmd =paste0('sample_order <- as.character(array_design[with(array_design, order(', target_factor, ')), "Sample"])')
    print(cmd)

    eval(parse(text = cmd))
    # print(sample_order)
    print(setdiff(colnames(array_dat), sample_order))
    array_dat<-array_dat[, sample_order]

    
    ## coef
    target_coef <- coef_ls[i]
    array_toptable <- topTable(array_eb_fit, 
                               coef = target_coef,
                               number = Inf) # get all the probes
    
    ## process top table
    toptable_df <- cbind(rownames(array_toptable),
                         array_toptable[, c('P.Value', 'adj.P.Val', 'logFC', 't')])
    table_title <- c('ProbeName', 'pValue', 'qValue', 'LogFC','TestStat')
    colnames(toptable_df) <- table_title
    ## add gene symbol and gene names
    toptable_df <- noWarnings(left_join(toptable_df, gene_annotation))
    
    #*****************************#
    # 2. split limma 2-sided pvalue into 1-sided pvalue of up and down-regulated genes and plot p values
    #*****************************#
    ## split the result for one sided p-value(not FDR) based on test stat (p-value for up or)
    ## - to avoid when a gene is up in one study and down in another
    ## 1-sided pvalues and adj pvalues are added to the toptable
    
    ## down regulated genes(for each probe) and BH adjust the 1 sided p value
    splitPvalUp <- function(pval, teststat){
      # 1 sided pvalue for upregulation probes
      # return 1 sided p value
      if(teststat>0){
        temp <- pval/2
      }else{
        temp <- 1-pval/2
      }
      return(temp)
    }
    
    splitPvalDown <- function(pval, teststat){
      # 1 sided pvalue for downregulation probes
      # return 1 sided p value
      if(teststat < 0){
        temp <- pval/2
      }else{
        temp <- 1-pval/2
      }
      return(temp)
    }
    
    print(paste0('Calculating 1-sided p-values for ', target_coef))
    toptable_df$up <- NA
    toptable_df$down <- NA
    for (i in 1:nrow(toptable_df)){
      toptable_df$up[i] <- splitPvalUp(toptable_df$pValue[i], toptable_df$TestStat[i])  
      toptable_df$down[i] <- splitPvalDown(toptable_df$pValue[i], toptable_df$TestStat[i])  
    }
    
    ## use adj_method to correction for the 1 sided pvalues
    toptable_df$up_padj <- p.adjust(toptable_df$up,adj_method)
    toptable_df$down_padj <- p.adjust(toptable_df$down,adj_method)
    
    ## plot 2-sided p values, and 1-sided up and down p values
    for (i in c('pValue', 'up', 'down')){
      print(paste0("Plot ", i, " p values"))
      df <- as.data.frame(toptable_df[, i])
      colnames(df) <- target_coef
      
      ## change target name for ggplot (to qoute the target, so that ggplot works for special chars)
      (tmp_target <-paste0("`", target_coef,"`" ))
      
      ## P value histo plot for individual factors
      filename <- paste0(plot_pre,i,'_', target_coef, '-vs-baseline_', baseline,'.png')
      (plt_title <- paste0(dataset_p, ' ', i, '\n', target_coef))
      p <- ggplot(df, aes_string(x = tmp_target)) + 
          geom_histogram(alpha = 0.5, binwidth = 0.005) +
          theme_bw() +
          ggtitle(plt_title)
          
      ggsave(filename = filename, plot=p, width = 6, height =5, units = "in")  
    }
    
    #*****************************#
    # 3. save results of DE genes, and 1 sided pvalues
    #*****************************#
    # - 3.1_save DE genes (2-sided p values) with FDR under threshold, and plot heatmap of the expression of these DE genes
    # - 3.2_add comment log for the tables to be saved
    # - 3.3_save the 1 sided p values (not adjusted) of all probesets (for genotype and other factors)
    # - 3.4_save the 1 sided p values (not adjusted) of all genes (for genotype only) for metaanalysis
    # - 3.5 save the toptable, which contains 2-sided pvalues, 1-sided pvalues, adj 1-sided p values, fdr, logFc, t etc
    
    ## get the timepoint if available
    if("Timepoint" %in% colnames(array_design)){
        timepoint <- levels(array_design$Timepoint)
        if(length(timepoint) >1){
            print(paste0("more than 1 timepoint is included. Use the first timepoint: ", levels(array_design$Timepoint)[1]))
            timepoint <- levels(array_design$Timepoint)[1]
        }
    }else{
        timepoint <- ""
    }
    
    ###
    ### 3.1_save DE genes (2-sided p values) with FDR under threshold, and plot heatmap of the expression of these DE genes
    ###
    df_t <- toptable_df[which(toptable_df$qValue <= q_threshold), ]
    
    ## save DE under threshold
    if(write_table){
        file_name <- paste0(result_pre, 'limma_', target_coef, '-vs-baseline_', baseline,'-q',q_threshold, '_', timepoint,'.tsv')
        write.table(df_t, file = file_name, row.names = F,
                    sep ='\t', quote = F)
    }
    
    ## get the expression of DE probes
    probe_names <- as.character(df_t$ProbeName)
    DE_df <- array_dat[probe_names, ] %>% droplevels()
    
    ## up and down regulated
    gene_up <- df_t[which(df_t$LogFC >0), "ProbeName"]
    df_up <- array_dat[gene_up, ] %>% droplevels()
    
    gene_down <- df_t[which(df_t$LogFC <0), "ProbeName"]
    df_down <- array_dat[gene_down, ] %>% droplevels()
    
    if(plot_heatmap){
      print(paste0('Plot heatmaps for DE genes:', target_coef, ', q_threshold: ', q_threshold))
      heatmapForTop200(DE_df,width, height,gene_annotation, array_design, target_coef, baseline,plot_pre, q_threshold,prefix = "all", timepoint=timepoint)
      heatmapForTop200(df_up, width, height,gene_annotation, array_design,target_coef, baseline,plot_pre,q_threshold,prefix = "up", timepoint=timepoint)
      heatmapForTop200(df_down, width, height,gene_annotation, array_design,target_coef, baseline,plot_pre,q_threshold,prefix = "down", timepoint=timepoint)
    }
    
    ###
    ### 3.2_add comment log for the tables to be saved
    ###
    ## log
    (msg_coef <- paste0('Factor: ', target_factor, '; baseline: ', baseline, 
                        ', coefficient: ', target_coef, 
                        '\n  DE genes: ', nrow(df_t), ', up: ', nrow(df_up), 
                        ', down: ', nrow(df_down), '\n'))
    msg_coef_all <-paste0(msg_coef_all, msg_coef)
    msg_pval_adj <- paste0("")
    #msg_pval_adj <- paste0('\n 1-sided pValue is adjusted by ', adj_method)
    
    ###
    ### 3.3_save the 1 sided p values (not adjusted) of all probesets (for genotype and other factors)
    ###
    
    ## Save the up and down pvalue for all probesets

    ## message headers
    tmp_msg <- paste0(preprocess_msg,'\n#-------------------------------#',
           '\n## LIMMA DE ANALYSIS\n#-------------------------------#\n',
           Sys.Date(), '----Dataset: ', dataset, '\n', 
           'Imported R objects: ', f_data, '\n',
           msg_subset, msg_model, '\n', ncol(array_dat), ' samples, ', nrow(array_dat), ' probes',
           "\n", "coefficients: ", paste(coef_ls, collapse = ", "),
           "\n", "q value threshold:", q_threshold, '\n',
           msg_coef)
    
    tmp_msg <- paste(unlist(strsplit(tmp_msg, split = '\\\n')), collapse="\n#")
    tmp_msg <- paste0("#", tmp_msg)
    
    
    if(write_table){
    # order by 1-sided p value and append table UP REGULATED
    file_name <- paste0(result_pre, 'probeset_up_regulation_1_sided_pval_', target_coef, '-vs-baseline_', baseline, '_', timepoint,'.tsv')
    # write the comments
    sink(file_name, type="output")
    writeLines(tmp_msg)
    sink()
    tmp_df_up <- toptable_df[order(toptable_df$up),c('ProbeName', 'GeneSymbols', 'up')]
    noWarnings(write.table(tmp_df_up, file = file_name, row.names = F,
                           sep ='\t', quote = F, append = T))
    
    # order by 1-sided p value and append table DOWN REGULATED
    file_name <- paste0(result_pre, 'probeset_down_regulation_1_sided_pval_', target_coef, '-vs-baseline_', baseline, '_', timepoint,'.tsv')
    # write the comments
    sink(file_name, type="output")
    writeLines(tmp_msg)
    sink()
    tmp_df_down <- toptable_df[order(toptable_df$down), c('ProbeName', 'GeneSymbols', 'down')]
    noWarnings(write.table(tmp_df_down, file = file_name, row.names = F,
                           sep ='\t', quote = F, append = T))
    }
    
    
    ###
    ### 3.4_save the 1 sided p values (not adjusted) of all genes (for genotype only) for metaanalysis
    ###
    # - if multiple probes are mapped to the same gene, the p values of these probesets are adjusted by adj_method (default is BH)
    #   and the minimum p is selected
    if(write_table){if(grepl("Genotype", target_coef)){
        print(paste0("Aggregate to gene p values: ", target_coef))
        timepoint <- as.character(levels(array_design$Timepoint))
        if(length(timepoint >1)){
            timepoint <- ''
        }else{timepoint <- paste0("_",timepoint)}
        
        ## for UP regulation
        x <- aggregateGenePvals(tmp_df_up, dataset, timepoint, target_coef, adj_method, plot = F, verbose =F)
        file_name <- paste0(result_pre, 'gene_up_regulation_1_sided_pval_', target_coef, '-vs-baseline_', baseline, '_', timepoint,'.tsv')
        df1<- x[[1]]
        msg_meta <- x[[2]]
        tmp_msg_meta <- paste0(tmp_msg, "\n#", msg_meta)
        # write the comments
        sink(file_name, type="output")
        writeLines(tmp_msg_meta)
        sink()
        noWarnings(write.table(df1, file = file_name, row.names = F,
                               sep ='\t', quote = F, append = T))
        
        ## for DOWN regulation
        x <- aggregateGenePvals(tmp_df_down, dataset, timepoint, target_coef, adj_method, plot = F, verbose =F)
        file_name <- paste0(result_pre, 'gene_down_regulation_1_sided_pval_', target_coef, '-vs-baseline_', baseline, '_', timepoint,'.tsv')
        df1<- x[[1]]
        msg_meta <- x[[2]]
        tmp_msg_meta <- paste0(tmp_msg, "\n#", msg_meta)
        # write the comments
        sink(file_name, type="output")
        writeLines(tmp_msg_meta)
        sink()
        noWarnings(write.table(df1, file = file_name, row.names = F,
                               sep ='\t', quote = F, append = T))
    }}
    
    
    ###
    ### 3.5 save the toptable, which contains 2-sided pvalues, 1-sided pvalues, adj 1-sided p values, fdr, logFc, t etc
    ###
    if(write_table){
        file_name <- paste0(result_pre, 'limma_toptable_', target_coef, '-vs-baseline_', baseline, '_', timepoint,'.tsv')
        # write the comments
        sink(file_name, type="output")
        writeLines(tmp_msg)
        sink()
        # write table
        noWarnings(write.table(toptable_df, file = file_name, row.names = F,
                               sep ='\t', quote = F, append = T))}
  } # end of loop for each coef
  
  ########################################
  #  Save R objects?? - may not be necessary
  #######################################
  ## save logs as R objects
  msg_limma <- paste0(preprocess_msg,'\n#-------------------------------#',
                      '\n## LIMMA DE ANALYSIS\n#-------------------------------#\n',
                      Sys.Date(), '----Dataset: ', dataset, '\n', 
                      'Imported R objects: ', f_data, '\n',
                      msg_subset, msg_model, '\n', ncol(array_dat), ' samples, ', nrow(array_dat), ' probes',
                        "\n", "coefficients: ", paste(coef_ls, collapse = ", "),
                       "\n", "q value threshold:", q_threshold, msg_pval_adj,
                      '\n\n#-------------------------------#\n DE ANALYSIS\n#-------------------------------#\n',
                       msg_coef_all)
  #print(msg_limma)
  
  ########################################
  #  write log
  #######################################
  sink(paste0(result_pre, "limma_log.txt"), type="output")
  writeLines(msg_limma)
  sink()
  
  ## save r objects
  save(array_dat, array_design, annotation, array_eb_fit, dataset, 
       gene_annotation, msg_limma, design_matrix, platform,model_info,model_formular,
       file = paste0(result_pre, "limma_objects.Rdata"))
  print(paste0("R objects saved, ", paste0(result_pre, "limma_objects.Rdata")))
  
  ## save a stat (table log), have a over all summary and also indi for each genotype
  df_stat <- data.frame(dataset = dataset, platform = platform, probes = nrow(array_dat), 
             genes = length(unique(annotation$GeneSymbols)),
             sample_size_all = nrow(array_design),
             batch_correction = model_info$batch_corrected,
             expr_min = min(array_dat),
             expr_max = max(array_dat),
             expr_median = median(as.matrix(array_dat)),
             limma_model = model_formular)
  
     for (i in c("Genotype","Timepoint", "OrganismPart", "Gender")){
         if (i %in% colnames(array_design)){
             (tmp <- table(array_design[, i]))
             tmp <- paste0(paste0(names(tmp), "(", tmp, ")"), collapse="; ")
             eval(parse(text = paste0("df_stat$", i , "_all <- tmp")))
         }else{
             eval(parse(text = paste0("df_stat$", i , "_all <- NA")))
         }
     }
  ## for indi genotypes
  case_geno <- setdiff(levels(array_design$Genotype), 'WT') ## mutant genotypes
  df_stat_indi <- NULL
  for (case_geno_i in case_geno){
      array_design_indi <- filterContain(array_design, 'Genotype', paste0(case_geno_i, '|WT'))
      ## count the genotype numbers
      (tmp <- table(array_design_indi[, 'Genotype']))
      control_n <- as.numeric(tmp['WT'])
      case_n <- as.numeric(tmp[case_geno_i])
      
      if ('Gender' %in% colnames(array_design_indi)){
          (tmp <- table(array_design_indi[, 'Gender']))
          gender_count <- paste0(paste0(names(tmp), "(", tmp, ")"), collapse="; ")
      }else{
          gender_count <- NA
      }
      if ('Timepoint' %in% colnames(array_design_indi)){
          time_count <- levels(array_design_indi$Timepoint)
      }else{
          time_count <- NA
      }
      
      df_stat_tmp <- data.frame(dataset  = dataset, 
                                Genotype = case_geno_i,
                                Timepoint = time_count,
                                control_n = control_n,
                                case_n    = case_n,
                                gender    = gender_count)
      df_stat_indi <- rbind(df_stat_indi, df_stat_tmp)
  }
  
  df_stat <- noWarnings(left_join(df_stat_indi, df_stat))
  
  
    file_name <- paste0(result_pre, "limma_table_log.tsv")
    sink(file_name, type="output")
    writeLines(paste0("#", Sys.Date(), '----Dataset: ', dataset))
    sink()
    # write table
    noWarnings(write.table(df_stat, file = file_name, row.names = F,
                         sep ='\t', quote = F, append = T))
}


#---------------------------------------------------------------------------#
# Part 5: LOOP
#---------------------------------------------------------------------------#

## load model info
model_info_all <- read.delim(f_model,check.names = F, comment.char = "#")

## default indecies
if(!exists("start_row")){start_row <- 1}
if(!exists("to_do")){to_do <- c(start_row:nrow(model_info_all))}
if(exists("dataset_todo")){
  to_do <- which(model_info_all$eeName %in% dataset_todo)}


# 1. set a start_row to loop the experiment list, default start from first row
if(!exists("start_row")){start_row <- 1}
# 2. can select a list of experiments by row index (over-write default to_do) 
if(!exists("to_do_ls")){to_do <- c(start_row:nrow(model_info_all))}
# 3. pick a specific dataset by GSE id (over-write default to_do) 
if(!exists("to_do_ls") & exists("dataset_todo")){
    to_do <- which(model_info_all$eeName %in% dataset_todo)}
if(exists("to_do_ls")){
    to_do <- to_do_ls}



for (i in to_do){
  print(i)
  dataset <- as.character(model_info_all$eeName[i])
  model_info <- model_info_all[i, ]%>% droplevels()  # subset the model_info
  
  ## load the R objects (no-batch correction, or batch corrected)
  if (tolower(as.character(model_info$batch_corrected)) == "yes"){
    ## array with batch corrected (outlier removed, subset of samples if previously processed)
    print(paste0("Loading ", dataset, " , batch corrected"))
    f_data <- paste0(object_dir, dataset, '/results/', dataset, '_objects_after_batch_correction.Rdata')
  } else {
    ## array without batch corrected (outlier removed, subset of samples if previously processed)
    print(paste0("Loading ", dataset, " , not batch corrected"))
    f_data <- paste0(object_dir, dataset, '/results/', dataset, '_objects.Rdata')
  }
  load(f_data)
  print(paste0('Loading R objects, ', f_data))
  mainLimmaDEA(array_dat, array_design, annotation, model_info, q_threshold,preprocess_msg = msg, adj_method=adj_method,
               plot_heatmap=plot_heatmap,combined_genotype = combined_genotype, write_table = write_table)
}
  