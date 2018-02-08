#' created 2015-11-?
#' update 2016-05-19
#' ################################################################################################
#' helper functions for 
#' combine_experimental_design.R
#' explore_data_general.R
#' limma_DEA.R
#' and others
#' 
#' dependencies
#' library('HelperFunctions') # for table functions and heatmap
#' Function list:
#' printSampleCountsByFactor() #OUTPUT: a msg (str) for log: count of sample size per factor level
#' 
#' plotHeatmap() # general function to produce heatmap
#' plotGeneCor() # plot gene-gene correlation (on randomly selected 3k probes)
#' 
#' plotFunctions() # produce sample corelation heatmap and PCA plots, gene-gene correlation plots
#'
#' modelTest() # A stat test: if the array matrix is associated with the response (e.g. different geneypes, batches etc)
#' rmSamples()  # remove a sample by name in data matrix and design file
#' subsetSamples() # subset data matrix and design file by given factors and levels to keep
#' 
#' getGemmaLogInfo() # History and Technical log files of Gemma uploads
#' 
#' processDesign() # get the specified dataset design from a file contain all designs of all datasets
#' processDesignOriginal() # read design file individually, used for combining multiple design files
#' processPlatformAnnotation() # read and process platform annotation
#' 
#' sanityCheckHeatmap() # plot the heatmap for selected genes for sanity check, sample are reordered by genotype (if recorded)
#' topGeneHeatmap() # plot the heatmap for a given gene list for input study
#' 
#' ################################################################################################


## load packages
package_ls <- c('plyr', 'dplyr','ggplot2', 'limma', 'RColorBrewer', 'reshape2', 'globaltest', 'grid', 'scales','gtable', 'VennDiagram','pheatmap')
for(i in package_ls){
  cmd <- paste0("if('package:", i, "' %in% search() == FALSE){library(", i, ")}")
  eval(parse(text =cmd))
}

rm(package_ls)
rm(cmd)
rm(i)

## source 
 
library('HelperFunctions')

## default values that needs update
# x_col: factor columns in design file
x_col <- c("Genotype", "Treatment",  "OrganismPart", "Organism.part", "Agent", 
           "DevelopmentalStage", "SamplingTimePoint", "Treatment.1", 
           "Collection.of.material", "Phenotype", "Batch","Gender","Timepoint",
           "GeneticBackgroud", "Sample",
           "Study", "Disease_stage")
## do not include "Original_genotype", GSE64398.1 will have trouble batch correct



#**************************************************************************#
#**************************************************************************#
#********** handy helper functions
#**************************************************************************#
#**************************************************************************#


saveSampleInfo <- function(file_dir, design_f){
    # for data preprocessing
    #' save a SampleInfo.txt with the first column list of  all CEL files in the dir and the rest are extra info from the 
    #' expermental design
    #' 
    #' design_f: design file for all datasets
    #' NOTES:
    #' cannot have any white space in the contents
    
    # grep files
    (file_ls <- grep('GSM.*\\.CEL', list.files(file_dir), value=T))  # grep CEL files only
    # if no CEL files, try to grep GSM
    if(length(file_ls) == 0){
        (file_ls <- grep('^GSM', list.files(file_dir), value=T))
    }
    if(length(file_ls) == 0){
        (file_ls <- grep('.CEL', list.files(file_dir), value=T, ignore.case = T))
    }
    
    (file_ls <- data.frame(FileName = file_ls))
    file_ls$ExternalID <-grep("^GSM", unlist(strsplit(as.character(file_ls$FileName), split = "\\.|\\_")), value = T) ## for files start with GSM
    
    
    ## get extra info from gemma design files
    (dataset <- getGSEID(file_dir))
    array_design <- processDesign(dataset, design_f)
    #df <- noWarnings(left_join(file_ls, array_design[, c("ExternalID", "Sample")]))
    df <- noWarnings(left_join(file_ls, array_design))
    df$label <- paste0(df$ExternalID,",", df$Sample)  
    
    sample_file <- paste0(file_dir, "/SampleInfo.txt")
    write.table(df, file = sample_file, row.names = F, sep ='\t', quote = F, col.names=T)
    print(paste0("return file ", sample_file))
}


getGSEID <- function(file_dir){
    ## get the GSE(str) from the file dir,for label purpose, only the GSE1234 str
  (dataset <- grep("^GSE[0-9]+$|^GSE.*[0-9]+$|^GSE.*l$|^GSE.*e$", unlist(strsplit(as.character(file_dir), split = "/")), value= T))
  if(length(dataset) == 0){
    (dataset <- grep("^GSE", unlist(strsplit(as.character(file_dir), split = "/")), value= T))
    if(length(dataset) > 1){dataset <- dataset[1]}
  }
    return(dataset)
}


printSampleCountsByFactor <- function(array_design){
  #' INPUT:
  #'  - array_design(df)
  #' OUTPUT: a msg (str) for log: count of sample size per factor level
  #' 
  design_factors <- intersect(x_col, colnames(array_design))  # design_factors(char list): a list of factors
  msg2 <- NULL
  for(i in design_factors){
    if(i != "Sample"){
      (temp_m <- levels(array_design[, i]))
      x <- table(array_design[, i])
      x <- paste(names(x), x, sep=": n=")
      (msg_t <- paste0("Factor: ", i, " (", length(temp_m), " levels): ", paste0(x,collapse =", ")))
      
      msg2 <- paste0(msg2, msg_t, collapes ="\n")
    }
  }
  return(msg2)
}

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

#**************************************************************************#
#**************************************************************************#
#********** heatmap function
#**************************************************************************#
#**************************************************************************#
plotHeatmap <- function(x, title = "", size_r = 5, size_c = 15, 
                        show_rownames = T,
                        show_colnames = T,
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE,
                        annotation_row = F,
                        array_design_df = array_design,
                        color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
                        default_border_color = NA,
                        NA_color = "grey",
                        scale_row = F,
                        legend = NULL, ...){
  ## input an ordered array matrix and plot the heatmap. heatmap function is 
  ## INPUT:
  ##   - x: is an ordered matrix (the col names matched with row names of the design data )
  ##   - array_design: must be defined for legend, row names are matched with design matrix column names
  ##   - show_col, or rownames (default = T): whether to display sample names in column and rows
  ##   - title(string): plot title
  ##   - size_r, size_c are the font size for row and column
  ##   - annotation_row (default = F), whether to show annotation bar on the side (annotation bar on the top is always shown)
  ##   - color_scheme(string): choose a brewer scheme, e.g. 'BuPu', 'YlOrRd', default is 'YlOrRd'
  ##   - default_border_color (default = NA): whether to have a border drawn for the cells, 
  ##       default is NA (no borders), pheamap default is border_color = "grey60"
  ##   - NA_color: define color for NA value (in pheatmapNA() function)
  ##   - scale_row: use z score to normalize the rows (default = F)
  ##   - ... more options in pheatmapNA(), e.g. breaks, scale etc.
  ## OUTPUT:
  ##   - a pretty heatmap
  
  #heatmap set legend
    row.names(array_design_df) <- array_design_df$Sample ## must set the rownames same as the input matrix for heatmap to map the annotation
    if(is.null(legend)){
        (legend <- setdiff(colnames(array_design_df), c("Bioassay","ExternalID", "Sample", "Dataset", "ExpMatrixName")))
    }
  (annotation <- array_design_df[legend])
    ## for annotation, 

    if(annotation_row){
        annotation_col = annotation
        annotation_row = annotation
        annotation = NA
    } else{
        annotation_col = NA
        annotation_row = NA
        annotation = annotation
    }
    if(scale_row){
        title <- paste(title," (z-score)")
        scale = 'row'
    } else{
        scale = 'none'
    }
    pheatmap(x,color = color,
             na_col = NA_color,
             cluster_rows = cluster_rows, 
             cluster_cols = cluster_cols,
             fontsize_row = size_r,
             fontsize_col = size_c,
             show_rownames = show_rownames,
             show_colnames = show_colnames,
             main = title,
             border_color = default_border_color,
             annotation = annotation,
             annotation_col = annotation_col,
             annotation_row = annotation_row,
             scale = scale,
             ...)
#   if(annotation_row){
#     pheatmapNA(x, colour_scheme = colour_scheme,
#                NA_color = NA_color,
#              cluster_rows = cluster_rows, cluster_cols = cluster_cols, # turn on or off dendrogram
#              annotation_col = annotation,
#              annotation_row = annotation,
#              fontsize_row = size_r,
#              fontsize_col = size_c,
#              show_rownames = show_rownames,
#              show_colnames = show_colnames,
#              main = title,
#              border_color = default_border_color, ...)
#   } else{
#     pheatmapNA(x, colour_scheme = colour_scheme,
#                NA_color = NA_color,
#            cluster_rows = cluster_rows, cluster_cols = cluster_cols, # turn on or off dendrogram
#            annotation = annotation,
#            fontsize_row = size_r,
#            fontsize_col = size_c,
#            show_rownames = show_rownames,
#            show_colnames = show_colnames,
#            main = title, 
#            border_color = default_border_color, ...)
#   }
}




#**************************************************************************#
#**************************************************************************#
# plot gene-gene correlation
#**************************************************************************#
#**************************************************************************#
plotGeneCor <- function(array_dat, dataset){
  ## sample 3k probes, or all probes if <3k
  ## plot correlation density of these selected genes
  index <- sample(1:nrow(array_dat), min(3000, nrow(array_dat)))
  df<- array_dat[index, ]
  df_cor <- cor(t(df))
  d <- density(df_cor, na.rm = T) # returns the density data
  plot(d, main = paste0("gene-gene correlation of ", 
                        dataset, "(", ncol(array_dat), " samples)"))
  abline(v=0, lty=2)
}


#**************************************************************************#
#**************************************************************************#
# plot heatmap and PCA, gene-gene correlations
#**************************************************************************#
#**************************************************************************#

plotFunctions <- function(array_dat, plot_pre, dataset, 
                          array_design,prefix ="", gene_cor =T, width = 1600,
                          height = 1000,simple_heatmap=F, ...){
  ## plot sample correlation heatmaps and PCA
  ## INPUT:
  ##   -array_dat (dataframe) is processed array_dat
  ##   -plot_pre (string): prefix of the plot dir
  ##   -dataset: experiment id (string). e.g. GSE1556
  ##   -array_design is the processed metat data, must match with the sample numbers of array_dat
  ##   -prefix(string): for the plot titles only
  ##   -gene_cor; output the gene-gene correlaiont plot (default = T)
  ##   - simple_heatmap: if just output a simple heatmap without column labels
  ##   -... is extra parameters for plotHeatmap() 
  ## OUTPUT: expression density for samples, heatmaps, PCA plots, gene-gene-correlation plots
    
    #'@notes: array design must has rownames that match with array_dat sample names
  
  ###############################
  ###### expression density
  ###############################
  main_title <- paste0('Sample Density: ', dataset, 
                       '\n', ncol(array_dat), ' samples, ', nrow(array_dat), ' probes')
  png(filename= paste0(plot_pre,'exp_density_', ncol(array_dat), "_samples.png"), width = width, height = height)
  plotDensities(array_dat, legend=F, main = main_title)
  dev.off()
  ###############################
  ###### heatmaps, correlations
  ###############################
  ## reorder samples by x (for heatmap)
  x <- c('Genotype', 'OrganismPart', 'Timepoint', 'Gender','Study','Disease_stage')
  x <- intersect(colnames(array_design), x)
  
  if(all(c('Study','Disease_stage') %in% x)){ # this is for mega analysis of all sample from all studies
      x <- c('Disease_stage', 'Study',setdiff(x, c('Study','Disease_stage')))
  }
  
  if(length(x) != 0){
    # get the sample names in the order by elements in x
    (cmd <- paste0("sample_order <- as.character(array_design[with(array_design, order(", 
                   paste(x, collapse=', '),
                   ")), 'Sample'])"))
    eval(parse(text =cmd))
    array_dat <- array_dat[, sample_order]
    
    ## reorder the design rows too
    (cmd <- paste0("array_design <- array_design[with(array_design, order(", 
                   paste(x, collapse=', '),
                   ")), ]"))
    eval(parse(text =cmd))
  }
  
  ## get sample to sample correlation, diag cells should be coloured as missing value
  sample_cor <- cor(array_dat)
  ## replace diag cells with NA value
  diag(sample_cor) <- NA

  if(simple_heatmap){ ## simple heatmap
      png(filename= paste0(plot_pre, 'sample_corr.png'), width = width, height = height)
      pheatmapNA(sample_cor, default_border_color = NA)
      dev.off()
  }else{ #complex heatmap (default)
      ## save plot
      png(filename= paste0(plot_pre, 'sample_corr.png'), width = width, height = height)
      plotHeatmap(sample_cor, paste0("Pearson correlation of samples_",dataset, prefix), 
                  size_r = 15, size_c = 15,
                  cluster_rows = F, 
                  cluster_cols = F,
                  annotation_row = T, 
                  array_design_df = array_design)
      dev.off()
      
      ## plot a heatmap by batch if batch is in the factors
      if("Batch" %in% colnames(array_design)){
          sample_order <- as.character(array_design[with(array_design, order(Batch)), "Sample"])
          array_dat_tmp<-array_dat[, sample_order]
          sample_cor <- cor(array_dat_tmp)
          diag(sample_cor) <- NA
          png(filename= paste0(plot_pre, 'sample_corr_order_by_batch.png'), width = width, height = height)
          plotHeatmap(sample_cor, paste0("Pearson correlation of samples (by batch)_",dataset, prefix), size_r = 15, size_c = 15,
                      cluster_rows = F, 
                      cluster_cols = F,
                      annotation_row = T, 
                      array_design_df = array_design, ...)
          dev.off()
      }
      
      ## plot a heatmap by Study if study is in the factors (a mixed model)
      if("Study" %in% colnames(array_design)){
          sample_order <- as.character(array_design[with(array_design, order(Study)), "Sample"])
          array_dat_tmp<-array_dat[, sample_order]
          sample_cor <- cor(array_dat_tmp)
          diag(sample_cor) <- NA
          png(filename= paste0(plot_pre, 'sample_corr_order_by_Study.png'), width = width, height = height)
          plotHeatmap(sample_cor, paste0("Pearson correlation of samples (by Study)_",dataset, prefix), size_r = 15, size_c = 15,
                      cluster_rows = F, 
                      cluster_cols = F,
                      annotation_row = T, 
                      array_design_df = array_design, ...)
          dev.off()
      }
      
  }

  
  
  ###############################
  ###### PCAs
  ###############################
  #prcomp() use the singular value decomposition (SVD).
  ## remove constant/zero column to unit variance
  row_variance <- apply(array_dat, 1, var)
  index <- which(row_variance ==0)
  if(length(index) >0){
      df <- array_dat[-index, ]
  }else{
      df <- na.omit(array_dat)
  }

  
  array_pca <- prcomp(t(na.omit(df)), scale = T)  ## must have a matrix with probes as columns, and samples as rows
  #quick plot for each loadings # normally PC1 should not > 0.8, a high PC1 may indicate an strong outliers (which
  # seperate itself from the rest of the group)
  proportion_of_variance <- summary(array_pca)$importance[2, ]
  png(filename= paste0(plot_pre, 'PCA.png'), width = width, height = height)
  barplot(proportion_of_variance, main = dataset, ylab = "proportion of variance")
  dev.off()
  
  
  ## get PC1 and PC2 (u) values for each sample, the formular is taken from 
  #'https://github.com/vqv/ggbiplot/blob/master/R/ggbiplot.r
  #'
#   pcobj <- prcomp(t(array_dat), scale = T)
#   (nobs.factor <- sqrt(nrow(pcobj$x) - 1))
#   (d <- pcobj$sdev)
#   (u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*'))
#   df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*')) #contain choices=1:2 PC of each sample
  
  (nobs.factor <- sqrt(nrow(array_pca$x) - 1))
  (d <- array_pca$sdev)
  (u <- sweep(array_pca$x, 2, 1 / (d * nobs.factor), FUN = '*'))
  (v <- array_pca$rotation)

  pc_n = 1:4      # choose the first few pcs (must <= the number of samples)
  obs.scale = 1
  df.u <- as.data.frame(sweep(u[, pc_n], 2, d[pc_n]^obs.scale, FUN='*')) # datafame of the first 2 PCs for the samples

  #append rotations to the array_design
  rownames(array_design) <- array_design$Sample
  array_prin_comp <- mapBindAllColumns(array_design, df.u)
  # set the color group (Genotype as default)
  color_col <- intersect(x_col, colnames(array_prin_comp))
  
  # plot the first two PCs color by each factor
  # labels are right aligned hjust = 0 does what you want. (hjust = 0 left-justified, 0.5: centered, 1: right-justified.)
  for (i in color_col){
    if(i != "Sample"){
      ggplot(array_prin_comp, aes_string("PC2","PC1", label = "Sample", color = i, hjust = 1.1)) +
        geom_point() +
        theme_bw() +
        xlab(paste0("PC2 (", proportion_of_variance[2] *100, "%)"))+  #add proportion of variance for PC2
        ylab(paste0("PC1 (", proportion_of_variance[1] *100, "%)"))+ #add proportion of variance for PC1
        geom_text(size=2) +
        ggtitle(paste0("First two principal components: ", dataset, prefix)) 
      ggsave(filename = paste0(plot_pre, 'PCA_colored_by_',i, '.png'), width = 8, height=5)
    }
  }
  
  ## plot the PCA1, 2,3, for the factors, require(rshape2)
  df <- array_prin_comp[, c(color_col, "PC1", "PC2", "PC3")] %>% droplevels()
  df <- melt(df, id=c("PC1", "PC2", "PC3", "Sample"))
  colnames(df)[5:6]<- c("factor", "factor_levels")
  for (i in c("PC1", "PC2", "PC3")){
    title <- paste0("")
    ggplot(df, aes_string("factor_levels", i, col = "factor")) + 
      geom_point()+
      theme_bw() +
      ggtitle(paste0(i, " and factors: ", dataset, prefix)) +
      theme(title = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1)) 
    ggsave(filename = paste0(plot_pre, 'PCA_',i, '_and_facotrs.png'), width = 8, height=5)
  }
  ###############################
  ###### gene-gene correlation
  ###############################
  if(gene_cor){
    png(filename= paste0(plot_pre, 'gene_corr.png'), width = width, height = height)
    plotGeneCor(array_dat, dataset)
    dev.off()
  }
  
}

#**************************************************************************#
#**************************************************************************#
#********** Model testing of factor and expression data
#**************************************************************************#
#**************************************************************************#
modelTest <- function(array_dat, array_design, response_ls){
  #' this is a wrapper for the gt function in {globaltest}
  #' A stat test: if the array matrix is associated with the response (e.g. different geneypes, batches etc)
  #' Null hypothesis is no association between response (a factor in the design file) and expression
  #' The test can be used for detecting batch affect. also see ?gt for aims to find which subset of genes is most 
  #' associated with the response
  #' 
  #' INPUT:
  #'  - array_dat(df), array_design(df)
  #'  - response_ls(a list of strings): design factors to be tested one by one
  #' OUTPUT:
  #'  - a summary results (string)
  
  ######## sanity checks ########
  if(class(array_dat) != "data.frame"){
    stop(paste0("Input array_dat is a ", class(array_dat), ". A data.frame is required."))
  }
  if(class(array_design) != "data.frame"){
    stop(paste0("Input array_design is a ", class(array_design), ". A data.frame is required."))
  }
  
  ## check if all input factors are in the design file
  d_ls <- setdiff(response_ls, colnames(array_design))
  has_ls <- intersect(response_ls, colnames(array_design))
  if(length(has_ls) == 0){
    stop(paste0("Input factor list: ", paste0(response_ls, collapse=','), 
                " are not found in the design factors"))
  }
  if(length(d_ls) > 0){
    print(paste0('Factors: ', paste0(d_ls, collapse=', '), 
                 ' are not in the design file. These factors are removed from list'))
    response_ls <- intersect(response_ls, colnames(array_design))
  }
  
  ######## stat test for each factor ########
  alt_array_matrix <-t(as.matrix(array_dat))
  msg <- ""
  for (response in response_ls){
    (cmd<- paste0("factor_levels <- length(levels(array_design$", response,"))")) # check the factor levels
    eval(parse(text = cmd))
    
    
    if(factor_levels == 1){
      print(paste0(response, " contains only 1 level, skipped."))
    }else{
      #tmp <- gt(response = array_design[, response], alternative = alt_array_matrix, trace =F)
      (cmd <- paste0("tmp <- gt(response = array_design[, \"", response, 
                     "\"], alternative = alt_array_matrix, trace =F)"))
      eval(parse(text = cmd))
      
      ## capture the outputs
      #x <- capture.output(summary(tmp))
      x <- t(summary(tmp))
      x <- paste(rownames(x), x[,1], sep = ": ")
      x <-paste(x, collapse="; ")
      x <- paste0('\n## Model test for factor: ', response, '\n', x)
      msg <- paste0(msg, x)
    }
  }
  return(msg)
}


#**************************************************************************#
#**************************************************************************#
#********** remove samples
#**************************************************************************#
#**************************************************************************#
rmSamples <- function(to_rm_ls, array_dat, array_design, to_log =T, rm_notes=""){
  ## remove sample from array design and array dat by a list of sample names
  ## INPUT:
  ##   - to_rm_ls: a list of sample names to be removed from analysis, sep by "; "
  ##   - array_dat(df), array_design(df)
  ##   - notes(string): some notes on why the samples are removed (only for log)
  ## OUTPUT:
  ##   - a new array_dat and array_design and log message
  
  ######## sanity checks ########
  if(class(array_dat) != "data.frame"){
    stop(paste0("Input array_dat is a ", class(array_dat), ". A data.frame is required."))
  }
  if(class(array_design) != "data.frame"){
    stop(paste0("Input array_design is a ", class(array_design), ". A data.frame is required."))
  }
  if(class(to_rm_ls) != "character"){
    stop(paste0("Input samples to be removed is a ", class(to_rm_ls), ". A list of strings is required."))
  }
  
  
  ######## get the rm sample list ########
  (to_rm_ls <- unlist(strsplit(to_rm_ls, split = "; ")))
  (to_rm_notes <- unlist(strsplit(rm_notes, split = "; ")))
  #print(to_rm_notes)
  
  ######## rm samples ########
  ## check if all of the input is in the samples
  missing_sample <- setdiff(to_rm_ls, colnames(array_dat))
  if(length(missing_sample) != 0){
    stop(paste0("sample(s): ", paste(missing_sample, collapse=", "), " are not in the expression data."))
  }
  missing_sample <- setdiff(to_rm_ls, rownames(array_design))
  if(length(missing_sample) != 0){
    stop(paste0("sample(s): ", paste(missing_sample, collapse=", "), " are not in the design file."))
  }
  
  ## return a df of samples removed and remove notes (if specified)
  rm_index <- NULL
  outlier_df <- NULL
  for (i in to_rm_ls){
      if(is.null(rm_index)){
          rm_index <- which(as.character(array_design$Sample) %in% i)
      }else{
          rm_index <- c(rm_index,  which(as.character(array_design$Sample) %in% i))
      }
  }
  
  if(length(rm_index) > 0){
      print('get outlier')
      outlier_df <- array_design[rm_index, ] %>% droplevels()
      outlier_df$category <- 'outlier'
      if(length(to_rm_ls) == length(to_rm_notes)){
          outlier_df$notes <- to_rm_notes
      }
  }
  
  ## rm sample from data matrix
  index <- which(colnames(array_dat) %in% to_rm_ls)
  array_dat <- array_dat[, -index] %>% droplevels()
  
  ## rm sample from design:
  index <- which(rownames(array_design) %in% to_rm_ls)
  array_design <- array_design[-index, ] %>% droplevels()
  

  
  if(!all.equal(rownames(array_design), colnames(array_dat))){
    stop(paste0("Not all samples names matched in both design file and expression matrix, need to check."))
  }
  
  ######### log  ########
  if(to_log){
    if(rm_notes != ""){
      notes <- paste0('\nReason for removal: ', rm_notes)
    }
    msg <- paste0("\n#-------------------------------#\n## REMOVED SAMPLES\n#-------------------------------#\n", 
                  length(to_rm_ls), " sample(s) removed: ", paste0(to_rm_ls, collapse = ", "),
                  notes,
                  "\nnumber of samples after removal: ", ncol(array_dat))
    msg2 <- printSampleCountsByFactor(array_design)

    msg <- paste0(msg, '\n', msg2)
  } else {
    msg <- ""
  }
  
  ######### return dataframe ########
  returnlist <- list(array_dat, array_design, msg,outlier_df)
  return(returnlist)
}


#**************************************************************************#
#**************************************************************************#
#********** subset array design and data matrix
#**************************************************************************#
#**************************************************************************#

subsetSamples <- function(df, dataset, array_dat, array_design, to_log =T, df_subset =T, model_test =T){
  ## subset data matrix and design file by given factors and levels to keep
  ## INPUT:
  ##    - df(df) contains columns c("eeName", "subset", "subset_by", "keep_subset"), these are required for subsetting
  ##    - array_dat(df): processed expression matrix
  ##    - array_design(df): design file
  ## OUTPUT:
  ##    - new array data, design and log message
  
  ######## sanity checks ########
  if(class(array_dat) != "data.frame"){
    stop(paste0("Input array_dat is a ", class(array_dat), ". A data.frame is required."))
  }
  if(class(df) != "data.frame"){
    stop(paste0("Input to_process is a ", class(df), ". A data.frame is required."))
  }
  if(class(array_design) != "data.frame"){
    stop(paste0("Input array_design is a ", class(array_design), ". A data.frame is required."))
  }
  if(class(dataset) != "character"){
    stop(paste0("Input dataset is a ", class(dataset), ". A string is required, e.g. 'GSE1234'."))
  }
  ## check if df contains the required columns to proceed
  need_col <-c("eeName", "subset", "subset_by", "keep_subset")
  if(!all.equal(need_col, intersect(need_col, colnames(df)))){
    stop(paste0("df has one or more missing columns: ", paste0(need_col, collapse=", ")))
  }
    
    ######## get the subset info #######
  index <- which(df$eeName == dataset)
  if(length(index)!=1){
    stop(paste0(length(index), " record(s) match to the dataset in the to_process file. Only 1 match is allowed"))
  }
  df <- df[index, ] %>% droplevels()
  array_design_old <- array_design
  ######## subset the data if required #######
  # only if yes to subest and content filled in subset_by
  if(tolower(as.character(df$subset)) == 'yes' & as.character(df$subset_by) != ""){
    subset_factor_ls <- unlist(strsplit(as.character(df$subset_by), split ="; "))
    subset_keep_factor_ls <- unlist(strsplit(as.character(df$keep_subset), split ="; "))    
    if(length(subset_factor_ls) != length(subset_keep_factor_ls)){
      stop("Error: Subset_by and keep_subset are not of the same length")
    }
    msg_original <- paste0('Original sample size: ', nrow(array_design))
    
    ## subset the design
    for (i in 1:length(subset_factor_ls)){
      (subset_by <- subset_factor_ls[i])    
      (keep_subset <- subset_keep_factor_ls[i])
      ## subset the design file
      (cmd <- paste0('array_design <- array_design[which(array_design$', subset_by, 
                     ' %in% ', keep_subset, '), ] %>% droplevels()'))
      eval(parse(text = cmd))
    }
    
    ## subset the matrix or reorder the samples
    array_dat <- array_dat[, as.character(array_design$Sample)] %>% droplevels()
    
    subset_rm_df <- NULL
    if(df_subset){
        ## get the removed subseted samples
        rm_index <- which(as.character(array_design_old$ExternalID) %in% as.character(array_design$ExternalID))
        if(length(rm_index) >0 ){
            subset_rm_df <- array_design_old[-rm_index, ]
            subset_rm_df$notes <- 'Removed due to subset'
            subset_rm_df$category <- 'removed by subset'
        }
        
    }

    
    ## test the factor and expression matrix
    response_ls <- intersect(x_col, colnames(array_design))
    response_ls <- response_ls[-length(response_ls)] # remove the last element, Sample
    if(model_test){msg_model <- modelTest(array_dat, array_design, response_ls)}else{msg_model <-'MODEL TEST NOT DONE'}

    ## log
    msg_subset <- paste0('\n#-------------------------------#\n## SUBSET SAMPLES\n#-------------------------------#\n',
                         '\nData is subseted. ', 
                         paste(paste(subset_factor_ls,subset_keep_factor_ls, sep = ": "), collapse=", and "),
                         ' are kept for analysis. \nNew dataset contains ',
                         nrow(array_design), ' samples. ', msg_original)
    
    msg2 <- printSampleCountsByFactor(array_design)
    msg_subset <- paste0( msg_subset, '\n', msg2, msg_model)
    
  }
  
  
  ######### return dataframe ########
  returnlist <- list(array_dat, array_design, msg_subset, subset_rm_df)
  return(returnlist)
}




#**************************************************************************#
#**************************************************************************#
#********** process genmma log history and technical log
#**************************************************************************#
#**************************************************************************#

getGemmaLogInfo <- function(curation_folder="/home/nlim/Curation/Output/"){
  ## return a dataframe combined both History and Technical log files of Gemma uploads
  
  # get the file dir
  f_gemma_history <- paste0(curation_folder, grep("history",list.files(curation_folder), value=T))
  f_gemma_history <- f_gemma_history[length(f_gemma_history)]  # get the last file
  f_gemma_tech <- paste0(curation_folder, grep("Technical",list.files(curation_folder), value=T))
  f_gemma_tech <- f_gemma_tech[length(f_gemma_tech)]  # get the last file
  
  print(paste0("history log: ", f_gemma_history, "\n technical log: ", f_gemma_tech))
  
  # read data
  gemma <- read.delim(f_gemma_history, comment.char = "#")
  gemma_tech <- read.delim(f_gemma_tech, comment.char = "#")
  
  # combine both columns
  gemma <- noWarnings(dplyr::left_join(gemma, gemma_tech))
  gemma$eeName <- as.factor(gemma$eeName)
  
  return(gemma)
}

#**************************************************************************#
#**************************************************************************#
#********** process genmma data design
#**************************************************************************#
#**************************************************************************#

processDesign <- function(dataset, s_file){
  ## read the shorten form, subset the data by dataset, remove NA columns, 
  ## sort rows by sample name, assign sample as row names
  ## return array_design
  ## dataset is the GSE id (string)
  ## design folder is the dir for all gemma design files
  ## s_file = '/home/bzhuang/AD_mouse_model_project/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv'
  if(!file.exists(s_file)){
    stop(paste0("Error: design file", s_file, " is not found."))
  }
  if(class(dataset) != "character"){
    stop(paste0("Input dataset is a ", class(dataset), ". A string is required, e.g. 'GSE1234'."))
  }

  ## read design file
  array_design <- read.delim(s_file,sep="\t", comment.char = "#",stringsAsFactors = T)
  
  ## fix the "FALSE" in the gender (when an exp is only F, it will import as False)
  if("Gender" %in% colnames(array_design)){
    array_design$Gender[which(array_design$Gender ==F)] <- "F"
  }
  
  array_design <- array_design %>% droplevels()
  
  ## subset from the dataset
  #df <- filterContain(array_design,"Dataset", dataset)
  df<- array_design[which(array_design$Dataset ==dataset), ]%>%droplevels()
  
  if(nrow(df) == 0){
    stop(paste0("Error: design file", s_file, " does not contain design info of ", dataset))
  }
  
  
  ## remove NA columns
  df <- df[, colSums(is.na(df)) != nrow(df)]
  ## reorder by sample name
  array_design <- df[with(df, order(Sample)), ]
  row.names(array_design) <- array_design[, "Sample"]
  return(array_design)
}

################################################
################################################
processDesignOriginal <- function(dataset, design_folder){
  ## read design file directly, used for combining multiple design files
  ## dataset is the GSE id (string)
  ## design folder is the dir for all gemma design files (string)
  ## output df is unsorted
  if(!file.exists(design_folder)){
    stop(paste0("Error: Dir for design files ", design_folder, " is not found."))
  }
  if(class(dataset) != "character"){
    stop(paste0("Input dataset is a ", class(dataset), ". A string is required, e.g. 'GSE1234'."))
  }
  
   x <- grep(paste0("_",dataset,"_"), 
             list.files(design_folder), value=T)
   if(length(x) == 0){
     stop(paste0("Error: design file of ", dataset, " is not found."))
   }
   if(length(x) != 1){
     stop(paste0("Error: more than 1 design files found: ", paste0(x, collapse = ", ")))
   }
  (design_file <- paste0(design_folder, "/", x))

  ## read design file
  array_design <- read.table(design_file,sep="\t",
                             header=TRUE,
                             row.names=NULL,fill=T,quote="", check.names = F)
  
  ## if the gender is included (fix the problem when gender is only, which imported as FALSE)
  (temp_gender <- intersect(c('gender', 'Gender'), colnames(array_design)))
  if(length(temp_gender) ==1){
    if(levels(factor(array_design[, temp_gender]))[1] == 'FALSE'){
    array_design[, temp_gender] <- as.factor("F")
    }
  }
  ## simplfy samples names for data design (Sample)
  x <- as.character(array_design$Bioassay)
  new_col <- NULL
  for(i in 1:length(x)){
    new_col <- c(new_col, unlist(strsplit(x[i], split = "Name="))[2])
  }
  array_design$Sample <- new_col
  #sample_order <- as.character(array_design$Sample)
  
  return(array_design)
}

#**************************************************************************#
#**************************************************************************#
#********** process platform
#**************************************************************************#
#**************************************************************************#
processPlatformAnnotation <- function(platform, platform_folder){
  ## read and process platform annotation
  ## return a list of 2 (x), annotation dataframe: as.data.frame(x[1]), process message: as.character(x[2])
  ## platform is id of the microarry platform (GPL...) as string
  ## platfrom_folder is the dir for all annotation files (string)
  if(!dir.exists(platform_folder)){
    stop(paste0("Error: Dir for annotation files ", platform_folder, " is not found."))
  }
  if(class(platform) != "character"){
    stop(paste0("Input platform is a ", class(platform), ". A string is required, e.g. 'GPL1261'."))
  }
  
  x <- grep(paste0(platform, "."), 
            list.files(platform_folder), value=T)

  if(length(x) == 0){
    stop(paste0("Error: annotation file of ", platform, " is not found."))
  }
  if(length(x) != 1){
    stop(paste0("Error: more than 1 platform files found: ", paste0(x, collapse = ", ")))
  }

  anno_file <- paste0(platform_folder, "/", x)
  print(paste0("Platform annotation file", anno_file))
  ## Read platform annotation file into R, gz file or unzipped file
  annotation <- read.table(anno_file,sep="\t",
                           header=TRUE,stringsAsFactors=FALSE,
                           row.names=NULL,fill=T,quote="", check.names = F)
  ## remove NA rows
  temp <- which(is.na(annotation) == TRUE)
  old_rows <- nrow(annotation)
  annotation <- data.frame(na.omit(annotation), check.names=F)   ## omit NA
  msg1 <- paste0('\n#-------------------------------#\n## GENERAL INFO\n#-------------------------------#\n', 
                 "Annotation file is ", anno_file, "\n")
  if(length(temp) > 0){
    rm_rows <- old_rows - nrow(annotation)
    msg1 <- paste0(msg1, "Annotation contain NAs, ", rm_rows, " out of ", old_rows, " removed.")
  }else{
    msg1 <- paste0(msg1, "No probes removed from annotation.")
  }
  print(msg1)
  
  ## remove duplicated probenames, and turn the probenames into rownames (dont column after)
  annotation$ProbeName <- as.character(annotation$ProbeName)  ## must change probe names to string
  annotation <- subset(annotation, !duplicated(annotation[,1]))
  row.names(annotation)<-annotation[,1]
  
  ## return 2 values
  returnlist <- list(annotation, msg1)
  return(returnlist)
}




#**************************************************************************#
#**************************************************************************#
#********** sanity heatmap checks
#**************************************************************************#
#**************************************************************************#
sanityCheckHeatmap <- function(array_dat, f_gene_list, array_design, annotation, 
                               plot_pre, height, width, dataset ="", prefix = "_sanity_check", plt_title = "", 
                               rm_na = T, top_200 = T, scale_row = F, ...){
    ## plot the heatmap for selected genes for sanity check, sample are reordered by genotype (if recorded)
    ## required functions: filterContain()
    ## INPUT:
    ##    - array_dat(df): processed expression matrix
    ##    - f_gene_list(str): a file contains list of gene names in column geneSymbol (e.g. c("Kdm5d", "Xist")), or a list of genes
    ##    - array_design(df): design file
    ##    - annotation(df): platform annotation
    ##    - plot_pre(string): plot dir and prefix, e.g. "/home/folder/prefix"
    ##    - height, width(numeric, in pch): plot dimensions
    ##    - prefix: a string for plot name and saved plot dir
    ##    - rm_na: rm NA of the probe exp, default = T
    ##    - top_200: plot only the first top 200 probes if df is more than 200, default is T
    ##    - ... more options in plotHeatmap(), e.g. breaks, scale etc. 
    
    if(file.exists(f_gene_list)[1]){
        ## load gene list
        gene_list <-read.delim(f_gene_list, comment.char = "#")
        if(!"geneSymbol" %in% colnames(gene_list)){
            stop(paste0("Input gene list ", f_gene_list, " doen't contain column geneSymbol"))  
        }
        gene_ls <- as.character(gene_list$geneSymbol)
    } else {
        print("Input is a list of genes")
        gene_ls <- f_gene_list
    }
    
    ## sanity checks
    if(class(gene_ls) != "character"){
        stop(paste0("Input gene_ls is a ", class(gene_ls), ". A list of strings is required, e.g. c('Xist', 'Kdm5d')."))
    }
    if(class(plot_pre) != "character"){
        stop(paste0("Input plot_pre is a ", class(plot_pre), ". A string is required, e.g. '/home/folder/prefix'."))
    }
    if(class(annotation) != "data.frame"){
        stop(paste0("Input Annotation is a ", class(annotation), ". A data.frame is required."))
    }
    if(class(array_dat) != "data.frame"){
        stop(paste0("Input array_dat is a ", class(array_dat), ". A data.frame is required."))
    }
    if(class(array_design) != "data.frame"){
        stop(paste0("Input array_design is a ", class(array_design), ". A data.frame is required."))
    }
    
    ## replace space in the prefix (remove space for file name)
    prefix <- gsub(" ","_", as.character(prefix))
    
    ## main function
    gene_annotation <- annotation[, c('ProbeName', 'GeneSymbols')]
    df_gene <- filterContain(gene_annotation, column="GeneSymbols", gene_ls)
    # get probe names
    probe_ls <- as.character(df_gene$ProbeName)
    df <- array_dat[probe_ls, ] %>% droplevels()
    if(rm_na) {df <- na.omit(df)}  # remove NA
    # add gene and probe names to row names (for heatmap display)
    gene_names_df <- data.frame(ProbeName = rownames(df))
    gene_names_df <- noWarnings(left_join(gene_names_df, gene_annotation))
    rownames(df) <- paste(gene_names_df$GeneSymbols, gene_names_df$ProbeName, sep = ':')
    
    ## reorder the array_dat samples by Genotype
    if ('Genotype' %in% colnames(array_design)){
        (sample_order <- as.character(array_design[with(array_design, order(Genotype)), "Sample"]))
        df<-df[, sample_order]
    }
    
    ## plot top 200 if top_200 = T
    if(top_200 & nrow(df)>200){
        df_heatmap <- df[1:200, ] %>% droplevels()
        heatmap_title <- paste0(dataset,plt_title, ' top 200 genes ', prefix)
        f_heatmap <- paste0(plot_pre, plt_title, prefix, '_top200.png')
    }else{
        df_heatmap <- df
        heatmap_title <- paste0(dataset,plt_title, prefix)
        f_heatmap <- paste0(plot_pre, plt_title, prefix, '.png')
    }
    
    ## plot heatmap for probes of the genes
    (heat_height <- nrow(df_heatmap)*15)
    (heat_width <- ncol(df_heatmap)*35)
    if(heat_height < height){heat_height <- height}
    if(heat_width < width){heat_width <- width}
    
    png(filename= f_heatmap, width = heat_width, height = heat_height)
    plotHeatmap(df_heatmap, title = heatmap_title, size_r=12, size_c =15,
                cluster_cols = F,
                annotation_row = F,
                array_design_df = array_design, 
                scale_row=scale_row, ...)
    dev.off()
}


#**************************************************************************#
#**************************************************************************#
#********** heatmap for top genes
#**************************************************************************#
#**************************************************************************#
topGeneHeatmap <- function(array_dat, f_gene_list, array_design, annotation, 
                           plot_pre, height, width, dataset ="", prefix = "", plt_title = "", 
                           rm_na = F, top_200 = F, scale_row = F, write_df =T, best_probe =T, auto_w_d = T,
                           size_r=12, size_c =15,array_eb_fit =NULL,display_full_name=T,
                           heat_height = 15, heat_width =35, return_df =F, legend = NULL,
                           cluster_rows =F, cluster_cols =F,verbose =F,
                           ...){
    ## plot the heatmap for selected genes for meta and jack top genes, sample are reordered by genotype (if recorded)
    ## required functions: filterContain()
    ## INPUT:
    ##    - array_dat(df): processed expression matrix
    ##    - f_gene_list(str): a file contains list of gene names in column geneSymbol (e.g. c("Kdm5d", "Xist")), or a list of genes
    ##    - array_design(df): design file
    ##    - annotation(df): platform annotation
    ##    - plot_pre(string): plot dir and prefix, e.g. "/home/folder/prefix"
    ##    - height, width(numeric, in pch): plot dimensions
    ##    - prefix: a string for plot name and saved plot dir
    ##    - rm_na: rm NA of the probe exp, default = F
    ##    - top_200: plot only the first top 200 probes if df is more than 200, default is F
    ##    - write_df: output the expression data with pvalue and logFC
    ##    - scale_row: scale the matrix by row (probe or gene), apply clipping [-2, 2] 
    ##                (ie  if >2, all -to 2, and <-2 all to -2 )
    ##    - best_probe: bool, whether to choose the probe with the best p value to represent the gene; if F, all 
    ##                  probes mapped to the gene will be plotted
    ##    - auto_w_d: auto adjust the width and height of the plot if the required is more the preset width and height
    ##    - array_eb_fit: the limma object loaded from limma results. need this for pvalues if best_probe =T
    ##    - display_full_name: use when best probe =T, if T, rowname = gene and probes, otherwise rowname =gene
    ##                          if best_probe =F, rowname is always gene and probes 
    ##    - heat_height = 15, heat_width =35, the height of a row, and width of a column
    ##    - return_df = F: return the expression df of the genes/probes
    ##    - ... more options in plotHeatmap(), e.g. breaks etc. 
    
    
    ##################
    ## get a list of genes
    ##################
    if(file.exists(f_gene_list)[1]){
        ## load gene list
        gene_list <-read.delim(f_gene_list, comment.char = "#")
        if(!"geneSymbol" %in% colnames(gene_list)){
            stop(paste0("Input gene list ", f_gene_list, " doen't contain column geneSymbol"))  
        }
        gene_ls <- as.character(gene_list$geneSymbol)
    } else {
        #print("Input is a list of genes")
        gene_ls <- f_gene_list
    }
    ##################
    ## sanity checks
    ##################
    if(class(gene_ls) != "character"){
        stop(paste0("Input gene_ls is a ", class(gene_ls), ". A list of strings is required, e.g. c('Xist', 'Kdm5d')."))
    }
    if(class(plot_pre) != "character"){
        stop(paste0("Input plot_pre is a ", class(plot_pre), ". A string is required, e.g. '/home/folder/prefix'."))
    }
    if(class(annotation) != "data.frame"){
        stop(paste0("Input Annotation is a ", class(annotation), ". A data.frame is required."))
    }
    if(class(array_dat) != "data.frame"){
        stop(paste0("Input array_dat is a ", class(array_dat), ". A data.frame is required."))
    }
    if(class(array_design) != "data.frame"){
        stop(paste0("Input array_design is a ", class(array_design), ". A data.frame is required."))
    }
    ##################
    # data should only contain 1 genotype and 1 timepoint before plotting
    # the heatmap
    ##################
    if("Timepoint" %in% colnames(array_design)){
        if(length(levels(array_design$Timepoint)) >1){
            stop(paste0("Input array_design has more than 1 timepoint", 
                        paste0(levels(array_design$Timepoint), collapse = ", ")))
        }
    }
    if("Genotype" %in% colnames(array_design)){
        if(length(levels(array_design$Genotype)) > 2){
            stop(paste0("Input array_design has more than 2 Genotypes: ", 
                        paste0(levels(array_design$Genotype), collapse = ", ")))
        }
    }
    

    
    ## replace space in the prefix (remove space for file name)
    prefix <- gsub(" ","_", as.character(prefix))
    
    ##################
    ## get the expression of genes (all probes)
    ##################
    ## main function
    gene_annotation <- annotation[, c('ProbeName', 'GeneSymbols')]
    df_gene <- noWarnings(left_join(data.frame(GeneSymbols = gene_ls), gene_annotation))
    
    # get probe names
    probe_ls <- as.character(df_gene$ProbeName)
    array_dat_tmp <- array_dat
    array_dat_tmp$ProbeName <- row.names(array_dat)
    df <- noWarnings(left_join(df_gene, array_dat_tmp))
    if(rm_na) {df <- na.omit(df)}  # remove NA
    
    
    ##################
    ## optional: if require to select the best p value probe for the gene
    ##################

    if(best_probe & !is.null(array_eb_fit)){
        library("limma")
        (coef_ls <- dimnames(array_eb_fit$coefficients)[[2]])
        (target_coef <- intersect(coef_ls, paste0("Genotype", levels(array_design$Genotype))))
        if(verbose){print(target_coef)}
        if(length(target_coef) ==1){
            array_toptable <- topTable(array_eb_fit, 
                                       coef = target_coef,
                                       number = Inf) # get all the probes
            array_toptable$ProbeName <- row.names(array_toptable)
            ## get the raw p values for the selected probes
            df <- noWarnings(left_join(df, array_toptable[, c("logFC","P.Value", "ProbeName")]))
            df <- na.omit(df) # must remove NA first otherwise will make multiple matching of NA probes
            
            ## aggregate by min value (if NA will be removed... so need to go back the original gene_ls)
            df2 <-  aggregate(P.Value ~ GeneSymbols, df, function(x) min(x, na.rm = T))
            df2 <- noWarnings(left_join(data.frame(GeneSymbols = gene_ls), df2))
            
            ## get the expression back:
            df3 <- noWarnings(left_join(df2, df))
            if(write_df){
                write.table(df3, file = paste0(plot_pre, plt_title, 'expression_pvalue.tsv'), quote = F, sep ='\t',
                            row.names =F)
            }
            df<-df3
        }
    }
    
    ### for full probe:gene or not for row names
    if(display_full_name){
        rownames(df) <- paste(df$GeneSymbols, df$ProbeName, sep = ':')
    }else{rownames(df) <- df$GeneSymbols}
    
    
    # should be genesymbols(have dup) and probe

    if(best_probe & is.null(array_eb_fit)){
        if(verbose){print('best_probe =T but no array_eb_fit input. Heatmap is represented by all probes mapped to the genes')}
    }
    
    ##################
    ## get ready to plot
    ##################
    # add gene and probe names to row names (for heatmap display)
    
    ## reorder the array_dat samples by Genotype
    if ('Genotype' %in% colnames(array_design)){
        (sample_order <- as.character(array_design[with(array_design, order(Genotype)), "Sample"]))
        df <- df[, sample_order]
    }else{
        df <- df[, setdiff(colnames(df), c('ProbeName', 'GeneSymbols'))]  ## rm character columns, df must be all numbers
    }

    ## plot top 200 if top_200 = T
    if(top_200 & nrow(df)>200){
        df_heatmap <- df[1:200, ] %>% droplevels()
        heatmap_title <- paste0(dataset,plt_title, ' top 200 genes ', prefix)
        f_heatmap <- paste0(plot_pre, plt_title, prefix, '_top200.png')
    }else{
        df_heatmap <- df
        (heatmap_title <- paste0(dataset,plt_title, prefix))
        (f_heatmap <- paste0(plot_pre, plt_title, prefix, '.png'))
    }
    
    df_heatmap_prescale <- df_heatmap
    
    ##################
    ## (optional) zscore scale before plot the heatmap. if data contain NA, heatmap cannot scale this
    ## so first scale; if the range is less than [-2, 2], then the max and min are set to 2 and -2
    ##################
    if(scale_row){
        ## use t to transform by rows!
        if(verbose){print("SCALE by row")}
        df_heatmap <- t(scale(t(as.matrix(df_heatmap))))
        
        ## rescale to -2, 2
        df_heatmap <- scales::rescale(df_heatmap, to=c(-2, 2))
    }
    

    
    ## plot heatmap for probes of the genes
    if(auto_w_d){  # for auto adjust plot weight and height
        (heat_height <- nrow(df_heatmap)*heat_height)
        (heat_width <- ncol(df_heatmap)*heat_width)
        if(heat_height < height){heat_height <- height}
        if(heat_width < width){heat_width <- width}
    }else{
        heat_height <- height
        heat_width <- width
    }

    df_heatmap <- as.data.frame(df_heatmap)
    df_heatmap <- df_heatmap[, setdiff(colnames(df_heatmap), c('ProbeName', 'GeneSymbols','geneSymbol'))]%>%droplevels()  ## rm character columns, df must be all numbers
    
    df_m <- as.matrix(df_heatmap)

    png(filename= f_heatmap, width = heat_width, height = heat_height)
    if(verbose){print('Plotting')}
    plotHeatmap(df_m, title = heatmap_title, 
                size_r=size_r, size_c =size_c,
                cluster_rows = cluster_rows, 
                cluster_cols = cluster_cols,
                annotation_row = F,
                array_design_df = array_design, 
                scale_row=F, 
                legend = legend, ...) # rows are prescaled if requested. don't use the heatmap scale function
    dev.off()
    print(f_heatmap)
    if(return_df){
        return(df_heatmap_prescale)
    }
}






#**************************************************************************#
#**************************************************************************#
#********** match object dir with dataset info
#**************************************************************************#
#**************************************************************************#

# for function with limma toptable dir, see getDatasetLabels()
getRObjectLabels <- function(files, df_info){
    ## get the files (limma_objects.Tdata) in order, get the dataset labels, time labels as specified in df_info
    ## return a dataset with 'Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype', 'File'(object file dir), file_labels (for plotting dataset names)
    ## allow multiple matching (e.g 1 R object match to two genotypes of the same dataset and timepoint)
    
    # get the time point from the file names, and index
    timepoint<- grep("^Timepoint_", unlist(strsplit(as.character(files), split = "/")), value= T)
    timepoint <- gsub("Timepoint_", "", timepoint)
    timepoint <- gsub("months_.*", "months", timepoint)
    
    (timepoint_index<- grep("Timepoint_", files))
    # load the dataset info (contains dataset label, phase label, order of the datasets etc.)
    dataset_info <- read.delim(df_info, comment.char="#")
    
    # check if all the columns are there
    if(! all(c('Dataset', 'Timepoint', 'Phase', 'Order') %in% colnames(dataset_info))){
        stop(paste0("'Dataset', 'Timepoint', 'Phase', 'Order' must be in ", df_info))
    }
    
    # dataset with timepoint specified in the file name
    tmp_df <- data.frame(Dataset = getGSEID(files[timepoint_index]))
    tmp_df$File <- files[timepoint_index]
    tmp_df$Timepoint <- timepoint
    tmp_df <- noWarnings(left_join(tmp_df, dataset_info))
    
    ## datast without timepoint specified in the file name
    tmp_df2 <- data.frame(Dataset = getGSEID(files[-timepoint_index]))
    tmp_df2$File <- files[-timepoint_index]
    tmp_df2 <- noWarnings(left_join(tmp_df2, dataset_info))
    
    ## get all files and sort by order
    tmp_df <- rbind(tmp_df, tmp_df2)
    tmp_df$file_labels <- paste0(tmp_df$Dataset,"_", gsub("onths|eeks|ays", "", tmp_df$Timepoint), "_", tmp_df$Genotype)
        
    tmp_df <- tmp_df[order(tmp_df$Order), ]
    
    return(tmp_df)
}


#**************************************************************************#
#**************************************************************************#
#**********  limma toptable file names, extract the dataset, timepoint,file(file dir),  genotype
#**************************************************************************#
#**************************************************************************#
getTimeFileName <- function(x){
    # x is the full file name (e.g. toptable, R object, up_down gene), return the time label in the file name, return "" if no timelabel
    (tmp <- grep('^[0-9].*_days$|^[0-9].*_months$|^[0-9].*_weeks$|early|late', 
                 unlist(strsplit(gsub("GSE.[0-9,.]+_|_limma|Genotype|-vs-baseline_WT_|.tsv|_gene|_probeset|\\.e_|\\.l_", "#", x), split = "#")), 
                 value = T))
    if(length(tmp) ==0){
        return("")
    }else{
        return(tmp[1])
    }
}


getLimmaLabels <- function(files){
    # input the filenames of limma toptables and extract the dataset, timepoint,file(file dir),  genotype
    # and return a df
    df <- data.frame(Dataset = getGSEID(files))
    df$File <- files
    df$Timepoint <- sapply(files,  function(x) getTimeFileName(x))
    df$Genotype <- sapply(files,  function(x) unlist(strsplit(gsub("Genotype|-vs-baseline_WT_|.tsv", "#", x), split = "#"))[2])
    return(df)
}





#**************************************************************************#
#**************************************************************************#
#********** match limma toptable dir with dataset info
#**************************************************************************#
#**************************************************************************#

# for function with limma object dir, see getRObjectLabels()
getDatasetLabels <- function(files, df_info, match_all = T){
    ## get the files in order, get the dataset labels, time labels as specified in df_info
    ## files can be dir of toptables (with Genotype info)
    ## return a dataset with 'Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype', 'File'(toptable file dir), file_labels (for plotting dataset names)
    #' @param match_all, what' in limma must be in the dataset info too, if F, only select the ones in dataset info
    # limma toptable file names, extract the Dataset, Timepoint,File (file dir),  genotype
    df_limmaname <- getLimmaLabels(files)
    
    # load the dataset info (contains dataset label, phase label, order of the datasets etc.)
    dataset_info <- read.delim(df_info, comment.char="#")
    
    # check if all the columns are there
    required_col <-c('Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype')
    if(! all(required_col %in% colnames(dataset_info))){
        stop(paste0(paste0(required_col, collapse = ', '), " must be in ", df_info))
    }
    # check if df_limmaname and df_info has the same rows
#     if(nrow(df_limmaname) != nrow(dataset_info)){
#         stop(paste0("Toptables has ", nrow(df_limmaname), " rows. Doesn't match with dataset info of ", nrow(dataset_info), ' rows'))
#     }
    
    info_df <- noWarnings(left_join(df_limmaname, dataset_info))
    
    if(match_all){
        # if force to match all limma results and predefined dataset info
        # stop if info_df has NA in Order or Phase
        if(any(is.na(info_df[, c('Order', 'Genotype')]) ==T)){
            print(subset(info_df,is.na(Order)))
            stop(paste0('missing dataset info from dataset_info.tsv'))
        }
    }else{ ## only get the dataset info specified entry
        info_df <- subset(info_df,!is.na(Order))
    }

    
    ## reorder by Order:
    info_df <- info_df[order(info_df$Order), ]
    
    ## add dataset label (file labels), sub months -> m, weeks: w, days: d
    info_df$file_labels <- paste0(info_df$Dataset,"_", gsub("onths|eeks|ays", "", info_df$Timepoint), "_", info_df$Genotype)
    
    return(info_df)
}
