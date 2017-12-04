#' created 2016-06-28
#' Use this after compare_mm_meta.R, with the ordered up or down pvalue list of genes

#' 
#' functions:
#' input a gene name, array_dat, array_design, output model diagnostics, and expression plot

#' @examples 
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease ='AD'
# source('config/config_wrappers.R')
# source('mixed_models/top_genes_mixed_model.R')
# model = 'random_intercept'
# data_dir <- paste0(disease_dir, 'mixed_model/')
# phase <- 'early'
# regulation <- 'up'
# threshold = 5
# 
# mainTopMMGenes(data_dir, phase, regulation, threshold=threshold, 
#                model =model)


#####
 
source('mixed_models/mixed_model.R')

#----------------------------------------------------#
## functions
#----------------------------------------------------#
topGeneExp <- function(f_mm, array_dat, array_dat_not_qn, all_array_dat_bf_aggregation,
                 threshold = 50){
    #' must load expression.Rdata, and mixed_model_results.Rdata
    #' input a up or down value list of genes, set the threshold
    #' output array_dat of the top genes after and before quantile norlization
    #' output probe expressions of the genes from each dataset (a list)
    #' 
    #' 
    df <- read.delim(f_mm, comment.char = '#')
    
    ## sort by up p or down p
    (compare_col <- intersect(colnames(df), c('up_pval', 'down_pval')))
    df <- df[order(df[, compare_col]), ]
    
    ## take the top genes
    df <- df[1:threshold, ]%>%droplevels()
    top_genes <- as.character(df$geneSymbol)
    
    ## top gene expressions
    array_dat <- array_dat[top_genes, ]%>%droplevels()
    array_dat_not_qn <- array_dat_not_qn[top_genes, ]%>%droplevels()
    
    ## get the probe expression of each dataset and output a list
    probe_exp <- vector('list')
    for(i in 1:length(all_array_dat_bf_aggregation)){
        probe_df <- all_array_dat_bf_aggregation[[i]]
        ## get the probes correponding to the genes
        gene_col <- grep('GeneSymbols', colnames(probe_df), value = T)
        probe_df <- probe_df[which(probe_df[, gene_col] %in% top_genes), ]%>%droplevels()
        probe_exp[[i]] <- probe_df
        names(probe_exp)[i] <- names(all_array_dat_bf_aggregation)[i]
    }
    returnlist <- list(array_dat, array_dat_not_qn, probe_exp)
    return(returnlist)
}


plotChange <- function(gene, i, probe_exp, array_design,color_factor=c('Genotype', 'ProbeName')){
    # given a gene in the top threshold genes, and a list of prob expression from topGeneExp()
    # output a plot
    #'@param i: the number of list in the probe_exp, correponding to expression in a dataset
    #'@param color_facotr, color the dots by genotype or probes
    #'
    #'@return: a plot with expresion of wt and disease, colored by genotype
    
    color_factor=color_factor[1]
    
    ## array_design with the matching list name in probe_exp
    array_design$label <- paste0(array_design$Study, '_', gsub('onths|ays|eeks', '', array_design$Timepoint))

    ## subset the design
    array_design_probe <- filterContain(array_design, column = "label", value = names(probe_exp)[i])
    label <- levels(factor(array_design_probe$label))
    
    ## get the expression of the gene
    df <- probe_exp[[i]]
    colnames(df)[grep('GeneSymbols', colnames(df))] <- 'GeneSymbols'
    colnames(df)[grep('ProbeName', colnames(df))] <- 'ProbeName'
    df_exp <- filterContain(df, column = "GeneSymbols", value = gene)
    # make a long table
    df_exp <- melt(df_exp, value.name = 'Expression', variable.name = 'Sample',id.vars = c('ProbeName', 'GeneSymbols'))
    # join the design
    df_exp <- noWarnings(left_join(df_exp, array_design_probe))
    df_exp$Disease_stage <- factor(df_exp$Disease_stage, levels = c('WT', 'Disease'))
    
    ## if no gene, return an empty plot
    if(nrow(df_exp) ==0){
        #stop(paste0(dataset, " doesn't have probes for ", gene))
        df_tmp <- filterContain(df_exp, "Disease_stage", "WT")
        p <- ggplot(df_tmp, aes(Genotype, Expression)) + geom_blank()
    }else{
        p <- ggplot(df_exp, aes_string('Disease_stage', 'Expression', color = color_factor)) + 
            geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size =8)+
            #stat_summary(fun.y=mean, colour="grey", geom="line", aes(group = 1))
          geom_line(stat="summary", fun.y="mean", aes(group=color_factor),alpha = 0.5)
    }
    (title <- paste0(label, "\n", gene, "\n"))
    p <- p + ggtitle(title) +
        theme_bw() +
        theme(legend.position="none", 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title = element_text(size = 32), 
              text =element_text(size = 20), # theme for poster size
              axis.text=element_text(size = 20))
    return(p)
}



multipPlots <- function(gene, probe_exp, array_design,out_f,
                        w_per_data=400,
                        plt_h = 800,
                        color_factor=c('Genotype', 'ProbeName')){
    #' input a gene, 
    #' output the plot for this gene in all input datasets
    #' @param out_f: the output name of the image

    plotlist <- vector('list')
    for(i in 1:length(probe_exp)){
        p <- plotChange(gene, i, probe_exp, array_design, color_factor=color_factor)
        plotlist[[i]] <- p
    }
    
    ## output the plot
    plt_w <- length(plotlist)*w_per_data
    png(filename=out_f, width = plt_w, height = plt_h)
    multiplot(plotlist=plotlist, cols=length(plotlist), layout=NULL, order =F)
    dev.off()
        
    print(paste0("plot: ", out_f))
}


#----------------------------------------------------#
## Main function
#----------------------------------------------------#

mainTopMMGenes <- function(data_dir, phase, regulation, threshold=50, 
                           model ='random_intercept'){
    #' @param data_dir: dir that has r data and the mixed model results(pvalues)
    #' @param
    
    
    ## file dirs
    file_dir <- paste0(data_dir, '/',phase,'/')
    plot_gene_dir <- paste0(file_dir, '/significant_genes/', Sys.Date(), '/',regulation,'/')
    model_dir <- paste0(file_dir, '/model_diagnostics/', regulation,'/')
    dir.create(plot_gene_dir,recursive = T, showWarnings = F)
    dir.create(model_dir,recursive = T, showWarnings = F)
    
    ## get the rdata files and the mixed model results (pvalues)
    r_ls <- list.files(file_dir, recursive = T, full.names = T, pattern = 'expression.Rdata|mixed_model_results.Rdata')
    if(length(r_ls) <2){
        stop(print(paste0('expression.Rdata and mixed_model_results.Rdata not found in ', file_dir)))
    }
    (f_mm <- list.files(data_dir,full.names = T, pattern = paste0(phase,'_',regulation,'.*.tsv')))
    
    if(length(f_mm) !=1){
        stop(paste0(phase,'_',regulation,'.*.tsv not found in ', data_dir))
    }
    
    ## load r data
    for (i in r_ls){
        load(i)
    }
    
    ## get the top genes and correponding expressions
    x <- topGeneExp(f_mm, array_dat, array_dat_not_qn, all_array_dat_bf_aggregation,
                    threshold = threshold)
    array_dat <- x[[1]]
    array_dat_not_qn <- x[[2]]
    probe_exp <- x[[3]]
    
    gene_ls <- as.character(rownames(array_dat))
    # for each gene, plot probe expression, colored by genotype
    for(i in 1:length(gene_ls)){
        
        ## plot  gene expressions
        gene <- gene_ls[i]
        # colored by genotype
        out_f_g <- paste0(plot_gene_dir,'genotype_',i,'_',gene,'_', phase,'_', regulation,'.png')
        multipPlots(gene, probe_exp, array_design,out_f_g,
                    w_per_data=400,
                    plt_h = 800,
                    color_factor='Genotype')
        # colored by genotype
        out_f_p <- paste0(plot_gene_dir,'probename_',i,'_',gene,'_', phase,'_', regulation,'.png')
        multipPlots(gene, probe_exp, array_design,out_f_p,
                    w_per_data=400,
                    plt_h = 800,
                    color_factor='ProbeName')
        
        ## model diagnostic for each top gene
        df_t <- mixedModelStudy(array_dat, array_design, i, model=model, 
                                plot_dir = model_dir,
                                full_report = T,
                                to_plot= T)
        
    }
}