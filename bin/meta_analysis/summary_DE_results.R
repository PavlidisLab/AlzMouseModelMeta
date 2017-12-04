# updated 2016-03-21
# use 'wrapper.R', "low_exp_rm_analysis/wrapper_filter_low.R"
source("helper_functions.R")
library('HelperFunctions')

# 1. summary table for each DE, how many FDR <0.05, FDR<0.1
# - get the plots of DE (p value distributions)
# the wraper scripts are the end of the script

#------------------------------------------------#
# helper functions
#------------------------------------------------#

# input pvalues and output DE ratio and esitmated sd of DE ratio
#' @param df: the toptable df with pValue
DERatioEstimate <- function(df){
    library(qvalue) ## estimate pi0 and DE ratio: 1-pi0
    library(boot)  ## bootstrap
    #' pi0: Estimates the proportion of true null p-values, i.e., those following the Uniform(0,1) distribution
    #' The bootstrap method is described in Storey, Taylor & Siegmund
    #' (2004). A closed form solution of the bootstrap method is used in the package and is significantly
    #' faster.
    #' DE-ratio = 1- pi0
    #' repeat the process by bootstrap 200 time to estimate sd of the DE ratio
    
    # function to obtain regression weights
    DERatio <- function(pvalue, indices) {
        d <- pvalue[indices] # allows boot to select sample
        pi0_obj <- pi0est(d, pi0.method = 'bootstrap')
        (DE_ratio <- 1 - as.numeric(pi0_obj[1]))
        return(DE_ratio)
    }
    # bootstrapping with 200 replications ORDINARY NONPARAMETRIC BOOTSTRAP
    #' to get the standard errors
    pi0_results <- boot(data=df$pValue, statistic=DERatio,
                        R=200)
    DE_ratio <- pi0_results$t0
    DE_ratio_sd <- sd(pi0_results$t)
    returnlist <- list(DE_ratio, DE_ratio_sd)
    return(returnlist)
}



#------------------------------------------------#
# output functions
#------------------------------------------------#
#------------------------------------------------#
## a independent function (not included in the main function)
#------------------------------------------------#
getSampleSummary <- function(datadir, f_out,f_DE_summary = NULL,plot_DE=T){
    # input a datadir with limma results. all the limma_table_log.tsv will be searched (recursively) and combine to a big table
    # if f_DE_summary (summary of all DE) is available the add the new info
    # and output a table with all limma results
    # f_out: the outdir for output table and a DE ratio vs sample size plot
    (file_ls <- list.files(datadir, recursive=T, full.names=T, pattern = "limma_table_log.tsv"))
    f_tb_out <- paste0(f_out, Sys.Date(), '_input_data_summary.tsv')
    f_plot_out <- paste0(f_out, Sys.Date(), '_DE_ratio_and_sample_size.png')
    
    df_all <- NULL
    for (i in file_ls){
        df <- read.delim(i, comment.char ="#")
        if(is.null(df_all)){
            df_all <- df
        }else{
            df_all <- rbindAllColumns(df_all, df)
        }
    }
    df_all$Timepoint <- as.factor(gsub("\\([1-9].*", "", df_all$Timepoint)) ## remove the sample size after timepoint
    
    if(!is.null(f_DE_summary)){ # add the summary data
        df_summary <- read.delim(f_DE_summary, comment.char='#')
        df_all <- noWarnings(left_join(df_all, df_summary))
    }
    msg <- paste0('saved file: ', f_tb_out)
    
    ## plot sample size and DE ratio (optional)
    if(plot_DE){
        need_col <- c("sample_size_all","DE_ratio",'phase','DE_ratio_sd')
        if(all(need_col%in% colnames(df_all))){# must contain these columns
            print('PLOT DE ratio vs sample size')
            df_plot <- na.omit(df_all[, need_col]%>%droplevels())
            
            ## test correlation of sample size and DE ratio
            cor_test <- cor.test(df_plot$sample_size_all, df_plot$DE_ratio)
            cor_estimate <- as.numeric(cor_test["estimate"])
            cor_p <- as.numeric(cor_test["p.value"])
            df_all$DE_sample_size_cor <- cor_estimate
            df_all$DE_sample_size_cor_p <- cor_p

            ## add the sd value
            df_plot$ymin=df_plot$DE_ratio - df_plot$DE_ratio_sd
            df_plot$ymax=df_plot$DE_ratio + df_plot$DE_ratio_sd
            
            ## plot and save
            ggplot(df_plot, aes_string("sample_size_all","DE_ratio",group = 'phase', colour='phase')) +
                geom_point() +
                theme_bw() +
                xlab(paste0("Sample Size"))+ 
                ylab(paste0("DE ratio"))+
                geom_errorbar(aes_string(ymin= 'ymin', ymax = 'ymax'))+
                ggtitle(paste0("Sample size VS. DE ratio\n Pearson correlation: ", 
                               round(cor_estimate, 3),
                               ', p-value = ', round(cor_p, 3))) 
            ggsave(filename = f_plot_out, width = 8, height=5)
            
            msg <- paste0(msg, '\nsaved plot: ', f_plot_out)
        }
    }
    
    cat(msg)
    
    ## write file
    msg = paste0('# ', Sys.Date(), 
                 '\n# Input summary file: ', f_DE_summary)
    writeTable(df_all, f_out= f_tb_out, msg=msg)
}




#------------------------------------------------#
## all the counts are by probe/ probesets
#------------------------------------------------#
summaryToptable <- function (info_df, i){
    ## info_df: a df with 'Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype', 'File'(toptable file dir), file_labels (for plotting dataset names)
    ##          (outcome of getDatasetLabels())
    ##  i: the row index of info_df
  # files: list of toptable files
  # file_labels, label for the x lab e.g. GSE1234_4_m_APP_KO
  # time_labels, label for disease phase, e.g. early, late
  # the orders must match
    print(i)
    print(info_df$File[i])
    # i is the index
    (one_f <- info_df$File[i])
    (one_l <- info_df$file_labels[i])
    (one_t <- info_df$Phase[i])
    (one_timepoint <- info_df$Timepoint[i])
    
    ## read the toptables
    df <- read.delim(one_f, comment.char="#")
    (dataset <- getGSEID(one_f)[1])
    (filename <- one_f)
    
    ## estimate DE ratio: 1- pi0
    tmp <- DERatioEstimate(df)
    DE_ratio <- tmp[[1]]
    DE_ratio_sd <- tmp[[2]]
    
    ## summary of genes
    (index <- which(df$qValue <0.05))
    FDR5 <- length(index)
    FDR5genes <- paste0(setdiff(unique(df$GeneSymbols[index]), ""), collapse= "; ") # rm empty genes
    FDR10 <- length(which(df$qValue <0.1))
    pvalue5 <- length(which(df$pValue <0.05))
    (p_adj_up5 <- length(which(df$up_padj < 0.05)))
    (p_adj_down5 <- length(which(df$down_padj < 0.05)))
    (p_adj_up10 <- length(which(df$up_padj < 0.1)))
    (p_adj_down10 <- length(which(df$down_padj < 0.1)))
    
    ## get the genes (rm multi genes)
    df_gene <- excludeMatch(df, "GeneSymbols", "")
    index <- grep("\\|", df_gene$GeneSymbols)
    df_gene <- df_gene[-index, ]%>% droplevels()
    
    countGene <- function(i){
        length(unique(df_gene$GeneSymbols[i]))
    }
    
    
    (index <- which(df_gene$qValue <0.05))
    (FDR5_by_gene <- countGene(index))
    
    FDR10_by_gene <- countGene(which(df_gene$qValue <0.1))
    
    (p_adj_up5_by_gene <-  countGene(which(df_gene$up_padj < 0.05)))
    (p_adj_down5_by_gene <-  countGene((which(df_gene$down_padj < 0.05))))
    (p_adj_up10_by_gene <-  countGene(which(df_gene$up_padj < 0.1)))
    (p_adj_down10_by_gene <- countGene(which(df_gene$down_padj < 0.1)))
    
    total_probes <- nrow(df)
    total_genes <- length(unique(df_gene$GeneSymbols))
    
    ####
    ### make a df with the sig genes (for bar plot)
    ###
    ## get the sig up or down genes and return a df
    keyword <- "up_padj"
    
    if(p_adj_up5_by_gene == 0){
        df_tmp <- NULL
    }else{
        (df_tmp <- df_gene[which(df_gene[, keyword] < 0.05), c("GeneSymbols", keyword)] %>% droplevels)
        colnames(df_tmp) <- c("GeneSymbols", "adj_p")
        df_tmp$Dataset <- one_l
        df_tmp$Phase <- one_t
        df_tmp$Regulation <- "Up-regulated"
    }
    
    keyword <- "down_padj"
    if(p_adj_down5_by_gene ==0){
        df_tmp2 <- NULL
    }else{
        (df_tmp2 <- df_gene[which(df_gene[, keyword] < 0.05), c("GeneSymbols", keyword)] %>% droplevels)
        colnames(df_tmp2) <- c("GeneSymbols", "adj_p")
        df_tmp2$Dataset <- one_l
        df_tmp2$Phase <- one_t
        df_tmp2$Regulation <- "Down-regulated"
    }

    df_gene_bar <- rbind(df_tmp, df_tmp2)
  
    
    ## the summary table
    df_summary <- data.frame(dataset= dataset, 
                             dataset_label =one_l,
                             Timepoint = one_timepoint,
                             Genotype = info_df$Genotype[i],
                             phase = one_t,
                             DE_ratio = DE_ratio,
                             DE_ratio_sd = DE_ratio_sd,
                             file = filename, 
                             total_probes = total_probes, 
                             fdr0.05_probe = FDR5,  
                             pvalue0.05_probe = pvalue5, 
                             p_adj_up0.05_probe = p_adj_up5,
                             p_adj_down0.05_probe = p_adj_down5,
                             fdr0.1_probe = FDR10,
                             p_adj_up0.1_probe = p_adj_up10,
                             p_adj_down0.1_probe = p_adj_down10,
                             
                             total_genes = total_genes, 
                             fdr0.05_gene = FDR5_by_gene,  
                             p_adj_up0.05_gene = p_adj_up5_by_gene,
                             p_adj_down0.05_gene = p_adj_down5_by_gene,
                             fdr0.1_gene = FDR10_by_gene,
                             p_adj_up0.1_gene = p_adj_up10_by_gene,
                             p_adj_down0.1_gene = p_adj_down10_by_gene,
                             fdr0.05_genes = FDR5genes)
    #print(df_gene_bar)
    returnlist <- list(df_summary, df_gene_bar)
    return(returnlist)
}

## plot bar graph:
plotBarPoster <- function(df, key_col, plot_title, fo_s, poster_style =T, bar_text = 8, facet_phase =F){
    #' @param facet_phase, facet grid by disease_phase
    p <- ggplot(df, aes(x = dataset, fill = regulation)) +
        geom_bar(stat="identity", aes_string(y=key_col), position="dodge") +
        geom_text(aes_string(x="dataset", y=key_col,label=key_col), 
                  position = position_dodge(width=0.9), size= bar_text) +  # asign count text to each bar
        ggtitle(plot_title) +
        theme_bw() +   ## white background
        ylab("Count") +
        xlab("") +
        theme(plot.title = element_text(size = 32), 
              text =element_text(size = 20), # theme for poster size
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20),
              axis.ticks = element_blank(), axis.text.x = element_blank())
    if(!poster_style){p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))} ## add sample names
    if(facet_phase){p + facet_grid(disease_phase ~.)} ## 2 panels by phase, no sample names
p + ggsave(filename = fo_s, width = 14, height =6)




}




#------------------------------------------------#
# main function
#------------------------------------------------#

mainSummaryDE <- function(files, f_out, df_info, df_return=F, known_modifier ="", bar_text=7, match_all =T){
  #' files: list of files of toptables (full name), dont' have to be in order
  #' df_info: a table specify the order of the datasets, phase (early, late etc), timepoint, dataset label for plot etc
  #' known_modifier is the list obtained from Joerg, with geneSymbol
  #' OUTPUT:
  #' 2 bar plots (1 with x labels and the other not)
  #' a sample size vs DE ratio plot
    
    if(length(files) <1){
        stop('no toptable files listed')
    }
    dir.create(f_out,showWarnings=F, recursive=T)
    
    ## table and plot out dirs
    out_plt_p <- paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts_bar_plot_poster.png") # no label on x
    out_plt_l <- paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts_bar_plot.png") # with label on x
    out_plt_fc <- paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts_bar_plot_phase.png") # early and late panels
    f_out_count <-paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts.tsv")
    f_out_sig_gene <-paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_0.05_genes.tsv")
    
    ## get the files and labels in order
    info_df <- getDatasetLabels(files, df_info, match_all=match_all)
    file_labels <- info_df$file_labels
    time_labels <- info_df$Phase
    files <- info_df$File
      
    ## get the counts
    df_all <- NULL
    df_gene_bar_all <- NULL 
    for (i in 1:length(files)){
      #print(files[i])
        x <- summaryToptable(info_df, i)
        
        df <-x[[1]] 
        df_all <- rbind(df_all, df)
        
        df_gene_bar <- x[[2]]
        df_gene_bar_all <- rbind(df_gene_bar_all,df_gene_bar)
    }
    
    df_gene_bar_all$Dataset <- factor(df_gene_bar_all$Dataset, levels = file_labels)
    
    ## log
    msg_count <- paste0("# ", Sys.Date(), " Summary of each DEA")
    
    if(known_modifier != ""){
        df_known <- read.delim(known_modifier, comment.char="#")
        df_gene_bar_all <- noWarnings(left_join(df_gene_bar_all, df_known[, c("Gene_name", "GeneSymbols", "Type_of_modification", "Pubmed_ID")]))
        df <- na.omit(df_gene_bar_all)
        df$gene_modifier <- paste0(df$GeneSymbols, "(", df$Gene_name, ")")
        (msg_count <- paste0(msg_count, "\n# known modifier genes: ", length(unique(df_known$Gene_name)), 
                            "(", length(unique(df_known$GeneSymbols)) , " genesymbols)",
                       "; Found in DE: ", length(unique(df$gene_modifier)),
                       "\n# ", paste0(df$GeneSymbols, "(", df$Gene_name, "): ",  df$Dataset, collapse=", ")))
        
        df2 <- df[, c("gene_modifier", "Dataset")]
        df2 <- data.frame(aggregate(gene_modifier ~ Dataset, data=df2, function(x) paste0(x, collapse= ", ")))
        colnames(df2)[1] <- "dataset_label"
        df_all <- noWarnings(left_join(df_all, df2))

    }
    

    
    ## save results 
    # save counts
    sink(f_out_count, type="output")
    writeLines(msg_count)
    sink()
    noWarnings(write.table(df_all, file = f_out_count, row.names = F,
                sep ='\t', quote = F, append=T))
    #save gene names
    sink(f_out_sig_gene, type="output")
    writeLines(msg_count)
    sink()
    noWarnings(write.table(df_gene_bar_all, file = f_out_sig_gene, row.names = F,
                sep ='\t', quote = F, append=T))
    
    
    ## prep df for plots
    df_plot_up <- df_all[, c("dataset_label", "p_adj_up0.05_gene", "phase")]%>% droplevels()
    colnames(df_plot_up) <- c("dataset", "significant_genes", "disease_phase")
    df_plot_up$regulation <- "up"
    
    df_plot_d <- df_all[, c("dataset_label", "p_adj_down0.05_gene", "phase")]%>% droplevels()
    colnames(df_plot_d) <- c("dataset", "significant_genes", "disease_phase")
    df_plot_d$regulation <- "down"
    
    df_plot <- rbind(df_plot_up, df_plot_d)
    df_plot$dataset <- factor(df_plot$dataset, levels = file_labels)
    df_plot$regulation <- factor(df_plot$regulation, levels = c("up", "down"))
    
    ## plots
    key_col <- "significant_genes"
    plot_title <- "Significant genes (FDR <0.05)"
    plotBarPoster(df_plot, key_col, plot_title, out_plt_p, bar_text= bar_text)
    plotBarPoster(df_plot, key_col, plot_title, out_plt_l, poster_style =F, bar_text= bar_text)
    plotBarPoster(df_plot, key_col, plot_title, out_plt_fc, poster_style =T, bar_text= bar_text, facet_phase = T)
    
    print(paste0("Results are in ", f_out))
    
    if('disease_phase' %in% colnames(df_plot)){
        for(phase in as.character(levels(df_plot$disease_phase))){
            plot_title <- paste0(toupper(phase), " phase: significant genes (FDR <0.05)")
            df_tmp <- filterContain(df_plot, column = 'disease_phase', value = phase)
            plt_p =paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts_bar_plot_poster_", phase, ".png") # no label on x
            plt_l =plt_p =paste0(f_out, "/", Sys.Date(), "_limma_DE_FDR_counts_bar_plot_", phase, ".png")
            plotBarPoster(df_tmp, key_col, plot_title, plt_p, bar_text= bar_text)
            plotBarPoster(df_tmp, key_col, plot_title, plt_l, poster_style =F, bar_text= bar_text)
        }
    }
    
    
    
    
#     returnlist <- list(df_plot, key_col, plot_title, out_plt_p)
#     return(returnlist)
    
    if(df_return){
        return(df_all)
    }
    
}

