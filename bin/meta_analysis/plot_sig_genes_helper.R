#' 2016-03-02
#' updated 2016-03-03
#' with helper functions to extract expression data and plot functions


library(HelperFunctions)
library(dplyr)
library(ggplot2)
 

#####################################
# helper functions
#####################################
#' given a list of genes and return a df and a list of probes
fromGeneToProbes <- function(gene, annotation){
    df_gene <- filterContain(annotation, column = "GeneSymbols", value = gene)[, 1:2]
    (probe_ls <- as.character(df_gene$ProbeName))
    
    returnlist <- list(df_gene, probe_ls)
    return(returnlist)
}

# given some numbers, return the value with max abs
findMaxAbs <- function(x){
    index <- which(abs(x) == max(abs(na.omit(x))))[1]
    return(x[index])
}

# convert log2FC value to FC
convertLogFC <- function(x){
    y <- 2^(abs(x))
    if(x < 0){
        y <- -1 *y
    }
    return(y)
}

#' extract of all the samples of the specified probes, return NULL if no match
#' Get the pvalue and other from limma expr (for the probes) (same probe for the same genotypes are the same, only non WT
#' has p value appended)
#' @param file_dir: the limma dir
#' 
#'@return a df with "Sample"      "ProbeName"   "Expression"  "Genotype"    "Timepoint"   "Dataset"     "GeneSymbols"
#'@param gene str/str ls. list of gene names
getGeneExp <- function(gene, annotation, array_dat, array_design,file_dir){
    #**********
    # Get the gene exp
    #**********
    x <- fromGeneToProbes(gene, annotation)
    df_gene <- x[[1]]
    probe_ls <- x[[2]]
    if(nrow(df_gene) ==0){
        return(NULL)
    }
    ## get probe exp
    temp <- array_dat[probe_ls, ]
    ## transpose as a tall table
    temp2 <- as.data.frame(t(temp))
    ## add the rownames sample name as the new col 
    temp3 <- temp2 %>%
        mutate(ProbeName = rownames(temp2))
    ## melt the table by ProbeName and update colnames
    temp4 <- reshape2::melt(temp3, id=c("ProbeName"))
    colnames(temp4) <- c("Sample", "ProbeName", "Expression")
    ## dplyr::left_join to add the metadata
    df <- noWarnings(left_join(temp4, array_design[, c("Sample", "Genotype", "Timepoint", "Dataset")]))
    
    ## get the gene symbols
    df <- noWarnings(left_join(df, df_gene))
    
    
    #**********
    # Get the pvalue and other from limma expr (for the probes)
    #**********
    ## grep the top tables
    (top_table_ls <- grep ("toptable_Genotype", list.files(file_dir, recursive=T, full.names=T), value=T))
    print(top_table_ls)
    
    
    ## get all  mouse models from array_design
    (genotype_ls <- setdiff(levels(df$Genotype), "WT"))
    
    df_toptable_all <- NULL
    for(i in 1:length(genotype_ls)){
        (genotype <- genotype_ls[i])
        print(paste0("genotype is ", genotype))
        ## get the top table of that genotype
        (top_table <- grep(genotype, top_table_ls, value=T))
        print(paste0("top_table is ", top_table))
        df_toptable <- read.delim(top_table, comment.char = "#")
        ## only get the probe list
        df_toptable <- filterContain(df_toptable,column="ProbeName", value = probe_ls)
        ## new column for the genotype for matching
        df_toptable$Genotype <- genotype
        df_toptable_all <- rbind(df_toptable_all, df_toptable)
    }
    ## must convert all probenames to char (some probes names are numbers, need to convert to char)
    df_toptable_all$ProbeName <- as.character(df_toptable_all$ProbeName)
    
    ## join the pvalues and the expression value
    df <- noWarnings(left_join(df, df_toptable_all))
    
    return(df)
}

#' one dataset, many genes
#' plot the expression of input gene list df: x, for 1 dataset
#' x is from x <- getGeneExp(gene, annotation, array_dat, array_design,file_dir)
#' 
plotExpOneDataset <- function(x){
    (title <- paste0(levels(x$Dataset), ": ", 
                     paste0(levels(x$Genotype), collapse = " vs."),
                     " (", levels(x$Timepoint), ")"))
    p <- ggplot(x, aes(Genotype, Expression, color = Genotype)) + 
        geom_point(position = position_jitter(width = 0.1))+
        stat_summary(fun.y = mean, geom = "point", shape = 4, size = 5) + 
        stat_summary(fun.y=mean, colour="grey", geom="line", aes(group = 1)) +
        ggtitle(title) +
        facet_wrap( ~ GeneSymbols)
    return(p)
}

#' one dataset, one gene, WT and 1 AD mouse genotype
#' give a gene and a dataset, plot the expr for that gene in the dataset
#' if the dataset doesnt have the gene expr, return an empty plot
#' labels: e.g. "GSE14499_7_months"     "GSE1556_12_months"
#' some datasets has multiple time points
plotExpOneGene <- function(df_all, label, gene){
    df <- filterContain(df_all, column = "label", value = label)
    df <- filterContain(df, column = "GeneSymbols", value = gene)
    
    ## if no gene, return an empty plot
    if(nrow(df) ==0){
        #stop(paste0(dataset, " doesn't have probes for ", gene))
        df_tmp <- filterContain(df_all, "Genotype", "WT")
        (title <- paste0(label, "\n", gene, "\n"))
        p <- ggplot(df_tmp, aes(Genotype, Expression)) + geom_blank()
    }else{
        (qvalue <- round(min(na.omit(unique(df$qValue))),3))
        (pvalue <- round(min(na.omit(unique(df$pValue))),3))
        (logFC <- round(findMaxAbs(df$LogFC),3))
        FC <- round(convertLogFC(logFC),3)
        #(title <- paste0(label, "\n", gene, "\nBest p=",pvalue, "(q=", qvalue, ")", "\n LogFC: ", logFC ))
        (title <- paste0(label, "\n", gene, "\nBest p=",pvalue, "(q=", qvalue, ")", "\n FC: ", FC ))
        p <- ggplot(df, aes(Genotype, Expression, color = ProbeName)) + 
            geom_point(position = position_jitter(width = 0.1), alpha = 0.8, size =8)+
            #stat_summary(fun.y = mean, geom = "point", shape = 4, size = 5) + 
            stat_summary(fun.y=mean, colour="grey", geom="line", aes(group = 1))
    }
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


#' to plot a gene expr of AD and controls for all selected datasets
#' return two plots: 1. the gene for early dataset, 2. the gene for late dataset and save plot
#' @param one_gene: specify a gene (must from a_gene_list, or an empty plot is returned), and plot the exp of AD vs. control of each dataset
#' @param plot_out_dir: dir to save the plot
#' @param return_early_late (bool), return plots for early and late phases separately
#' @param return_all_plt_in_one (bool), return both phases in 1 plot
#' @param APLP2_KO (bool), if genotype APLP2_KO is included (default is included)
plotSigGenes <- function(df_all, one_gene, plot_out_dir, gene_order="",return_early_late = T, 
                         return_all_plt_in_one=F, APLP2_KO =T){
    #*************************************************
    ## plot a specific gene for early samples 1-5m
    #*************************************************
    # "GSE48622_2_months", "GSE36237_5_months", "GSE63617.1_3_months"
    
    i=1
    ## for GSE48622: APP_KO  WT
    label <- "GSE48622_2_months"

    df_tmp <- filterContain(df_all, "Genotype", c("APP_KO","WT"))
    (obj <- paste0("for_multi_plot_early", i))
    assign(obj, eval(parse(text = "plotExpOneGene(df_tmp, label, one_gene)")))
    early_plot_ls <- obj
    
    ## for GSE48622:  APLP2_KOWT
    if(APLP2_KO){
        label <- "GSE48622_2_months"
        df_tmp <- filterContain(df_all, "Genotype", c("APLP2_KO","WT"))
        i=i+1
        (obj <- paste0("for_multi_plot_early", i))
        assign(obj, eval(parse(text = "plotExpOneGene(df_tmp, label, one_gene)")))
        (early_plot_ls <- paste0(early_plot_ls, ", ", obj))
    }

    
    ## for "GSE63617.1_3_months", "GSE36237_5_months"
    label_1_early <- c("GSE63617.1_3_months", "GSE36237_5_months")
    for (label in label_1_early){
        i=i+1
        (obj <- paste0("for_multi_plot_early", i))
        assign(obj, eval(parse(text = "plotExpOneGene(df_all, label, one_gene)")))
        (early_plot_ls <- paste0(early_plot_ls, ", ", obj))
    } 
    
    
    #*************************************************
    ## plot a specific gene for late samples 6-15 m and GSE52022 5XFAD 4m
    #*************************************************
    i=0
    label_1_late <- c("GSE52022_4_months","GSE63617.2_6_months", "GSE14499_7_months", "GSE1556_12_months", "GSE63617.2_15_months","GSE50521_14-15_months")
    plot_ls <- ""
    for (label in label_1_late){
        i=i+1
        (obj <- paste0("for_multi_plot_late", i))
        assign(obj, eval(parse(text = "plotExpOneGene(df_all, label, one_gene)")))
        if(plot_ls==""){
            plot_ls <- obj
        }else{
            plot_ls <- paste0(plot_ls, ", ", obj)
        }
    } 
    
    
    ## asign plots
    assign("plotlist_early", eval(parse(text = paste0("list(", early_plot_ls, ")"))))
    assign("plotlist_late", eval(parse(text = paste0("list(", plot_ls, ")"))))
    
    
    #*************************************************
    ## save plots
    #*************************************************
    
    w_per_data <- 400
    plt_h <- 800
    
    if(return_early_late){
        ## early
        plt_w <- length(plotlist_early)*w_per_data
        out_f <- paste0(plot_out_dir, "early_phases_",gene_order, one_gene, ".png")
        png(filename=out_f, width = plt_w, height = plt_h)
        multiplot(plotlist=plotlist_early, cols=length(plotlist_early), layout=NULL, order =F)
        dev.off()
        
        print(paste0("plot: ", out_f))
        
        ## late
        plt_w <- length(plotlist_late)*w_per_data
        out_f <-paste0(plot_out_dir, "late_phases_",gene_order, one_gene, ".png")
        png(filename=out_f, width = plt_w, height = plt_h)
        multiplot(plotlist=plotlist_late, cols=length(plotlist_late), layout=NULL, order =F)
        dev.off()
        print(paste0("plot: ", out_f))
    }
    
    if(return_all_plt_in_one){
        plt_w <- length(plotlist_early)*w_per_data + length(plotlist_late)*w_per_data
        out_f <- paste0(plot_out_dir, "all_phases_",gene_order, one_gene, ".png")
        png(filename=out_f, width = plt_w, height = plt_h)
        multiplot(plotlist=c(plotlist_early, plotlist_late), cols=length(plotlist_early)+ length(plotlist_late), layout=NULL, order =F)
        dev.off()
        print(paste0("plot: ", out_f))
        
    } 
}




    

#' 
#' give a list of genes and plot the expr of the probes for each of the selected dataset
#' the object dir are hard coded, need to change if updated
#' @example
# setwd("/home/bzhuang/git/bin/mouse_dataset_process")
# plot_out_dir <-("/home/bzhuang/AD_mouse_model_project/results/significant_genes/test/")
# dir.create(plot_out_dir)
# a_gene_list <- c("Trem2", "Urm1")
# rm_multi_gene_probe <- T ## remove probes that map to multiple genes
# return_early_late <- F # not return early and late phase separately 
# return_all_plt_in_one <- T

#' @param a_gene_list: given this list genes and return a df_all with all expression levels from all selected datasets of the genes
#' @param one_gene: specify a gene (must from a_gene_list), and plot the exp of AD vs. control of each dataset
#' @param file_dir_ls: list of dirs where the R object (contain array_dat, array_design, gene annotation etc) of the selected datasets
#' @param with_gene_order (bool): if use a numbered order to name the plots
#' @param rm_multi_gene_probe (bool): remove probes that map to multiple genes

getExpForPlot <- function(file_dir_ls, a_gene_list, rm_multi_gene_probe = T, plot_out_dir ="", 
                          with_gene_order=F, 
                          return_early_late=T,
                          return_all_plt_in_one =F,
                          APLP2_KO =T){
    #**************************
    ## get expression level of the genes and save
    #**************************
    df_all <- NULL
    for (file_dir in file_dir_ls){
        
        all_files <- list.files(file_dir, full.names = T, recursive = T)
        (object_dir <- grep("objects", all_files, value = T))
        print(object_dir)
        load(object_dir)
        print("object loaded")
        df <- getGeneExp(a_gene_list, annotation, array_dat, array_design,file_dir)
        df_all <- rbind(df_all, df)
        print(nrow(df_all))
    }
    
    df_all <- df_all[!duplicated(df_all), ]
    df_all$label <- paste0(df_all$Dataset, "_", df_all$Timepoint)
    df_all$label <- as.factor(df_all$label)
    
    if(rm_multi_gene_probe){
        ## remove multiple genes
        index <- grep('\\|', df_all$GeneSymbols)
        if(length(index) != 0){df_all <- df_all[-index, ] %>%droplevels()}
    }
    
    ## relevel df_all (WT as first factor)
    df_all$Genotype <- as.factor(df_all$Genotype)
    geno_ls <- levels(df_all$Genotype)
    (geno_ls <- c("WT", setdiff(geno_ls, "WT")))
    df_all$Genotype <- factor(df_all$Genotype, levels = geno_ls)
    
    ## update the gene list (rm genes that are not in the output gene symbols)
    a_gene_list <- intersect(a_gene_list, unique(df_all$GeneSymbols))
    
    print("df_all returned")
    write.table(df_all, file = paste0(plot_out_dir, "gene_expression.tsv"), sep = "\t",
                row.names = F,  quote = F)
    
    
    #**************************
    ## plot each gene (loop)
    #**************************
    counter <- 0
    for (one_gene in a_gene_list){
        counter <- counter +1
        if(with_gene_order){
            gene_order <- paste0("_", counter, "_")
        }else{
            gene_order <- ""
        }
        plotSigGenes(df_all, one_gene, plot_out_dir, gene_order=gene_order,return_early_late, return_all_plt_in_one, APLP2_KO)
    }
    
    
    return(df_all)
}


#' @param file_dir_ls: list of dirs where the R object (contain array_dat, array_design, gene annotation etc) of the selected datasets
#' @param with_gene_order (bool): if use a numbered order to name the plots
#' @param rm_multi_gene_probe (bool): remove probes that map to multiple genes
#' @examples
# file_dir_ls <- c("/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE14499/", 
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE36237/", 
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-15/GSE48622/",
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-02-15/GSE1556/", 
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.1/", 
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_6_months_OrganismPart_Hippocampus/",
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-02-23/GSE63617.2/subset/Timepoint_15_months_OrganismPart_Hippocampus/",
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/gemma_prioritized_limma/2016-03-07/GSE50521/",
#                 "/home/bzhuang/AD_mouse_model_project/limma_DE/hippocampus/CEL_prioritized_limma/2016-04-05/GSE52022/subset/Timepoint_4_months/")
# 
# plot_out_dir <- '/home/bzhuang/AD_mouse_model_project/meta_analysis/meta_CEL_and_gemma/test/prioritized_genes/plots_up/'
# f_gene_list <- f_out_up
# plotExpWithGeneList(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel=T, with_gene_order=T)

plotExpWithGeneList <- function(f_gene_list, file_dir_ls, plot_out_dir, return_one_panel, with_gene_order,
                                APLP2_KO =T){
    ## give a table with geneSymbol, plot the expression(from limma DE expression) of each gene
    dir.create(plot_out_dir, recursive= T)
    rm_multi_gene_probe <- T ## remove probes that map to multiple genes
    if (return_one_panel){
        return_early_late <- F # all phases in 1 panel
        return_all_plt_in_one <- T
    }else{
        return_early_late <- T # return early and late phases separately 
        return_all_plt_in_one <- F
    }
    df <- read.delim(f_gene_list, comment.char="#")
    ## plot expression
    (a_gene_list <- as.character(df$geneSymbol))
    getExpForPlot(file_dir_ls, a_gene_list, rm_multi_gene_probe,plot_out_dir = plot_out_dir, 
                  with_gene_order = with_gene_order, 
                  return_early_late,
                  return_all_plt_in_one,
                  APLP2_KO)

}

####################################################
## final wrapper
####################################################
#' use the jackknife as input, plot the expression of the top genes
#' plot expression of cell type specific genes (not related to jackknife)
#' plot expression of gender genes (not related to jackknife)
#' @param jack_dir: folder containing jackknife files
#' @param plot_out: folder for plots
#' @param file_dir_ls: folder for the Rdata for gene expression
#' @param phase_ls, list of phases to be included
#' @param return_one_panel: return all datasets in 1 panel
wrapperPlotGeneExp <- function(jack_dir, plot_out, file_dir_ls, plt_jack =T, 
                               plt_cell =F, plt_gender =F, phase_ls = c("early", "late", "med"),
                               return_one_panel = T,
                               APLP2_KO =T){
    
    #---------------------------------------------------------------------------#
    # PART1A: PLOT GENE EXPRESSION of Jackknife genes
    #---------------------------------------------------------------------------#
    # after meta and jackknife, plot some top genes
    #*************************
    ## plot top 10 or FRD < 0.1 jackknife genes
    #*************************
    if(plt_jack){
        print("PLOTTING SIGNIFICANT JACKKNIFE GENES")
        (search_phase <- paste0("(", paste0(phase_ls, collapse="|"), ").*jackknife"))
        (jack_files <- grep(search_phase, list.files(jack_dir, full.names=T), value=T))
        (jack_label <- gsub("_regulation_jackknife.tsv", "", grep(search_phase, list.files(jack_dir), value = T)))
        dir.create(plot_out, recursive= T)
        plot_out_sig <- paste0(plot_out, "/Jackknife_sig_genes/")
        dir.create(plot_out_sig, recursive= T)
        
        with_gene_order <- T  ## give the plot names with an order
        rm_multi_gene_probe <- T ## remove probes that map to multiple genes
        if (return_one_panel){
            return_early_late <- F # all phases in 1 panel
            return_all_plt_in_one <- T
        }else{
            return_early_late <- T # return early and late phases separately 
            return_all_plt_in_one <- F
        }
        
        ## for loop
        for (x in 1: length(jack_files)){
            print(x)
            print(paste0("processing ", jack_label[x]))
            ## create output dir
            plot_out_dir <- paste0(plot_out_sig, jack_label[x], "/")
            dir.create(plot_out_dir, recursive=T, showWarnings=F)
            
            ## read jack files
            (jack <- jack_files[x])
            df <- read.delim(jack, comment.char="#")
            
            ## get genes under the threshold or the top 10
            index <- which(df$adj_combined_max_p <0.1)
            if(length(index) > 10){
                df <- df[index, ]
            }else{
                df <- df[1:10, ]
            }
            
            # save the file
            write.table(df, file= paste0(plot_out_dir, jack_label[x], "_top_genes.tsv"), sep = "\t", quote =F, row.names = F)
            
            ## plot expression
            (a_gene_list <- as.character(df$geneSymbol))
            getExpForPlot(file_dir_ls, a_gene_list, rm_multi_gene_probe, plot_out_dir, 
                          with_gene_order, 
                          return_early_late,
                          return_all_plt_in_one,
                          APLP2_KO)
        }
    }
    
    
    
    #---------------------------------------------------------------------------#
    # PART2: PLOT CELL POPULATION gene exp
    #---------------------------------------------------------------------------#
    if(plt_cell){
        print("PLOTTING CELL TYPE SPECIFIC GENES")
        
        # use Ogan's gene list of pyramidalDEEP/hippocampus genes
        plot_out_cell <-paste0(plot_out, "/Cell_type_gene_exp/")
        dir.create(plot_out_cell, showWarnings=F, recursive=T)
        
        cell_pop_dir <-"./brainCellTypeSpecificGenes-master/analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Hippocampus/"
        (cell_pop_files <- list.files(cell_pop_dir, full.names=T))
        (cell_pop_label <- list.files(cell_pop_dir))
        
        with_gene_order <- F  ## give the plot names with an order
        rm_multi_gene_probe <- F ## do not remove probes that map to multiple genes (some of the cell type specific markers are probes mapped to multiple genes)
        return_early_late <- F
        return_all_plt_in_one <- T # return all phases in 1 plot
        
        ## for loop
        for (x in 1: length(cell_pop_files)){
            print(x)
            print(paste0("processing ", cell_pop_label[x]))
            ## create output dir
            (plot_out_dir <- paste0(plot_out_cell, cell_pop_label[x], "/"))
            dir.create(plot_out_dir, recursive=T, showWarnings=F)
            
            ## read cell_pop files
            cell_pop <- cell_pop_files[x]
            df <- read.delim(cell_pop, comment.char="#", header =F)
            
            # save the file
            write.table(df, file= paste0(plot_out_dir, cell_pop_label[x], "_cell_type_specific_genes.tsv"), sep = "\t", quote =F, row.names = F)
            
            ## plot expression
            (a_gene_list <- as.character(df[,1]))
            getExpForPlot(file_dir_ls, a_gene_list, rm_multi_gene_probe, plot_out_dir, 
                          with_gene_order, 
                          return_early_late,
                          return_all_plt_in_one,
                          APLP2_KO)
            
            # save the file
            f_input_gene_list <- paste0(plot_out_dir, cell_pop_label[x], "_input_gene_list.tsv")
            write.table(a_gene_list, file= f_input_gene_list, sep = "\t", quote =F, row.names = F, col.names=F)
            print(paste0("input gene list: ", f_input_gene_list))
            
        }
        
    }
    
    #---------------------------------------------------------------------------#
    # PART3: plot sex genes
    #---------------------------------------------------------------------------#
    if(plt_gender){
        print("PLOTTING GENDER GENES")
        plot_out_gender <-paste0(plot_out, "/gender_genes/")
        dir.create(plot_out_gender,showWarnings=F)
        a_gene_list <- c("Xist", "Kdm5d")
        with_gene_order <- F  ## give the plot names with an order
        rm_multi_gene_probe <- T ## remove probes that map to multiple genes
        return_early_late <- F # not return early and late phase separately 
        return_all_plt_in_one <- T
        getExpForPlot(file_dir_ls, a_gene_list, rm_multi_gene_probe, plot_out_gender, 
                      with_gene_order, 
                      return_early_late,
                      return_all_plt_in_one)
    }
    
}
