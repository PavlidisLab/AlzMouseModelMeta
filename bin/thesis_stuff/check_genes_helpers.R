#' helper functions for check_genes.R
#+++++++++++++++++++
# plot pvalue distribution
#+++++++++++++++++++
##############
## thesis plot - histo plot for pvalues before and after cell type correction
##############
library(HelperFunctions)
library(dplyr)
library(ggplot2)


plotPvalue <- function(df, one_plot_font_size =14, f_out){
    library(ggplot2)
    p <- ggplot(df, aes_string(x = 'pvalue')) + 
        geom_histogram(alpha = 0.8, binwidth = 0.005) +
        theme_bw() +
        facet_grid(phase ~ correction,drop = TRUE,scale="free")+
        ggtitle('')
    
    
    p <-   p+xlab('P values') +
        theme(text=element_text(family = 'Arial'),
              #legend.position="none",
              axis.title.x=element_text(size = one_plot_font_size),
              axis.title.y=element_text(size = one_plot_font_size),
              title = element_text(size = one_plot_font_size, colour = 'black'),
              axis.text.x = element_text(size = 10), ## vjust 0.5 put x labels in the middle
              axis.text.y=element_text(size = 10),
              strip.text.x = element_text(size = one_plot_font_size),
              strip.text.y = element_text(size = one_plot_font_size))
    ggsave(filename = f_out, plot=p, width = 6, height =5, units = "in")  
    return(p)
}

#+++++++++++++++++++
# check correlations of ranks before and after 
#+++++++++++++++++++
compareAdjCor <- function(disease, phase, regulation){
    df1 <- all_ranks[which(all_ranks$disease == disease & all_ranks$phase == phase), ]
    df1 <- df1[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df1)[2] <- 'up_jack_mm'
    
    df2 <- all_ranks_adj[which(all_ranks_adj$disease == disease & all_ranks_adj$phase == phase), ]
    df2 <- df2[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df2)[2] <- 'up_jack_mm_adj'
    
    df <- noWarnings(na.omit(left_join(df1, df2)))
    # (correlation <- cor(df[,2], df[,3], method = 'spearman'))
    x <- cor.test(df[,2], df[,3], method = 'spearman', exact = F)
    correlation <- x$estimate
    pvalue <- x$p.value
    df_o <- data.frame(disease=disease, phase = phase, regulation = regulation, correlation = correlation, pvalue = pvalue)
    print(paste(disease,phase, regulation, correlation, collapse = ', '))
    return(df_o)
}



#+++++++++++++++++++
# check correlations of ranks early and late
#+++++++++++++++++++
compareELCor <- function(all_df, disease, prefix, regulation ='up', gene_names=F){
    df1 <- all_df[which(all_df$disease == disease & all_df$phase == 'early'), ]
    df1 <- df1[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df1)[2] <- 'early_mm'
    
    df2 <- all_df[which(all_df$disease == disease & all_df$phase == 'late'), ]
    df2 <- df2[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df2)[2] <- 'late_mm'
    
    df <- noWarnings(na.omit(left_join(df1, df2)))
    x <- cor.test(df[,2], df[,3], method = 'spearman', exact = F)
    correlation <- x$estimate
    pvalue <- x$p.value
    if(gene_names){gene_names <- as.character(paste0(sort(df$geneSymbol), collapse = '; '))}else{
        gene_names =""
    }
    
    df_o <- data.frame(disease=disease,correlation = correlation, prefix = prefix, genes = nrow(df), pvalue=pvalue, gene_names = gene_names)
    print(paste(disease,prefix, correlation, collapse = ', '))
    return(df_o)
}


#+++++++++++++++++++
# check correlations of diseases of early and late
#+++++++++++++++++++
compareADHDCor <- function(all_df, disease, prefix='', phase,regulation ='up'){
    df1 <- all_df[which(all_df$disease == 'AD' & all_df$phase == phase), ]
    df1 <- df1[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df1)[2] <- 'early_mm'
    
    df2 <- all_df[which(all_df$disease == 'HD' & all_df$phase == phase), ]
    df2 <- df2[, c('geneSymbol', paste0(regulation, '_jack'))]
    colnames(df2)[2] <- 'late_mm'
    
    df <- noWarnings(na.omit(left_join(df1, df2)))
    x <- cor.test(df[,2], df[,3], method = 'spearman', exact = F)
    correlation <- x$estimate
    pvalue <- x$p.value
    
    df_o <- data.frame(disease='between AD and HD',correlation = correlation, prefix = paste0(prefix, '_', phase), 
                       genes = nrow(df), pvalue=pvalue)
    print(paste('combined',phase, correlation, collapse = ', '))
    return(df_o)
}


#+++++++++++++++++++
##### check the # of genes up or down regulated within padj<threshold
#+++++++++++++++++++
#' HD late has large # of genes up or down regulated

upDownCount <- function(df, threshold){
    df_count <- NULL
    for(reg in c('up', 'down')){
        if(reg=='up'){
            threshold_index = which(df$pvalue_adj<= threshold & df$Estimate >0)
        }else{
            threshold_index = which(df$pvalue_adj<= threshold & df$Estimate <0)
        }
        
        df_tmp <- data.frame(table(df[,c('disease', 'phase') ]))
        colnames(df_tmp)[3] <- 'total'
        df2 <- noWarnings(left_join(data.frame(table(df[threshold_index,c('disease', 'phase') ])), df_tmp))
        df2$ratio <- df2$Freq/df2$total
        df2$threshold <- threshold
        df2$regulation <- reg
        if(is.null(df_count)){
            df_count= df2
        }else{
            df_count=rbind(df_count, df2)
        }
    }
    return(df_count)
}

#+++++++++++++++++++
##### check the # of genes up or down regulated within padj<threshold
#' write a table and do the thesis bar plot for each disease
#+++++++++++++++++++
##############
## thesis plot - bar plot for DE gene count
##############
plotBarDE <- function(df, one_plot_font_size=10, f_out, x_col = 'regulation',
                      y_col = 'Freq', x_label = 'Regulation', y_label = 'Count',
                      return_p = F, 
                      facet_phase = T,
                      save_p =T){

    
    p <- ggplot(df, aes_string(x = x_col, fill = 'correction')) +
        geom_bar(stat="identity", aes_string(y= y_col), position="dodge") +
        geom_text(aes_string(x=x_col, y=y_col, label='Count'),
                  position = position_dodge(width=0.9), size =3.5) +

        xlab(x_label) +
        ylab(y_label) +
        theme_bw() +
        theme(text=element_text(family = 'Arial'),
              #legend.position="none",
              axis.title.x=element_text(size = one_plot_font_size),
              axis.title.y=element_text(size = one_plot_font_size),
              title = element_text(size = one_plot_font_size, colour = 'black'),
              axis.text.x = element_text(size = one_plot_font_size,colour = 'black'), ## vjust 0.5 put x labels in the middle
              axis.text.y=element_text(size = one_plot_font_size,colour = 'black'),
              strip.text.x = element_text(size = one_plot_font_size),
              strip.text.y = element_text(size = one_plot_font_size))
    
    
    if(facet_phase){
        p <- p + facet_grid(phase ~ .,drop = TRUE,scale="free")
    }
    
    if(save_p){
        print(f_out)
        ggsave(filename = f_out, plot=p, width = 6.45, height =3.15, units = "in")  
    }
    
    if(return_p){
        return(p)}
}


#+++++++++++++++++++
# overlap of genes early vs. late
#+++++++++++++++++++


overlapGenes <- function(df1, disease1, phase1, rank1, df2, disease2, phase2, rank2, threshold){
    tmpFun <- function(df1, disease1, rank1, phase1){
        df1 <- filterContain(df1, 'disease', disease1)
        df1 <- filterContain(df1, 'phase', phase1)%>%droplevels()
        df1 <- orderCol(df1,rank1)
        return(as.character(df1$geneSymbol[1:threshold]))
    }
    
    return(intersect(tmpFun(df1, disease1, rank1, phase1), tmpFun(df2, disease2, rank2,phase2)))
}



#+++++++++++++++++++
## get which DE genes are also cell type markers before and after correction
#+++++++++++++++++++
## thesis tables /results/ND_results/DE_genes_markers/'


getCellMarkers <- function(disease, phase, threshold=NULL, fdr = NULL){
    need_col <- c('disease', 'phase', 'geneSymbol', 'cell_type', 'P_adj', 'Estimate') 
    
    full_name <- c('Oligo', 'Astrocyte', "Microglia","DentateGranule",'GabaSSTReln', 'StriatumCholin','Pyramidal', 'ForebrainCholin', 'Spiny', 'Microglia_deactivation', 'Microglia_activation')
    new_name <- c('Oligodendrocytes','Astrocytes',"Microglia", "Dentate granule cells" ,'GABAergic cells','Cholinergic neurons', 'Pyramidal cells','Cholinergic neurons', 'Medium spiny neurons','Microglia_deactivation', 'Microglia_activation')
    
    
    
    tmpFun <- function(df, prefix, fdr){
        df_marker <- NULL
        for(regulation in c('up', 'down')){

            ## if by fdr
            if(!is.null(fdr)){
                print('count by FDR')
                df1 <- df[which(df$disease == disease & 
                                    df$phase == phase &
                                    df$P_adj < fdr),
                          c(need_col, paste0(regulation,'_jack'))]
                if(regulation =='up'){
                    df1 <- df1[which(df1$Estimate >0), ]
                }else{
                    df1 <- df1[which(df1$Estimate <0), ]
                }
            }else{
                ## if sort by gene ranking
                df1 <- df[which(df$disease == disease & 
                                    df$phase == phase &
                                    df[, paste0(regulation,'_jack')] <= threshold),
                          c(need_col, paste0(regulation,'_jack'))]
            }

            df1[, 'cell_type'] <- plyr::mapvalues(df1[, 'cell_type'], from = full_name, 
                                                             to = new_name, warn_missing = F)
            colnames(df1)[5] <-'Rank'
            df1 <- na.omit(df1)
            df_c <- as.data.frame(table(df1$cell_type))
            
            if(nrow(df_c) == 0){
                df1 <- NULL
            }else{
                colnames(df_c) <- c('cell_type', 'count')
                
                df_tmp <- aggregate(geneSymbol~cell_type, df1, function(x) paste0(sort(x), collapse = ', '))
                df1 <- noWarnings(left_join(df_c,df_tmp))
                df1$regulation = regulation
                
                if(is.null(df_marker)){
                    df_marker <- df1
                }else{
                    df_marker <- rbind(df_marker,df1)
                }
            }
            
        } # end of regulation
        rownames(df_marker) <- paste0(df_marker$cell_type, df_marker$regulation)
        
        colnames(df_marker) <- paste0(colnames(df_marker), prefix)
        return(df_marker)
    } # end of tmpFun
    
    df <- all_ranks
    df_marker <- tmpFun(df, '_before',fdr)
    
    df <- all_ranks_adj
    df_marker_adj <- tmpFun(df,'_after',fdr)
    
    df1 <- mapBindAllColumns(df_marker, df_marker_adj)
    
    df1$disease = disease
    df1$phase = phase
    df1 <- orderCol(df1, c('cell_type_before'))
    #df1$cell_type_before <- df1$cell_type_after[which(is.na(df1$cell_type_before))]
    
    return(df1)
}

enrichedGOTerms <- function(outdir, GO_adj,phase_ls =c('early', 'late')){
    
    ## input the GO_adj rdata, output a clean table with top 20 significant GO terms
    ## and the associated top genes
    
    dir.create(outdir,recursive = T, showWarnings = F)
    df <- GO_adj
    df <- df[, c(1:3,19, 5:9, 11, 15:16)]
    ## get the ranking of the top genes in the gene sets
    
    countTop <- function(x, threshold=100,return_count=T, return_gene=F){
        x <- as.character(x)
        x_char <- unlist(strsplit(x, ', '))
        x <- strsplit(x_char, '\\(|\\)')
        y=NULL
        for(i in 1:length(x)){
            y <- c(y, as.numeric(x[[i]][2]))
        }
        y_len <- length(which(y <= threshold))
        
        if(return_count) return(y_len)
        if(return_gene) return(paste0(x_char[1:y_len],collapse = ', '))
        
    }
    
    
    df$top_count_20 <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =20))
    df$top_count_20_gene <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =20, return_count = F, return_gene = T))
    df$top_20 <- paste0(df$NumGenes, ' (',df$top_count_20, ')')
    
    
    ## remove pathways if there's no top 20 genes involved
    df <- df[which(df$top_count_20 >=2), ]%>%droplevels()
    
    
    ## get the repeated top 20 genes
    overlap_term <- data.frame(top_count_20_gene=df$top_count_20_gene[duplicated(as.character(df$top_genes))])
    overlap_term$overlap_top_20 = 'overlap_20'
    overlap_term <- rmDup(overlap_term)
    
    df_out <- noWarnings(left_join(df,overlap_term))
    
    
    
    ## if top genes are overlapped, maybe select the lowest multifunctionality
    tmp <- df_out[,c('Multifunctionality', 'Name','overlap_top_20', 'top_count_20_gene', 'NumGenes' )]
    tmp$agg <- paste0(tmp$top_count_20_gene, tmp$overlap_top_20)
    
    df2 <- aggregate(Multifunctionality ~ agg, tmp[,c(1:4,6)], function(x) min(x, na.rm = T))
    df2$min_multi <- 'yes'
    tmp <- noWarnings(left_join(tmp, df2))
    
    ## get the pathways with min gene set number
    df3 <- aggregate(NumGenes ~ agg, tmp[,c(1:3,5:6)], function(x) min(x, na.rm = T))
    df3$min_n_genes <- ''
    tmp <- noWarnings(left_join(tmp, df3))
    
    df_out <- noWarnings(left_join(df_out,tmp))
    
    ## final filter:: remove pathways has the same top 20 genes with higher multifunctionality
    df_min <- df_out[which(!is.na(df_out$min_multi)),]
    df_min <- orderCol(df_min, c('disease', 'phase', 'regulation', 'Multifunctionality','min_n_genes'))
    ## replace NA with ''
    index <- which(is.na(df_min$min_n_genes))
    df_min$min_n_genes[index] <- 'No'
    index <- which(is.na(df_min$overlap_top_20))
    df_min$overlap_top_20[index] <- ''
    
    df_min$regulation = paste0(df_min$phase, ': ', df_min$regulation, '-regulated')
    df_min$tmp <- as.factor(paste0(df_min$disease,df_min$regulation))
    
    index <- NULL
    for(keep in levels(df_min$tmp)){
        index <- c(index, which(df_min$tmp == keep)[1])
    }
    
    
    
    ## save tables
    df_min$disease <- as.character(df_min$disease)
    df_min$disease[setdiff(1:nrow(df_min),index)] <- ''
    df_min$regulation <- as.character(df_min$regulation)
    df_min$regulation[setdiff(1:nrow(df_min),index)] <- ''
    
    
    need_col <- c('disease','regulation','phase','ID','Name','P_adj','Multifunctionality','top_20',
                  'top_count_20_gene')
    new_col <- c('Disease','Regulation','phase', 'GO Term ID','GO Term Description','FDR','Multifunctionality Score','Number of Genes',
                 'Top Hits')
    
    
    df_min <- df_min[,need_col]
    colnames(df_min) <- new_col
    
    f <- paste0(outdir, Sys.Date(), "_sig_go_terms.tsv")
    writeTable(df_min, f_out = f)
    print(paste0('File out: ', f))
    return(df_min)
}

