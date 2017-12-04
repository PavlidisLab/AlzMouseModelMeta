#' 2016-08-09
#' this script is to compare the results from mixed model random intercept and random slope
#' get the spearman correlation for all genes and the top genes and plot the rankings
## must run after compare_mm
# see wrapper

(outdir = paste0(slope_f, 'intercept_vs_slope/'))
dir.create(outdir,recursive = T, showWarnings = F)

stat_df = NULL

for(phase in phase_ls){
    for(regulation in regulation_ls){
        f_out_p <- paste0(outdir, phase, '_',regulation,'_intercept', model_keyword, '_slope', model_keyword)
        
        ## get the files
        (i_f <- list.files(intercept_f, full.names = T, pattern = paste0(phase, '_', regulation, '.*.tsv')))
        (s_f <- list.files(slope_f, full.names = T, pattern = paste0(phase, '_', regulation, '.*.tsv')))
        
        ## read intercept
        df <- read.delim(i_f, comment.char = '#')
        rownames(df) <- df$geneSymbol
        
        ## read slope
        df_s <- read.delim(s_f, comment.char = '#')
        rownames(df_s) <- df_s$geneSymbol
        
        ## get the overlap genes
        (index_col<- intersect(colnames(df), c('up_pval','down_pval')))
        genes <- data.frame(geneSymbol = intersect(df$geneSymbol, df_s$geneSymbol))
        genes <- noWarnings(left_join(genes, df[, c('geneSymbol', index_col)]))
        colnames(genes)[which(colnames(genes) == index_col)] <- 'intercept'
        
        genes <- noWarnings(left_join(genes, df_s[, c('geneSymbol', index_col)]))
        colnames(genes)[which(colnames(genes) == index_col)] <- 'slope'
        
        
        genes$rank_intercept <- rank(genes[, 'intercept'], na.last = 'keep')
        genes$rank_slope <- rank(genes[, 'slope'], na.last = 'keep')
        
        
        ## do spearman correlation for all 
        (sp_cor_all <- cor(genes[, 'rank_intercept'], genes[, 'rank_slope'], method = 'spearman'))
        print(paste0(phase, '_',regulation,': intercept', model_keyword, ' vs. slope', model_keyword, " spearman cor: ", sp_cor_all))
        
        ## compare the top genes (after re-rank all common genes)
        (top_meta <- as.character(genes$geneSymbol[order(genes[, 'rank_intercept'])][1:threshold]))
        #print(top_meta)
        (top_MM <- as.character(genes$geneSymbol[order(genes[, 'rank_slope'])][1:threshold]))
        
        df_overlap <- data.frame(geneSymbol = intersect(top_MM, top_meta))
        df_overlap <- noWarnings(left_join(df_overlap, genes))
        df_overlap$avg_rank <- (df_overlap$rank_intercept + df_overlap$rank_slope)/2
        df_overlap <- df_overlap[order(df_overlap$avg_rank), ]
        writeTable(df_overlap, f_out=paste0(f_out_p,  '_overlap_genes.tsv'),
                   msg = paste0('# ', Sys.Date(), ': Top ', threshold, ' overlapped genes',
                                '\n# ', i_f, 
                                '\n# ', s_f))
        
        ## plot and save
        plotXY(genes[, c('rank_intercept', 'rank_slope')], x= 'rank_intercept', y = 'rank_slope',
               title = paste0(phase, ' ', disease, ', ',regulation, 
                              '-regulated (intercept vs. slope)',
                              '\nspearman correlation = ', round(sp_cor_all, 3),', ', nrow(genes), ' genes' ))
        ggsave(filename = paste0(f_out_p, '.png'), width = 8, height=8)
        
        
        ## get the correlation of the top genes
        (outdir_top <- paste0(outdir, '/top_', threshold, '/'))
        dir.create(outdir_top, recursive = T, showWarnings = F)
        genes_top <- genes[which(genes[, 'rank_slope'] <=threshold), ] %>% droplevels()
        
        #top genes: do spearman correlation
        (sp_cor <- cor(genes_top[, 'rank_intercept'], genes_top[, 'rank_slope'], method = 'spearman'))
        print(paste0(phase, '_',regulation,': intercept', model_keyword, ' vs. slope', model_keyword, ": top ", threshold,
                     " spearman cor: ", sp_cor))
        
        
        # top genes: plot and save
        plotXY(genes_top[, c('rank_intercept', 'rank_slope')], x= 'rank_intercept', y = 'rank_slope',
               title = paste0(phase, ' ', disease, ', ',regulation, 
                              '-regulated (intercept vs. slope) ',
                              '\nspearman correlation = ', round(sp_cor, 3),', top ', nrow(genes_top), ' genes' ))
        ggsave(filename = paste0(f_out_p, '_top_', nrow(genes_top), '.png'), width = 8, height=8)
        ## get the stat
        df_t <- data.frame(gene_count = nrow(genes),
                           disease_phase = phase,
                           regulation = regulation,
                           spearman_corr_all_genes = sp_cor_all, 
                           spearman_corr_top_genes = sp_cor, 
                           NA_included = model_keyword,
                           top_genes_overlap = paste0(nrow(df_overlap), '/', threshold),
                           input_intercept = i_f,
                           input_slope = s_f)
        ## summary of all phase and regulation
        stat_df <- rbind(stat_df, df_t)
    }
}






