# 2016-08-04
#' this script plot a big figure with all samples from early or late phase
#' with the input expression from mixed model quantile normalized
#' the heatmap includes top up and down genes
# topGeneHeatmap from helper_function.R
# plot name is paste0(plt_out, plt_title, '.png')



dataset = ''
row_width =25

for(keyword in keyword_ls){ # loop1
    print(paste0('start ', keyword))
    ## get the meta or jack files or mm up and down regulations
    (f_gene_list_all <- list.files(jack_meta_folder, full.names=T,
                                   pattern=paste0('_regulation_',keyword, '.tsv')))
    print(f_gene_list_all)
    
    ##************************##
    ## plot output dir
    ##************************##
    plotdir <- paste0(plt_out, Sys.Date(), '/not_clustered/')
    print(paste0('plotdir is ', plotdir))
    dir.create(plotdir, showWarnings=F, recursive=T)
    
    for (disease_phase in phase_ls){ #loop2 for different phases
        ##************************##
        ## load the expression data and get annotations for each phase
        ##************************##
        mm_phase <- paste0(mm_dir,'/', disease_phase,'/')
        rdata <- grep(rdata_keyword, list.files(mm_phase, full.names = T), value = T)
        if(length(rdata) !=1){
            stop(paste0(mm_phase, ' HAVE ', rdata))
        }
        load(rdata)
        
        ## get annotation
        annotation <- data.frame(ProbeName=row.names(array_dat), GeneSymbols=row.names(array_dat))
        if(reoder_samples){
            ## reorder the samples to wt and disease
            array_design <- array_design[, c('Sample', 'Disease_stage')]
            array_design <- array_design[with(array_design, order(array_design$Disease_stage, decreasing = F)), ]
            array_design$Disease_stage <- factor(array_design$Disease_stage, levels = c('WT', 'Disease'))
        }

        array_dat <- array_dat[, array_design$Sample]
        
        ##************************##
        # get the top up and down gene list
        ##************************##
        gene_list <- NULL
        for(regulation in c('up', 'down')){ ## loop 3get up and down top genes
            (f_gene_list <-grep(regulation, f_gene_list_all, value = T))  ## get the right regulation
            (f_gene_list <-grep(disease_phase, f_gene_list, value = T)) ## get the right disease phase
            #***************************
            #' loop for each disease phase and correcponding studies
            ## f_gene_list: a meta or jack file
            ## threshold: top # of genes in the f_gene_list
            ## info_df: a df with File(the limma object files), Phase, Extra_phase
            #***************************
            print(f_gene_list)
            
            ## get the genes from jack or meta file
            df <- read.delim(f_gene_list,comment.char = "#")
            
            ## reorder by the score column
            (score_col <- intersect(c('adj_combined_max_p','Fisher', 'up_pval', 'down_pval'), colnames(df)))
            df <- df[order(df[, score_col]), ]
            
            gene_list <- c(gene_list, as.character(df$geneSymbol[1:threshold]))
        }# loop3
        
        
        ##************************##
        # plot the heatmap per phase
        ##************************##
        #### for figure quality
        ## for pub not to show sample names, plots are the same dimentions regardless of sample size
        # width based on the sample size, height bease on the number of probes/genes
        plt_title <- paste0(keyword,'_', disease_phase)
        p_w <- row_width * ncol(array_dat)+ row_width*13
        p_h <- 25 * length(gene_list)
        df_tmp <- topGeneHeatmap(array_dat, gene_list, array_design, annotation, 
                             plot_pre=plotdir, 
                             height = p_h, 
                             width = p_w, 
                             prefix = '',
                             plt_title =plt_title,
                             scale_row =T,
                             show_colnames =F,
                             write_df =F,
                             best_probe =T,
                             auto_w_d= F,
                             display_full_name = F,
                             legend_breaks =  -2:2, 
                             legend_labels = -2:2,
                             size_r=24, return_df = T)
        ## write the df
        genes <- apply(as.data.frame(row.names(df_tmp)), 1, function(x) unlist(strsplit(x, split = ':'))[1])
        df_tmp_gene <- cbind(data.frame(GeneSymbols = genes), df_tmp)
        f_out <- paste0(plotdir, '/', disease_phase,'_top_', threshold, '_up_and_down_genes.tsv')
        msg <- paste0('# ', Sys.Date(), '\n# ', plt_title,'\n# gene expression are quantile normalized across all studies')
        writeTable(df_tmp_gene, f_out = f_out,msg = msg)
        
        print(f_out)
    } # end loop2-phase
}# end loop1-keyword