#' 2016-07-19
#'@plotTopMetaGenes() plot significant gene of meta and jack 
#'@multiplotGenesExp() plot a given list of genes with a given expression r data from mixed model
# the gene expression value is from all probes of the gene from each study (before mixed model quantile normalization)
# home_dir <- '/home/bzhuang/'

source('helper_functions.R')
source('mixed_models/top_genes_mixed_model.R') ## get the multipPlots() function

#----------------------------------------------------#
## functions
#----------------------------------------------------#
geneExp <- function(all_array_dat_bf_aggregation, gene_ls){
    #' must load expression.Rdata from mixed models, the array_dat has gene as row names
    #' input gene_list is string of genes
    #' output probe expressions of the genes from each dataset (a list)
    ## get the probe expression of each dataset and output a list
    probe_exp <- vector('list')
    for(i in 1:length(all_array_dat_bf_aggregation)){
        probe_df <- all_array_dat_bf_aggregation[[i]]
        ## get the probes correponding to the genes
        gene_col <- grep('GeneSymbols', colnames(probe_df), value = T)
        probe_df <- probe_df[which(probe_df[, gene_col] %in% gene_ls), ]%>%droplevels()
        probe_exp[[i]] <- probe_df
        names(probe_exp)[i] <- names(all_array_dat_bf_aggregation)[i]
    }
    return(probe_exp)
}



#----------------------------------------------------#
## Main functions
#----------------------------------------------------#
#' @example 
# plot_gene_dir <- "/home/bzhuang/HD_mouse_model_project/mixed_model/tmp/"
# r_obj <- "/home/bzhuang/HD_mouse_model_project/mixed_model/random_intercept/late/expression.Rdata"
# gene_ls <- 'Srm'
# multiplotGenesExp(gene_ls, r_obj, plot_gene_dir)

multiplotGenesExp <- function(gene_ls, r_obj, plot_gene_dir, prefix =''){
    #' r_obj is the expression.Rdata from mixed models, which has all_array_dat_bf_aggregation,
    #' and array_design
    #' gene_ls, a string of gene names to be plotted
    #' plot_gene_dir: outdir for plot
    #' prefix: prefix for the plot name
    dir.create(plot_gene_dir, recursive = T, showWarnings = F)
    load(r_obj)
    
    ## get probe expressions of the gene list
    probe_exp <- geneExp(all_array_dat_bf_aggregation,gene_ls)
    
    if(prefix !=''){
        prefix <- paste0(prefix,'_')
    }
    
    # for each gene, plot probe expression, colored by genotype
    for(i in 1:length(gene_ls)){
        ## plot  gene expressions
        gene <- gene_ls[i]
        # colored by genotype
        (out_f_g <- paste0(plot_gene_dir,prefix, 'genotype_',i,'_',gene,'.png'))
        multipPlots(gene, probe_exp, array_design,out_f_g,
                    w_per_data=400,
                    plt_h = 800,
                    color_factor='Genotype')
        # colored by genotype
        out_f_p <- paste0(plot_gene_dir,prefix,'probename_',i,'_',gene,'.png')
        multipPlots(gene, probe_exp, array_design,out_f_p,
                    w_per_data=400,
                    plt_h = 800,
                    color_factor='ProbeName')
    }
}


#----------------------------------------------------#
## Main function loop for meta or jack files
#----------------------------------------------------#
# ## for jackknife, meta, mixed_model, mixed_model 
# value_col <- c('adj_combined_max_p', 'Fisher','down_pval','up_pval')
#' input the jack_meta_folder, and Rdata from mixed model, plot the top genes by probe names and genotype
#' (only diseasa vs wt)
#' loop for each phase, up or down regulation, and meta or jackknife file for a disease
#' @param plot_out: plot dir, within the folder, a phase, and up or down regulation subfolder will be create
#' 


plotTopMetaGenes <- function(jack_meta_folder, mm_dir, plot_out, threshold = 100, 
                             phase_ls=c('early', 'late'),
                             regulation_ls=c('up', 'down')){

    for (phase in phase_ls){ ## loop 1
        
        ## get the expression r data from mixed model
        print(phase)
        datadir = paste0(mm_dir,'/', phase, '/')  ## get mixed model Rdata(intercept or slope are the same)
        mm_f <- grep('expression.Rdata', list.files(datadir, full.names = T), value = T)
        if (length(mm_f) != 1){
            stop(paste0('no expression.Rdata found in ', datadir))
        }
        print(paste0('load expression data: ', mm_f))
        
        for(regulation in regulation_ls){ ## loop 2
            print(regulation)
            for(meta_jack in c('meta', 'jackknife')){## loop 3
                ## select the meta or jack file
                (mj_f <- list.files(jack_meta_folder,full.names = T, 
                                    pattern = paste0(phase, '.*', regulation, '.*', meta_jack, '.*.tsv')))
                
                print(mj_f)
                ## check if the file is low_exp_rm input
                if(grepl('low_exp_rm', mj_f)){
                    low_exp_rm = '_low_exp_rm'
                }else{low_exp_rm = ''}
                
                ## read meta or jack files
                df <- read.delim(mj_f, comment.char = '#')
                (index_col<- intersect(colnames(df), c('Fisher',"adj_combined_max_p")))
                rownames(df) <- df$geneSymbol
                
                ## get the top meta/jack genes
                (gene_ls<- as.character(df$geneSymbol[order(df[, index_col])][1:threshold]))
                
                ## plot the genes
                (plot_gene_dir <- paste0(plot_out, '/', phase, '/', regulation, low_exp_rm,'/'))
                dir.create(plot_gene_dir, recursive = T, showWarnings = F)
                (plot_pre <- paste(phase, regulation, meta_jack, sep = '_'))
                multiplotGenesExp(gene_ls, r_obj=mm_f, plot_gene_dir=plot_gene_dir, prefix =plot_pre)
                
                ## log
                msg <- paste0('# ', Sys.Date(), 
                              '\n# Input for ', phase,' ', regulation, '-regulated',low_exp_rm,
                              '\n# Input expression Rdata: ', mm_f,
                              '\n# Input gene list file: ', mj_f,
                              '\n# Plot out dir: ', plot_gene_dir,
                              '\n# Top ', threshold,' genes')
                        
                writeTable(df=data.frame(geneSymbol = gene_ls), msg = msg, 
                           f_out = paste0(plot_gene_dir, plot_pre, '_input.tsv'))
                } ## loop 3 end
        }
    }
}

        