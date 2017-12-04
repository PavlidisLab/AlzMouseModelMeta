# 2016-9-02
#' updated: 2016-11-07
#' updated: 2016-11-8: add model types
#' heatmap plot for mixmdel results, expression is intersect subtracted and corrected for gender
#' each study will have a plot, if a study has 2 or mouse models, a sub plot for each model will be plotted too

###########

source('helper_functions.R')

####################################
####### helper functions###########
####################################
getGeneList <- function(opposite_phase,  f_gene_list_all, disease_phase, threshold){
    # return a gene list
    gene_list <- NULL
    f_gene_list_ls <- NULL
    for(regulation in c('up', 'down')){ ## loop 3 get up and down top genes
        if(opposite_phase){
            (f_gene_list <-grep(regulation, f_gene_list_all, value = T))  ## get the right regulation
            (op_phase <- setdiff(c('early', 'late'), disease_phase))
            msg = paste0('Input expression is from ', disease_phase, ' . Top genes are from ', op_phase)
            (f_gene_list <-grep(op_phase, f_gene_list, value = T)) ## get the right disease phase
        }else{
            (f_gene_list <-grep(regulation, f_gene_list_all, value = T))  ## get the right regulation
            (f_gene_list <-grep(disease_phase, f_gene_list, value = T)) ## get the right disease phase
            msg=''
        }
        
        #***************************
        #' loop for each disease phase and correcponding studies
        ## f_gene_list: a meta or jack file
        ## threshold: top # of genes in the f_gene_list
        ## info_df: a df with File(the limma object files), Phase, Extra_phase
        #***************************
        print(f_gene_list)
        f_gene_list_ls <- c(f_gene_list_ls, f_gene_list)
        
        ## get the genes from mixe model results (or jacknife mixed model results)
        df <- read.delim(f_gene_list,comment.char = "#")
        rownames(df) <- df$geneSymbol
        
        ## remove specified genes
        if(!is.null(rm_gene_list)){
            print(intersect(rownames(df), rm_gene_list))
            gene_rm <- which(rownames(df) %in% intersect(rownames(df), rm_gene_list))
            if(length(gene_rm)>0){
                df <- df[-gene_rm, ]
            }
        }
        
        ## reorder by the score column:avg_rank_percent is for jackknife
        (score_col <- intersect(c('adj_combined_max_p','Fisher', 'up_pval', 'down_pval','avg_rank_percent'), colnames(df)))
        df <- df[order(df[, score_col]), ]
        
        gene_list <- c(gene_list, as.character(df$geneSymbol[1:threshold]))
    }# loop3
    return(list(gene_list, msg,f_gene_list_ls))
}



plotGeneHeat <- function(df_exp,df_design, gene_annotation,gene_list,plotdir, 
                         disease_phase,plt_title,
                         show_rownames, row_width = 25,
                         legend = c('Genotype', 'Dataset'),cluster_rows =F, cluster_cols =F, ...){
    if(show_rownames){
        (p_w <- row_width * ncol(df_exp)+ max(nchar(gene_list))*15+116)
    }else{
        (p_w <- row_width * ncol(df_exp)+ 116)
    }
    
    p_h <- 25 * length(gene_list)
    df_tmp <- topGeneHeatmap(df_exp, gene_list, df_design, gene_annotation, 
                             plot_pre=plotdir, 
                             height = p_h, 
                             width = p_w, 
                             prefix = '',
                             plt_title =plt_title,
                             scale_row =T,
                             show_colnames =F,
                             show_rownames = show_rownames,
                             write_df =F,
                             best_probe =T,
                             auto_w_d= F,
                             display_full_name = F,
                             legend_breaks =  -2:2, 
                             legend_labels = -2:2,
                             size_r=24, 
                             legend = legend, ## must include genotype for label
                             return_df = T, 
                             cluster_rows = cluster_rows, 
                             cluster_cols = cluster_cols, ...) 
    
}

##########################
## MAIN function
##########################


plotThesisHeatGenes <- function(mm_dir,mm_rdata_dir,
                                plt_out,
                                phase_ls, 
                                threshold,
                                rdata_keyword = 'mixed_model_results_exp_corrected.Rdata',  # choose which rdata to get from: corrected or not, 
                                rm_gene_list = NULL,
                                top_genes =T, opposite_phase =F, 
                                gene_list =NULL,
                                row_width =25,
                                legend_ls= list(c('Genotype', 'Model_types')),
                                #legend_ls= list(c('Genotype', 'Study'), c('Genotype', 'Original_genotype', 'Study')),
                                plot_indi = T, 
                                cluster_rows =F, cluster_cols =F, 
                                geno_f =NULL, 
                                plot_suffix =NULL, 
                                plotdir = NULL,
                                ...){
    #' geno_f=paste0('../../doc/', disease, '_mouse_dataset_doc/dataset_info_genotypes.tsv')
    ## get the mm up and down regulations
    ## geno_f is the file with model_types
    #' ... more for heatmap settings
    #' plotdir = NULL, add disease phase to to plot outdir, else just use plt_out
    (f_gene_list_all <- list.files(mm_dir, full.names=T,
                                   pattern=paste0('regulation.*.tsv')))
    print('INPUT gene list')
    print(f_gene_list_all)
    
    
    for (disease_phase in phase_ls){ #loop2 for different phases
        ##************************##
        ## load the expression data and get gene_annotations for each phase
        ##************************##
        mm_phase <- paste0(mm_rdata_dir,'/', disease_phase,'/')
        if(is.null(plotdir)){
            (plotdir <- paste0(plt_out,'/', disease_phase,'/'))
        }else{
            (plotdir <- plt_out)
        }

        print(paste0('plotdir is ', plotdir))
        dir.create(plotdir, showWarnings=F, recursive=T)
        rdata <- grep(rdata_keyword, list.files(mm_phase, full.names = T), value = T)
        if(length(rdata) !=1){
            stop(paste0(mm_phase, ' HAVE ', rdata))
        }
        print(paste0('Input expression R data is ', rdata))
        load(rdata)
        
        ## get gene_annotation
        gene_annotation <- data.frame(ProbeName=row.names(array_dat), GeneSymbols=row.names(array_dat))
        
        ## order expression by sample
        array_dat <- array_dat[, array_design$Sample]
        
        ##************************##
        # OPTION1: get the top up and down gene_list (opposite_phase = T: plot expression of genes in early from the top late genes) of 
        #' the same regulation
        ##************************##
        if(top_genes){
            print('GET TOP GENES')
            x <- getGeneList(opposite_phase,  f_gene_list_all, disease_phase,threshold)
            gene_list <- x[[1]]
            msg <- x[[2]]
            f_gene_list <- x[[3]]
        }else{
            ##************************##
            # OPTION2: use the input gene_list
            ##************************##
            msg = '\nMANUAL INPUT GENE LIST\n'
            print(msg)
            print(gene_list)
            f_gene_list <- NULL
            
        }
        
        ## remove .1 .2 from study names
        array_design$Study <- gsub('\\.1|\\.2', '', as.character(array_design$Study))
        
        
        ##************************##
        ### plot for all studies ## used in thesis
        ##************************##
        print('Plot a heatmap for all samples')
        gene_list <- unique(intersect(gene_list, rownames(array_dat)))
        all_exp <- array_dat[gene_list, ]%>%droplevels()  ## get the genes
        all_design <- array_design
        all_design$Study <- as.factor(all_design$Study)
        all_design$Original_genotype <- all_design$Genotype
        all_design$Genotype <- as.character(all_design$Disease_stage)
        all_design$Genotype[which(all_design$Genotype == 'Disease')] <- disease
        all_design$Genotype <- as.factor(all_design$Genotype)
        all_design$Timepoint <- '1m'
        ## must remember to relevel
        all_design$Genotype <- factor(all_design$Genotype, levels = c('WT', disease))
        
        ##################################
        if(!is.null(geno_f)){
            print('Input the model types')
            df_geno <- read.delim(geno_f, comment.char = '#')
            ## get disease
            df_geno <- df_geno[which(df_geno$Phase == disease_phase), ]
            df_geno <- orderCol(df_geno, 'Order')
            
            all_design <- noWarnings(left_join(all_design, df_geno))
            all_design$Model_types <- as.character(all_design$Model_types)
            all_design$Model_types[which(all_design$Genotype == 'WT')] <- 'WT'
            
            ## get the model order first by study, then by models
            (target_order = as.character(unique(df_geno$Study)))
            all_design$Study <- factor(all_design$Study, levels = c(target_order))
            
            
            target_order = as.character(unique(df_geno$Model_types))
            all_design$Model_types <- factor(all_design$Model_types, levels = c('WT', target_order))
            
            all_design <- all_design %>% arrange(Study, Model_types)%>%droplevels()
            
            ## reorder the expr too
            all_exp <- all_exp[, all_design$Sample]
            
            #legend_ls <- c(legend_ls, list(c('Genotype', 'Model_types','Study'), c('Genotype', 'Model_types')))
            writeTable(all_design, f_out = paste0(plotdir, disease,'_',disease_phase,'_design.tsv'))
        }
        
        
        #################################
        
        ## plotting
        for(legend_factors in legend_ls){
            (plt_title <- paste0(c(disease,disease_phase, legend_factors, plot_suffix), collapse = '_'))
            plotGeneHeat(all_exp,all_design, gene_annotation,gene_list,plotdir, 
                         disease_phase,plt_title,
                         show_rownames = T, row_width = 10, legend = legend_factors,
                         cluster_rows = cluster_rows, 
                         cluster_cols = cluster_cols, ...)
        }
        
        msg <- paste0('#', Sys.time(), ': \n#', msg,'\n#\n#Input expression data: \n#', rdata, '\n#\n#Input gene ranking lists:\n#',
                      paste0(f_gene_list, collapse = '\n#'))
        all_exp_out <- cbind(data.frame(geneSymbol = row.names(all_exp)), all_exp)
        writeTable(df =all_exp_out, f_out = paste0(plotdir, disease,'_',disease_phase,plot_suffix,'_inputdata.tsv'), msg = msg)
        
        ##************************##
        ### plot for each study
        ##************************##
        if(plot_indi){
            ## add timepoint in the study if not there already
            if(!length(grep('_months', array_design$Study))>0){
                array_design$Study <- paste0(array_design$Study, '_', array_design$Timepoint)
                array_design$Study <- as.factor(array_design$Study)
            }
            for(i in 1:length(levels(array_design$Study))){
                study <- levels(array_design$Study)[i]
                print(study)
                df_design <- filterContain(array_design, column = 'Study', value = study)
                
                ## change the label of disease in the plot label
                df_design$Genotype <- as.character(df_design$Disease_stage)
                df_design$Genotype[which(df_design$Genotype == 'Disease')] <- disease
                df_design$Genotype <- as.factor(df_design$Genotype)
                
                
                ## must remember to relevel
                df_design$Genotype <- factor(df_design$Genotype, levels = c('WT', disease))
                
                df_exp <- array_dat[gene_list, df_design$Sample]%>%droplevels()  ## get the genes 
                ##check how many genotypes
                (geno_ls <- levels(df_design$Original_genotype))
                
                if(i==1){ ## the first plot with gene name
                    (plt_title <- paste0('rownames_', disease_phase, '_', study))
                    plotGeneHeat(df_exp,df_design, gene_annotation,gene_list,plotdir, 
                                 disease_phase,plt_title,
                                 show_rownames = T, legend = c('Genotype'), ...)
                }
                
                if(length(geno_ls) ==2){
                    ## plot
                    (plt_title <- paste0(disease_phase, '_', study))
                    plotGeneHeat(df_exp,df_design, gene_annotation,gene_list,plotdir, 
                                 disease_phase,plt_title,
                                 show_rownames = F, legend = c('Genotype'), ...)
                }else{ 
                    ######
                    ## if more than 2 mouse models: plot all and plot each mouse model against WT
                    ######
                    (mouse_model <- setdiff(geno_ls, 'WT'))
                    ## plot all
                    (plt_title <- paste0(disease_phase, '_', study,'_all'))
                    plotGeneHeat(df_exp,df_design, gene_annotation,gene_list,plotdir, 
                                 disease_phase,plt_title,
                                 show_rownames = F, legend = c('Genotype'), ...)
                    ## plot each genotype
                    for(j in mouse_model){
                        df_design_tmp <- df_design[which(df_design$Original_genotype %in% c('WT', j)), ]%>%droplevels()
                        df_exp_tmp <- df_exp[, df_design_tmp$Sample]%>%droplevels()  ## get the genes 
                        (plt_title <- paste0(disease_phase, '_', study,'_',j))
                        plotGeneHeat(df_exp_tmp,df_design_tmp, gene_annotation,gene_list,plotdir, 
                                     disease_phase,plt_title,
                                     show_rownames = F, legend = c('Genotype'), ...)
                    }
                }
            } 
        } ## end if plot individual plots
        
    } # end loop2-phase
    
}


