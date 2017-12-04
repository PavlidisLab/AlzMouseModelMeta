# 2016-08-23 for thesis figure
#' plot a figure with top up and down regulated genes from MM, the expression
#' the expression is from normalized mm model

ggColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## plot size in inches
plt_w=7
plt_h=7.5



for(keyword in keyword_ls){ # loop1
  print(paste0('start ', keyword))
  ## get the meta or jack files or mm up and down regulations
  (f_gene_list_all <- list.files(jack_meta_folder, full.names=T,
                                 pattern=paste0('_regulation_',keyword, '.tsv')))
  print(f_gene_list_all)
  
  
  for (disease_phase in phase_ls){ #loop2 for different phases
    ##************************##
    ## load the expression data and get annotations for each phase
    ##************************##
    mm_phase <- paste0(mm_dir,'/', disease_phase,'/')
    rdata <- grep('expression.Rdata', list.files(mm_phase, full.names = T), value = T)
    if(length(rdata) !=1){
      stop(paste0(mm_phase, ' HAVE ', rdata))
    }
    load(rdata)
    
    ## get annotation
    annotation <- data.frame(ProbeName=row.names(array_dat), GeneSymbols=row.names(array_dat))
    # ## reorder the samples
    # array_design <- array_design[, c('Sample', 'Disease_stage')]
    # array_design <- array_design[with(array_design, order(array_design$Disease_stage, decreasing = F)), ]
    # array_design$Disease_stage <- factor(array_design$Disease_stage, levels = c('WT', 'Disease'))
    # array_dat <- array_dat[, array_design$Sample]
    
    ##************************##
    # get the top up and down gene list, and do the plot for both up and down
    ##************************##
    gene_list=NULL
    for(regulation in regulation_ls){ ## loop 3get up and down top genes
      (f_gene_list <-grep(regulation, f_gene_list_all, value = T))  ## get the right regulation
      (f_gene_list <-grep(disease_phase, f_gene_list, value = T)) ## get the right disease phase
      #***************************
      #' loop for each disease phase and correcponding studies
      ## f_gene_list: a meta or jack file
      ## threshold: top # of genes in the f_gene_list
      ## info_df: a df with File(the limma object files), Phase, Extra_phase
      #***************************
      #print(f_gene_list)
      
      ## get the genes from jack or meta file
      df <- read.delim(f_gene_list,comment.char = "#")
      
      ## reorder by the score column
      (score_col <- intersect(c('adj_combined_max_p','Fisher', 'up_pval', 'down_pval'), colnames(df)))
      df <- df[order(df[, score_col]), ]
      
      (gene_list <- c(gene_list,as.character(df$geneSymbol[1:threshold])))
      #df_gene <- data.frame(geneSymbol=gene_list, rank= 1:length(gene_list))
      
    }# loop3
    
    ##************************##
    # plot per phase, for both up and down
    ##************************##
    
    array_gene <- array_dat[gene_list, ]%>%droplevels()
    array_gene$geneSymbol <- row.names(array_gene)
    
    ## change it to long table and join the design
    df_exp <- melt(array_gene, value.name = 'Expression', variable.name = 'Sample',id.vars = c('geneSymbol'))
    df_exp <- noWarnings(left_join(df_exp, array_design))
    
    ## relevel
    df_exp$geneSymbol <- factor(df_exp$geneSymbol, levels = gene_list)  ## order by gene list

    df_exp$Disease_stage[which(df_exp$Disease_stage == 'Disease')] <- disease
    df_exp$Disease_stage <- factor(df_exp$Disease_stage, levels = c('WT', disease))
    df_exp$Study <- gsub('onths|ays|eeks', '', df_exp$Study)
    df_exp$Study <- gsub('_m$', 'm', df_exp$Study)
    df_exp$Study <- gsub('_d$', 'd', df_exp$Study)
    df_exp$Study <- gsub('_w$', 'w', df_exp$Study)

  
    ## assign color
    df_exp$Study <- factor(df_exp$Study)
    study <- levels(df_exp$Study)
    print(study)
    (palette <- ggColorHue(length(study)))
    (names(palette) <- study)
    
    
    ## plot and save as svg


    dir.create(plotdir, showWarnings=F, recursive=T)
    
    (f_o <- paste0(plotdir, disease,'_', disease_phase, "_top",threshold,'_', 
                    paste0(regulation_ls,collapse = '_')))
    print(paste0('plotdir is ', f_o))
    
    title=""
    p <- ggplot(na.omit(df_exp), aes_string('Disease_stage', 'Expression', color = 'Study')) + 
      geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size =1.5)+
      #stat_summary(fun.y = mean, geom = "point", shape = 4, size = 5) + 
      ggtitle(title) +
      theme_bw() +
      theme(
        axis.title.x=element_blank(),
        plot.title = element_text(size = 1),
        text =element_text(size = 14), # theme for poster size
        axis.text=element_text(size = 12)) +
      facet_wrap(~ geneSymbol, ncol = 5)+ 
      geom_line(stat="summary", fun.y="mean", aes(group=as.character(Study)), size =1,alpha = 0.5)+
      scale_colour_manual(values = palette) +  # customize color
      scale_fill_manual(values = palette)
    
    ## for png
    ggsave(filename = paste0(f_o, '.png'), plot = p, width=plt_w, height = plt_h)
    ## for svg
    ggsave(filename = paste0(f_o, '.svg'), plot = p, width=plt_w, height = plt_h)
    
    
    
    ## for the zoom in version each gene has it's own y axis
    plotlist <- vector('list')
    for(i in 1:length(gene_list)){
      gene <- gene_list[i]
      df <- filterContain(df_exp, column = 'geneSymbol', value = gene)
      
      p <- ggplot(na.omit(df), aes_string('Disease_stage', 'Expression', color = 'Study')) + 
        geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size =1.5)+
        ggtitle(gene) +
        theme_bw() +
        geom_line(stat="summary", fun.y="mean", aes(group=as.character(Study)), size =1,alpha = 0.5)+
        scale_colour_manual(values = palette) +  # customize color
        scale_fill_manual(values = palette) +
        theme(legend.position="none", 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title = element_text(size = 9),
              text =element_text(size =9), # theme for poster size
              axis.text=element_text(size = 9))
      
      plotlist[[i]] <- p
    }
    
    
    ## output the plot (and leave room for legend)
    png(filename=paste0(f_o, '_zoom_in.png'), width = plt_w-1.5, height = plt_h, units = 'in', res=400)
    multiplot(plotlist=plotlist, cols = 5,order = F)  ## must order = f, from l to r, top to bottom
    dev.off()
    print(paste0(f_o, '_zoom_in.png'))
  } # end loop2-phase
}# end loop1-keyword