# 2017-03-01 # thesis plots
# plot the expression of a gene (or gene list) from each study
# with raw expression, study corrected or study and MGP corrected



 



library(dplyr)
library(ggplot2)
library(HelperFunctions)




arrayData <- function(disease, phase, kuhn =F, geno_f=NULL){
    
    #' load array data: raw: quantile normalized
    #' MGP: study and MGP are corrected and subtracted from expression
    #' stduy: study is corrected and subtracted from expression
    #' geno_f: if add info like model types and repeat length
    
    ## load R data (before correction)
    model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'  ## to specify which model folder
    rdata_keyword <- 'expression'
    (rdata <- paste0(home_dir, disease, '_mouse_model_project/mixed_model/random_intercept', 
                     model_keyword, '/', phase, 
                     '/',rdata_keyword,'.Rdata'))
    load(rdata)
    array_dat_raw <- array_dat
    
    
    ## load R data (after MGP and study correction)
    model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'  ## to specify which model folder
    rdata_keyword <- 'mixed_model_results_exp_corrected'
    
    (rdata <- paste0(home_dir, disease, '_mouse_model_project/mixed_model/random_intercept', 
                     model_keyword, '/', phase, 
                     '/',rdata_keyword,'.Rdata'))
    load(rdata)
    array_dat_MGP <- array_dat
    
    
    ## load R data (before correction)
    model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder
    rdata_keyword <- 'mixed_model_results_exp_corrected'
    (rdata <- paste0(home_dir, disease, '_mouse_model_project/mixed_model/random_intercept', 
                     model_keyword, '/', phase, 
                     '/',rdata_keyword,'.Rdata'))
    load(rdata)
    array_dat_study <- array_dat
    
    
    if(kuhn){
        array_design$Study <- as.character(array_design$Study)
        array_design$Study <- gsub('GSE10202|GSE7958|GSE9375|GSE9857','Kuhn et al.',array_design$Study)
        array_design$Study <- as.factor(array_design$Study)
    }
    
    
    ## if add info like model types and repeat length
    if(!is.null(geno_f)){
        print('add model types and repeat lengths')
        df_geno <- read.delim(geno_f, comment.char = '#')
        
        need_col <- intersect(c('Sample', 'Model_types', 'repeat_length'), colnames(df_geno))
        
        array_design <- noWarnings(left_join(array_design, df_geno[, need_col ]))
    }
    
    return(list(array_dat_raw, array_dat_MGP, array_dat_study,array_design))
}


processDF <- function(array_dat, gene_ls, disease,array_design){
    df <- NULL
    for(gene in gene_ls){
        df_t <- as.data.frame(t(array_dat[gene, ]))
        colnames(df_t) <- 'expression'
        df_t$Sample <- rownames(df_t)
        df_t$gene <- gene
        if(is.null(df)){
            df <- df_t
        }else{
            df <- rbind(df, df_t)
        }
    }
    
    ## remove NA and join designs
    df <- noWarnings(left_join(df, array_design)) %>% droplevels()
    
    ## remove NA expression
    df <- df[which(!is.na(df$expression)),] %>% droplevels()
    
    ## change disease to names
    df$Disease_stage <- as.character(df$Disease_stage)
    df$Disease_stage[which(df$Disease_stage =='Disease')] <- disease
    df$Disease_stage <- factor(df$Disease_stage, levels = c('WT', disease))
    ## order studies
    df$Study <- factor(df$Study, levels = sort(levels(df$Study)))
    ## order genes
    df$gene <- factor(df$gene, levels = gene_ls)
    
    
    return(df)
}



plotGeneExprMGP <- function(df, 
                         one_plot_font_size =10,
                         x_angle=0,
                         jitter_w = 0.1,
                         jitter_h = 0.1,
                         plt_title ='',
                         ylab_title = 'Expression',
                         y_axis_null =F,
                         facet_1,
                         facet_2,
                         box_size = 0.3,
                         dot_size = 0.8){
    

    theme_MGP = function(fontSize,x_angle){
        list(cowplot::theme_cowplot(fontSize),
             ylab(ylab_title),
             theme(axis.text.x = element_text(angle = x_angle),
                   axis.title.x = element_blank()),
             geom_violin( # color="#C4C4C4" ,
                 # fill="#C4C4C4"
             ),
             geom_boxplot(color = 'black',width=0.3,fill = 'lightblue',outlier.size = 0,alpha = 0.5,lwd = 0.3),
             geom_jitter(color = 'black',fill = 'black',size = 0.3,width = 0.1)
        )
    }
    #' box plot gene expression with input df, before or after MGP (can be multiple genes), separated by study or not
    #' box plot gene expression with input df, before and after MGP (only 1 gene), separated by study or not
    
    #facet_1,  "~", facet_2: combos: . ~ Study; .~ Model_types; mgp ~ Study; mgp ~ Model_types
    
    box_colour <- c('indianred1','dodgerblue')
    (names(box_colour) <- c(levels(df$Disease_stage)))
    
    p_one_plot <- ggplot(df,aes_string(x = 'Disease_stage', y = 'expression', colour='Disease_stage',fill = 'Disease_stage')) + 
        theme_MGP(one_plot_font_size,0)
    
    # p_one_plot <- ggplot(df,aes_string(x = 'Disease_stage', y = 'expression', colour='Disease_stage')) + 
    #     theme_bw() + 
    #     ylab(ylab_title)+
    #     # scale_y_continuous(breaks = seq(from = -2, to =14, by=1)) +
    #     theme(axis.title.x=element_blank(),
    #           axis.title.y=element_text(size = one_plot_font_size),
    #           title = element_text(size = one_plot_font_size, colour = 'black'),
    #           axis.text.x = element_text(size = one_plot_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
    #           axis.text.y=element_text(size = one_plot_font_size),
    #           strip.text.x = element_text(size = one_plot_font_size),
    #           strip.text.y = element_text(size = one_plot_font_size))
    if(y_axis_null){
        p_one_plot <- p_one_plot + theme(axis.title.y=element_blank())
    }

    ## plot 1: plain box plot with points
    p_box <- p_one_plot + # geom_boxplot(lwd=box_size, outlier.shape = NA)+
        # geom_jitter(width = jitter_w, height = jitter_h,size = dot_size, alpha = 0.5) +
        theme(legend.position="none")+
        scale_colour_manual(values = box_colour) + 
        scale_fill_manual(values= box_colour)
    
    
    
    
    cmd <- paste0("p_box_final <- p_box + facet_grid(", facet_1,  "~", facet_2, ",drop = TRUE,scale='free')")
    eval(parse(text = cmd))
    
    return(p_box_final)
}


mainPlot <- function(gene_ls, array_dat_raw, array_dat_MGP, 
                     disease, array_design,f_out_dir,one_plot_font_size =9,
                     y_axis_null =F,
                     return_df = F,
                     save_p = T,
                     ylab_title = 'Expression',
                     ...){
    

    # ... more arguments in plotGeneExprMGP()
    df_raw <- processDF(array_dat_raw, gene_ls, disease,array_design)
    df_mgp <- processDF(array_dat_MGP, gene_ls, disease,array_design)
    df_raw$mgp <- 'before MGP'
    df_mgp$mgp <- 'after MGP'
    df_gene <- rbind(df_raw, df_mgp)
    df_gene$mgp <- factor(df_gene$mgp, levels = c('before MGP', 'after MGP'))
    
    for(gene in gene_ls){
        print(gene)
        df <- df_gene[which(df_gene$gene == gene), ]%>%droplevels()
        
        
        # facet_1_ls = c('mgp', 'mgp')
        facet_1_ls = 'mgp'
        facet_2_ls = 'Study'
        # facet_2_ls = c('Study', 'Model_types')
        for(i in 1: length(facet_1_ls)){ # plot for different facet for mgp
            facet_1 <- facet_1_ls[i]
            facet_2 <- facet_2_ls[i]
            (p_box <- plotGeneExprMGP(df, plt_title = '',ylab_title = ylab_title, 
                                      one_plot_font_size =one_plot_font_size,
                                      y_axis_null = y_axis_null,
                                      facet_1 = facet_1,
                                      facet_2 = facet_2, ...))
            if(save_p){
                (f_plot <- paste0(f_out_dir,disease,'_', phase,'_',levels(df$gene),'_',facet_1, '_', facet_2,'.png'))
                ggsave(filename = f_plot, plot = p_box, width = 6.45, height =3.25 , units ='in')
            }
            
        } # end loop plot for different facet for mgp
        
        
        ## plot only before or after
        
        # facet_1_ls = c('.', '.')
        # facet_2_ls = c('Study', 'Model_types')
        facet_1_ls = '.'
        facet_2_ls = 'Study'
        
        for(x in c("after MGP", "before MGP")){
            if(x == 'after MGP'){ylab_title_x = 'Corrected Expression'
            }else{
                ylab_title_x = ylab_title 
            }
            df_x <- df[which(df$mgp == x), ]%>%droplevels()
            for(i in 1: length(facet_1_ls)){ # plot for different facet for after or before mgp only
                facet_1 <- facet_1_ls[i]
                facet_2 <- facet_2_ls[i]
                (p_box <- plotGeneExprMGP(df_x, plt_title = '',ylab_title = ylab_title_x, 
                                          one_plot_font_size =16,
                                          y_axis_null = y_axis_null,
                                          facet_1 = facet_1,
                                          facet_2 = facet_2, 
                                          box_size = 0.6,
                                          dot_size = 1, ...))
                if(save_p){
                    (f_plot <- paste0(f_out_dir,disease,'_', phase,'_',levels(df$gene),'_', facet_2,'_',x,'.png'))
                    if(nlevels(df_x[,facet_2])>5){
                        print('plot longer pic')
                        p_width = 12
                    }else{ p_width =8}
                    ggsave(filename = f_plot, plot = p_box, width = p_width, height =5.3 , units ='in')
                }
                
            } 
        }# end loop plot for different facet for after or before mgp only
        
        
        ### this is for HD, plot repeat length vs. expression
        
        
        if(disease == 'HD'){
            x= "after MGP"
            df_x <- df[which(df$mgp == x), ]%>%droplevels()
            df_x$Model_types <- as.character(df_x$Model_types)
            df_x$Model_types[which(df_x$Genotype =='WT')] <- 'WT'
            df_x$Model_types <- as.factor(df_x$Model_types)
            
            p <- ggplot(df_x,aes_string(x = 'repeat_length', y = 'expression')) + 
                theme_bw() + 
                geom_smooth() +
                xlab('Repeat Length')+
                ylab('Corrected Expression')+
                theme(title = element_text(size = one_plot_font_size, colour = 'black'),
                      axis.text.x = element_text(size = one_plot_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
                      axis.text.y=element_text(size = one_plot_font_size),
                      strip.text.x = element_text(size = one_plot_font_size),
                      strip.text.y = element_text(size = one_plot_font_size))+
                geom_jitter(aes_string(colour='Model_types'), width = 0.2, height = 0.2,size = 1.5, alpha = 0.8)
            
            (f_plot <- paste0(f_out_dir,disease,'_', phase,'_',levels(df$gene),'_repeat_length_',x,'.png'))
            ggsave(filename = f_plot, plot = p, width = 6.45, height =3.25 , units ='in') 
        }
        
    }# end of each gene
    
    if(return_df){
        return(df_gene)
    }
}


## run functions
# #############################################################################
# f_out_dir <- paste0('../../results/ND_results/gene_expr_before_after_MGP/', Sys.Date(), '/')
# dir.create(f_out_dir, showWarnings = F,recursive = T)
# 
# 
# #####################################
# disease <- 'AD'
# phase <- 'early'
# #####################################
# x <- arrayData(disease,phase)
# array_dat_raw <- x[[1]]
# array_dat_MGP <- x[[2]]
# array_dat_study <- x[[3]]
# array_design <- x[[4]]
# 
# ### compare before and after MGP ( 1 gene a plot) ## plot
# gene_ls <- c('Sqle','Msmo1')
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir)
# 
# 
# 
# 
# #####################################
# disease <- 'AD'
# phase <- 'late'
# #####################################
# x <- arrayData(disease,phase)
# array_dat_raw <- x[[1]]
# array_dat_MGP <- x[[2]]
# array_dat_study <- x[[3]]
# array_design <- x[[4]]
# 
# ### compare before and after MGP ( 1 gene a plot) ## plot
# gene_ls <- c('Trem2','Msmo1')
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir)
# # 
# # ####### other plots
# # gene_ls <- c('Trem2', 'Inppl1','Tyrobp')
# # gene_ls <- c('Man2b1', 'Idh1', 'Cd68', 'Msmo1', 'Sqle', 'Cidea', 'Nsdhl', 'Fdps', 'Wdr82', 'Trem2', 
# #              'Reep1', 'Slc22a4', 'Slc14a1', 'Twsg1', 'Tpst1', 'Brix1', 'Tyrobp', 'Fahd1', '0610007P14Rik', 'Rab6a')
# # 
# # gene_ls <- gene_ls[1:10]
# # 
# # gene_ls <- c('Msmo1','Sqle','Nsdhl') # AD late up
# # gene_ls <- c('Fdps')
# # 
# # gene_ls <- c('Trem2')
# # 
# # ## before correction with raw data
# # df_raw <- processDF(array_dat_raw, gene_ls, disease,array_design)
# # plotGeneExprMGP(df_raw, by_study = T, plt_title = 'raw', ylab_title = 'Expression')
# # plotGeneExprMGP(df_raw, by_study = F,plt_title = 'raw',ylab_title = 'Expression')
# # 
# # ### after study correction
# # df_study <- processDF(array_dat_study, gene_ls, disease,array_design)
# # plotGeneExprMGP(df_study, by_study = T, plt_title = 'after study',ylab_title = 'Corrected Expression')
# # plotGeneExprMGP(df_study, by_study = F, plt_title = 'after study',ylab_title = 'Corrected Expression')
# # 
# # 
# # ### after MGP and study correction
# # df_mgp <- processDF(array_dat_MGP, gene_ls, disease,array_design)
# # plotGeneExprMGP(df_mgp, by_study = T, plt_title = 'after MGP and study',ylab_title = 'Corrected Expression')
# # plotGeneExprMGP(df_mgp, by_study = F, plt_title = 'after MGP and study',ylab_title = 'Corrected Expression')
# 
# 
# 
# 
# 
# 
# #############################################################################
# 
# 
# 
# #####################################
# disease <- 'HD'
# phase <- 'early'
# #####################################
# x <- arrayData(disease,phase)
# array_dat_raw <- x[[1]]
# array_dat_MGP <- x[[2]]
# array_dat_study <- x[[3]]
# array_design <- x[[4]]
# 
# ### compare before and after MGP ( 1 gene a plot) ## plot
# gene_ls <- c('Coa3', 'Ndufs3')
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir)
# 
# #####################################
# disease <- 'HD'
# phase <- 'late'
# #####################################
# x <- arrayData(disease,phase,kuhn =T)
# array_dat_raw <- x[[1]]
# array_dat_MGP <- x[[2]]
# array_dat_study <- x[[3]]
# array_design <- x[[4]]
# 
# ### compare before and after MGP ( 1 gene a plot) ## plot
# gene_ls <- c('Isl1','Tac1', 'Ddit4l', 'Nrep')
# gene_ls <- c('Ddit4l', 'Nrep')
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,one_plot_font_size =9,y_axis_null =F)
# 
# 
# ####explore
# 
# gene_ls <- c('Enpp6') ## up overlap in early and late
# mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir, save_p = F, return_p = T)

# 
# 
# ####### plot
# 
# gene_ls <- c('Doc2b') # top up (but down before MGP)
# gene_ls <- c('Lrrn3','Isl1') # top up 1-2 agreement
# gene_ls <- c('Smoc1')  # top up, 5 agreement
# 
# 
# ## top 20 up
# gene_up <- c('Doc2b', 'Lrrn3', 'Isl1', 'Cbx8', 'Smoc1', 
#              'Klhl13', 'Nagk', 'P2ry1', 'Robo1', 'Meis1', 
#              'Sap130', 'Pcdhb22', 'Fam196b', 'Sertad4', 'Tac1', 
#              'Wbp5', 'Psme1', 'Hexim1', 'Dcc', 'Prkch')
# 
# 
# gene_up <- c('Doc2b', 
#              'Klhl13', 'P2ry1',
#              'Sap130', 'Pcdhb22', 'Sertad4', 'Tac1', 
#              'Dcc', 'Prkch') ## no agreement
# 
# gene_ls <- gene_up[1:3]
# gene_ls <- gene_up[4:6]
# 
# 
# 
# ## before correction with raw data
# df_raw <- processDF(array_dat_raw, gene_ls, disease,array_design)
# plotGeneExprMGP(df_raw, by_study = T, plt_title = 'raw', ylab_title = 'Expression')
# plotGeneExprMGP(df_raw, by_study = F,plt_title = 'raw',ylab_title = 'Expression')
# 
# ### after study correction
# df_study <- processDF(array_dat_study, gene_ls, disease,array_design)
# plotGeneExprMGP(df_study, by_study = T, plt_title = 'after study',ylab_title = 'Corrected Expression')
# plotGeneExprMGP(df_study, by_study = F, plt_title = 'after study',ylab_title = 'Corrected Expression')
# 
# 
# ### after MGP and study correction
# df_mgp <- processDF(array_dat_MGP, gene_ls, disease,array_design)
# plotGeneExprMGP(df_mgp, by_study = T, plt_title = 'after MGP and study',ylab_title = 'Corrected Expression')
# plotGeneExprMGP(df_mgp, by_study = F, plt_title = 'after MGP and study',ylab_title = 'Corrected Expression')
# 
# 
