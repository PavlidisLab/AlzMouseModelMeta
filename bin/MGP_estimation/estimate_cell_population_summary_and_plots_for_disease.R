#2016-10-31
#' updatecd 2016-11-03
#' 11-8: add model types
##
#' after all the estimate cell populations
# disease_ls = c('AD', 'HD', 'PD')
# phase_ls =c('early', 'late')
# source('estimate_cell_population_summary_and_plots_for_disease.R')
#

library('HelperFunctions')
library(dplyr)
library(ggplot2)

###########
#' helper functions
###########

ggColorHue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

assignSig <- function(df, adj_p){
    ## assign the significance only to mouse models samples
    #' <0.01: ***
    #' <0.05 : **
    #' <0.1:*  ## remove this
    df$significance =''
    # threshold = c(0.1, 0.05, 0.01)
    # stars = c('*', '**', '***')
    threshold = c(0.05, 0.01)
    stars = c('**', '***')
    for(i in 1: length(threshold)){
        df$significance[intersect(which(df$group !='WT'), which(df[, adj_p] <= threshold[i]))] = stars[i]
    }
    df$significance <- as.factor(df$significance)
    return(df)
}

processCellInput <- function(disease_ls, phase_ls, in_dir){
    df_disease_all = NULL
    
    for(disease in disease_ls){## loop 1 for disease
        
        ## get the estimated for early and late phases
        f_ls=NULL
        f_estimate_ls=NULL
        keyword_ls=NULL
        for(phase in phase_ls){ ## loop2 for phase
            print(phase)
            (folder_1 = grep(phase, list.dirs(in_dir), value=T))
            (folder = max(grep('201', list.dirs(folder_1, recursive=F), value = T))) ## only grep folder with 201*
            if(is.na(folder)){
                (folder = folder_1) ## grep any folder
            }
            (f = grep('_cell_proportion_estimation_scaled.tsv',list.files(folder, recursive=T, full.names = T), value=T))
            print(f)
            (f_estimate = grep('estimatefile.tsv',list.files(folder, recursive=T, full.names = T), value=T))
            f_ls=c(f_ls, f)
            f_estimate_ls=c(f_estimate_ls, f_estimate)
            keyword_ls =c(keyword_ls, phase)
            
            
            ## for phase in the disease
            df_all = NULL
            for(i in 1:length(f_ls)){
                (f= f_ls[i])
                print(i)
                print(keyword_ls[i])
                print(f)
                (f_e= f_estimate_ls[i])
                print(f_e)
                
                ## get the samples names
                df_e <- read.delim(f_e, comment.char = '#')
                df_e$Sample = rownames(df_e)
                df_e <- reshape2::melt(df_e, id.vars = c('Sample', 'groups'), variable.name="cell_type", value.name = 'PC1')
                
                df = read.delim(f, comment= '#')
                df = noWarnings(left_join(df, df_e))
                
                df$keyword = keyword_ls[i]
                if(is.null(df_all)){
                    df_all = df
                }else{
                    df_all =rbind(df_all, df)
                }
            }
        }## loop2 for phase
        df_all$disease = disease
        
        
        df_disease = df_all
        
        if(is.null(df_disease_all)){
            df_disease_all = df_disease
        }else{
            df_disease_all =rbind(df_disease_all, df_disease)
        }
    }## loop 1 for disease
    df_disease_all <- rmDup(df_disease_all) ## gseGSE64398.1_4_months WT will be duplicated for early and late
    # because they had been used in both phases
    
    
    
    ## assign the significance
    #' <0.01: ***
    #' <0.05 : **
    #' <0.1:*
    df_disease_all <- assignSig(df=df_disease_all, adj_p = 'adj_p_wt')
    
    
    ###
    ## add cell type category, and reorder cell types (and rename)
    ###
    full_name <- c('Oligo', 'Astrocyte', "Microglia","DentateGranule",'GabaSSTReln', 'StriatumCholin','Pyramidal', 'ForebrainCholin', 'Spiny', 'Microglia_deactivation', 'Microglia_activation', 
                   'A1','A2')
    new_name <- c('Oligodendrocytes','Astrocytes',"Microglia", "Dentate granule cells" ,'GABAergic cells','Cholinergic neurons', 'Pyramidal cells','Cholinergic neurons', 'Medium spiny neurons','Microglia_deactivation', 'Microglia_activation',
                  'A1_astrocyte','A2_astrocyte')
    
    df_disease_all[, 'cell_type'] <- plyr::mapvalues(df_disease_all[, 'cell_type'], from = full_name, 
                                                     to = new_name, warn_missing = F)
 
    new_type <- c(rep('glia', 3), rep('neurons', 6), rep('microglia_act', 2),rep('astrocyte_type', 2))
    
    df_disease_all$brain_cells <- df_disease_all[, 'cell_type']
    df_disease_all[, 'brain_cells'] <- plyr::mapvalues(df_disease_all[, 'brain_cells'], from = new_name, 
                                                       to = new_type, warn_missing = F)
    
    df_disease_all$plot_title <- paste0(df_disease_all$disease, ' ', df_disease_all$keyword)  ## make the timepoint in the second line
    df_disease_all$Study <- as.factor(grep('GSE', unlist(strsplit(df_disease_all$Sample, split = '_')), value =T))
    
    df_disease_all  <- as.data.frame(unclass(df_disease_all ))
    ### define outlier values of pc1 scaled:
    df_all=NULL
    for(disease in disease_ls){
        for(phase in phase_ls){
            for(cell_type in levels(df_disease_all$cell_type)){
                df_d <- df_disease_all[which(df_disease_all$cell_type == cell_type & df_disease_all$disease ==disease &df_disease_all$keyword == phase), ] %>%droplevels()
                for(group in levels(df_d$group)){
                    df <- df_d[which(df_d$group == group), ]
                    (maxv <- mean(df$PC1_scaled)+ 2*sd(df$PC1_scaled))
                    (minv <- mean(df$PC1_scaled)- 2*sd(df$PC1_scaled))
                    index <- which(df$PC1_scaled >= maxv)
                    (index <- c(index, which(df$PC1_scaled <=minv)))
                    if(length(index) >0){
                        df <- df[index, ]
                        if(is.null(df_all)){
                            df_all <- df
                        }else{
                            df_all <- rbind(df_all, df)
                        }
                    }
                } ## groups
            }## celltype
        } ## phase
    }## disease
    df_all$outlier <- df_all$Sample
    
    df_disease_all <- noWarnings(left_join(df_disease_all, df_all))
    ## change all to factors
    df_disease_all  <- as.data.frame(unclass(df_disease_all ))
    
    return(df_disease_all)
}



cellPopPlots <- function(disease_ls, phase_ls, in_dir, out_dir, x_angle = 0,one_plot_font_size = 12,
                         geno_f, # e.g. '/AD_HD_samples_model_types.tsv',
                         violin =F,
                         outlier_p =F,
                         plot_type = 'png',
                         poster =F,
                         poster_font_size =24,
                         jitter_w = 0.1,
                         jitter_h = 0,
                         outlier_rm_from_box = F){
    #' x_angle: the text angle of x labels
    #' one_plot_font_size: all the text font size
    #' outlier_rm_from_box, outliers are removed from the box plot for better seperation
    #####
    ## save disease all table
    #####
    
    ## get the processed df of cell pop info
    df_disease_all <- processCellInput(disease_ls, phase_ls, in_dir)
    
    
    #' do not remove duplicates, in AD early and late GSE64398, some WT are used in both early and late phases
    df_geno <- read.delim(geno_f, comment.char = '#') 
    
    model_order_AD <- c("Amyloid transgenic models","TAU transgenic models","KO models" ,"Other models")
    model_order_HD <- c("N-terminal transgenic", "CAG repeat KI","Human HTT exon 1 KI","Full-length transgenic")
    model_order <- c(model_order_AD, model_order_HD)
    df_geno$Model_types <- factor(df_geno$Model_types, levels = model_order)
    
    ## assign each model a unique color
    (model_palette <- c(ggColorHue(length(model_order_AD)), ggColorHue(length(model_order_HD))))
    (names(model_palette) <- model_order)
    
    df <- df_geno[c('Sample', 'Model_types', 'Phase')]%>%droplevels()
    colnames(df)[3] = 'keyword'
    
    df_disease_all <- noWarnings(left_join(df_disease_all, df))
    
    ## save table
    f_out_dir = paste0(out_dir, Sys.Date(), '/')
    dir.create(f_out_dir, recursive=T, showWarnings=F)
    (f=paste0(f_out_dir, 'cell_population_estimate_', Sys.Date(), '.tsv'))
    writeTable(df_disease_all, f_out=f)
    print(f)
    
    
    #####
    #plot
    #####
    for(disease in disease_ls){
        
        ## get the disease, and phase
        index <- Reduce(intersect, list(which(df_disease_all$disease == disease)))
        df_cell <- df_disease_all[index, ]%>%droplevels()
        
        df_cell$group <- as.factor(gsub('Disease', disease, df_cell$group))
        
        ## relevel
        (new_levels <- c('WT', setdiff(unique(df_cell$group), 'WT')))
        df_cell$group <- factor(df_cell$group, levels = new_levels)
        ## WT has 1 colour, others have the other
        df_cell$group_colour <- 'WT'
        df_cell$group_colour[which(df_cell$group !='WT')] <- 'model'
        df_cell$group_colour <- as.factor(df_cell$group_colour)
        
        
        #####
        ## for thesis, plot all studies, and all cell types in glia and neurons
        #####
        
        for(brain_cells in levels(df_cell$brain_cells)){ ## loop 2 for brain cell types
            df_brain <- df_cell[which(df_cell$brain_cells == brain_cells), ]%>% droplevels()
            df_brain$group_colour <- as.factor(df_brain$group_colour)
            print(paste0(disease, ' ', brain_cells))
            
            if(outlier_rm_from_box){
                print('### outlier are removed from the box plots')
                df_brain <- df_brain[which(is.na(df_brain$outlier)), ]
            }
            
            #### get the text label for significance
            #' see https://trinkerrstuff.wordpress.com/2012/09/01/add-text-annotations-to-ggplot2-faceted-plot/
            sig_label_df <- rmDup(df_brain[, c('cell_type', 'plot_title','significance','group')])
            sig_label_df$y = 1.05
            sig_label_df$x = sig_label_df$group
            
            #### define box color (for wt and diseaes)
            box_colour <- c('dodgerblue','indianred1')
            (names(box_colour) <- c(levels(df_brain$group_colour)))
            
            
            #### plot each cell type and each phase as a indi plot
            theme_MGP = function(fontSize,x_angle){
                list(cowplot::theme_cowplot(fontSize),
                     ylab('Marker Gene Profiles'),
                     coord_cartesian(ylim = c(-0.03, 1.10)),
                     theme(axis.text.x = element_text(angle = x_angle),
                           axis.title.x = element_blank()),
                     geom_violin( # color="#C4C4C4"# ,
                                  #fill="#C4C4C4"
                     ),
                     geom_boxplot(color= 'black',width=0.1,fill = 'lightblue')
                )
            }
            
            
            f_out_dir_indi <- paste0(f_out_dir, 'indi/')
            dir.create(f_out_dir_indi)
            indi_font_size = 14
            for(x1 in levels(df_brain$cell_type)){#loop for each cell type
                for(x2 in levels(df_brain$keyword)){ # for each phase
                    
                    df_x <- df_brain[which(df_brain$cell_type == x1 & df_brain$keyword == x2), ]%>%droplevels()
                    sig_label_df_x <- sig_label_df[which(sig_label_df$cell_type == x1 &
                                                             sig_label_df$plot_title == paste0(disease, ' ', x2) &
                                                             sig_label_df$x == disease), ]%>% droplevels()
                    
                    p_indi = ggplot(df_x,aes_string(x = 'group', y = 'PC1_scaled',color = 'group_colour', fill='group_colour')) + 
                        theme_MGP(fontSize = indi_font_size,
                                  x_angle = x_angle) + 
                        theme(legend.position="none") +  
                        scale_colour_manual(values = box_colour) + 
                        scale_fill_manual(values = box_colour) + 
                        geom_text(aes(x,y, label=significance), data=sig_label_df_x, size = 10,inherit.aes = FALSE) 
                    
                    # p_indi <- ggplot(df_x,aes_string(x = 'group', y = 'PC1_scaled', color='group_colour')) + 
                    #     ylab('Relative marker gene profiles')+
                    #     ylim(-0.1, 1.1) +
                    #     scale_y_continuous(breaks = seq(0,1,by = 0.2)) +
                    #     theme(axis.title.x=element_blank(),
                    #           axis.title.y=element_text(size = indi_font_size),
                    #           title = element_text(size = indi_font_size, colour = 'black'),
                    #           axis.text.x = element_text(size = indi_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
                    #           axis.text.y=element_text(size = indi_font_size),
                    #           strip.text.x = element_text(size = indi_font_size),
                    #           strip.text.y = element_text(size = indi_font_size))+
                    #     geom_boxplot(lwd=0.8, outlier.shape = NA)+
                    #     geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                    #     theme(legend.position="none")+
                    #     scale_colour_manual(values = box_colour) + 
                    #     geom_text(aes(x,y, label=significance), data=sig_label_df_x, size = 10, color = "black") ## add the significance dots
                    (f_plot <- paste0(f_out_dir_indi,disease,'_', x2,'_',x1, '.', plot_type))
                    ggsave(filename = f_plot, plot = p_indi, width = 3, height =4 , units ='in')
                    
                } # for each phase
                
            }#loop for each cell type


            #### plot cell type and phase in 1 plot
            p_one_plot = ggplot(df_brain,aes_string(x = 'group', y = 'PC1_scaled',color = 'group_colour', fill='group_colour')) + 
                theme_MGP(fontSize = one_plot_font_size,
                          x_angle = x_angle) + 
                scale_colour_manual(values = box_colour) + 
                scale_fill_manual(values = box_colour) + 
                facet_grid(cell_type~plot_title,drop = TRUE,scale="free", space='free') + 
                geom_text(aes(x,y, label=significance,group=NULL), data=sig_label_df, size = 10, color = "black",,inherit.aes = FALSE)
            
            # p_one_plot <- ggplot(df_brain,aes_string(x = 'group', y = 'PC1_scaled', colour='group_colour')) + 
            #     facet_grid(cell_type~plot_title,drop = TRUE,scale="free", space='free')+ 
            #     ylab('Marker gene profiles')+
            #     ylim(-0.1, 1.1) +
            #     theme(text=element_text(family = 'Arial'),
            #           #legend.position="none",
            #           axis.title.x=element_blank(),
            #           axis.title.y=element_text(size = one_plot_font_size),
            #           title = element_text(size = one_plot_font_size, colour = 'black'),
            #           axis.text.x = element_text(size = one_plot_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
            #           axis.text.y=element_text(size = one_plot_font_size),
            #           strip.text.x = element_text(size = one_plot_font_size),
            #           strip.text.y = element_text(size = one_plot_font_size)) + 
            #     geom_text(aes(x,y, label=significance,group=NULL), data=sig_label_df, size = 10, color = "black") ## add the significance dots
            
            p_one_poster <- ggplot(df_brain,aes_string(x = 'group', y = 'PC1_scaled', colour='group_colour')) + 
                facet_grid(cell_type~plot_title,drop = TRUE,scale="free", space='free') + 
                ylab('Marker gene profiles')+
                ylim(-0.1, 1.1) +
                theme(text=element_text(family = 'Arial'),
                      #legend.position="none",
                      axis.title.x=element_blank(),
                      axis.title.y=element_text(size = poster_font_size),
                      title = element_text(size = poster_font_size, colour = 'black'),
                      axis.text.x = element_text(size = poster_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
                      axis.text.y=element_text(size = poster_font_size),
                      strip.text.x = element_text(size = poster_font_size),
                      strip.text.y = element_text(size = poster_font_size)) + 
                geom_text(aes(x,y, label=significance,group=NULL), data=sig_label_df, size = 10, color = "black") ## add the significance dots

            ## plot 1: plain box plot with points
            p_box <- p_one_plot +
                # geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                theme(legend.position="none")
                # scale_colour_manual(values = box_colour) # customize color
            #             p_box <- p_one_plot + geom_boxplot(lwd=0.3, outlier.shape = NA)+
            #                 geom_jitter(aes(colour = Study),width = 0.2, height = 0.2,size = 0.8, alpha = 0.5)
            (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box.', plot_type))
            ggsave(filename = f_plot, plot = p_box, width = 6, height =9 , units ='in')
            
            ## plot 2: plain violin plot with points
            # if(violin){
            #     p_violin <- p_one_plot + geom_violin(lwd=0.3)+ geom_boxplot(width = 0.1,lwd=0.3, outlier.shape = NA)+
            #         geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5)+
            #         theme(legend.position="none")+
            #         scale_colour_manual(values = box_colour) # customize color
            #     
            #     (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','violin.', plot_type))
            #     ggsave(filename = f_plot, plot = p_violin, width = 6, height =9 , units ='in')
            # }
            
            if(outlier_p){
                ## plot 3. box plot with outliers labelled
                p_box <- p_one_plot + geom_boxplot(lwd=0.3, outlier.shape = NA)+
                    geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                    geom_text(aes(label=outlier),hjust=0.5, vjust=0.5)+
                    theme(legend.position="none")+
                    scale_colour_manual(values = box_colour) # customize color
                (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_outlier.', plot_type))
                ggsave(filename = f_plot, plot = p_box, width = 20, height =18 , units ='in')
            }

            
            ## plot 4
            ## assign color for studies
            study <- c(levels(df_brain$Study))
            print(study)
            (palette <- c(ggColorHue(length(study)), box_colour))
            
            (names(palette) <- c(study, levels(df_brain$group_colour)))
            p <- p_one_plot +#  geom_boxplot(lwd=0.3, outlier.shape = NA)+
                geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                geom_line(stat="summary", fun.y="median", aes(group=as.character(Study), colour = Study), 
                          size =0.7,alpha = 0.5) +
                scale_colour_manual(values = palette) # customize color
            (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_indi_studies.', plot_type))
            ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
            # colour by indi studies
            p <- p_one_plot + # geom_boxplot(lwd=0.3, outlier.shape = NA)+
                geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5, aes(colour= Study)) +
                scale_colour_manual(values = palette) # customize color
            (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_color_studies.', plot_type))
            ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
            if(poster){ ## plot for poster (thicker lines)
                p <- p_one_poster + # geom_boxplot(lwd=0.8, outlier.shape = NA)+
                    geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5, aes(colour= Study)) +
                    scale_colour_manual(values = palette) # customize color
                (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_color_studies_poster.', plot_type))
                ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
                
            }
            
            
            
            
            
            ## plot 5 by mouse model type
            (palette <- c(model_palette, box_colour))
            p <- p_one_plot + # geom_boxplot(lwd=0.3, outlier.shape = NA)+
                geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                geom_line(stat="summary", fun.y="median", aes(group=as.character(Model_types), colour = Model_types), 
                          size =0.7,alpha = 0.5) +
                scale_colour_manual(values = palette) # customize color for each model type
            (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_indi_models.', plot_type))
            ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
            # color dots by mouse models
            p <- p_one_plot + # geom_boxplot(lwd=0.3, outlier.shape = NA)+
                geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5,aes(colour = Model_types)) +
                scale_colour_manual(values = palette) # customize color for each model type
            (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_color_models.', plot_type))
            ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
            
            if(poster){ ## plot for poster (thicker lines)
                p <- p_one_poster + # geom_boxplot(lwd=0.8, outlier.shape = NA)+
                    geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                    geom_line(stat="summary", fun.y="median", aes(group=as.character(Model_types), colour = Model_types), 
                              size =0.7,alpha = 0.5) +
                    scale_colour_manual(values = palette) # customize color for each model type
                (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_indi_models_poster.', plot_type))
                ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
                
                
                p_box <- p_one_poster + # geom_boxplot(lwd=0.8, outlier.shape = NA)+
                    geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                    theme(legend.position="none")+
                    scale_colour_manual(values = box_colour) # customize color
                (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_poster.', plot_type))
                ggsave(filename = f_plot, plot = p_box, width = 6, height =9 , units ='in')
                
            }
            
        } ## end loop2  for brain cell types
        
        
        
        ## here to plot for each mouse types
        print('plot each mouse models')
        
        for(mouse_type in intersect(model_order, levels(df_cell$Model_types))){ ## for loop to creat individual plots for each mouse type (wilcox test redone)
            df_brain_model <- df_cell[which(df_cell$Model_types == mouse_type), ]%>% droplevels()
            
            ## redo the wilcox test by the models type and assign significance
            df <- df_brain_model
            df_p <- NULL
            for(c_t in levels(df$cell_type)){ #loop wilcox test
                disease_estimates <- df[which(df$cell_type == c_t & df$groups == 'Disease'), 'PC1_scaled']
                wt_estimates <- df[which(df$cell_type == c_t & df$groups == 'WT'), 'PC1_scaled']
                (p_v <- as.numeric(wilcox.test(disease_estimates, wt_estimates)[3]))
                df_tmp <- data.frame(cell_type = c_t, p=p_v)
                if(is.null(df_p)){df_p <- df_tmp} else {df_p <- rbind(df_p, df_tmp)}
            }# end loop wilcox text
            df_p$adj_p <- p.adjust(as.numeric(df_p$p), method = 'BH')
            
            df <- noWarnings(left_join(df, df_p))
            df_brain_model <- assignSig(df, 'adj_p')
            
            
            for(brain_cells in levels(df_brain_model$brain_cells)){ ## loop to plot for each mouse models
                df <- df_brain_model[which(df_brain_model$brain_cells == brain_cells), ]%>% droplevels()
                print(paste0(disease, ' ', brain_cells, ': ', mouse_type))
                
                #### get the text label for significance
                #' see https://trinkerrstuff.wordpress.com/2012/09/01/add-text-annotations-to-ggplot2-faceted-plot/
                sig_label_df <- rmDup(df[, c('cell_type', 'plot_title','significance','group')])
                sig_label_df$y = 1.05
                sig_label_df$x = sig_label_df$group
                
                p_one_plot <- ggplot(df,aes_string(x = 'group', y = 'PC1_scaled', colour='group_colour')) + 
                    facet_grid(cell_type~plot_title,drop = TRUE,scale="free", space='free') + 
                    ylab('Marker gene profiles')+
                    ylim(-0.1, 1.1) +
                    theme(text=element_text(family = 'Arial'),
                          #legend.position="none",
                          axis.title.x=element_blank(),
                          axis.title.y=element_text(size = one_plot_font_size),
                          title = element_text(size = one_plot_font_size, colour = 'black'),
                          axis.text.x = element_text(size = one_plot_font_size, angle = x_angle), ## vjust 0.5 put x labels in the middle
                          axis.text.y=element_text(size = one_plot_font_size),
                          strip.text.x = element_text(size = one_plot_font_size),
                          strip.text.y = element_text(size = one_plot_font_size)) + 
                    geom_text(aes(x,y, label=significance,group=NULL), data=sig_label_df, size = 10, color = "black") ## add the significance dots
                
                box_colour <- c('dodgerblue','indianred1')
                (names(box_colour) <- c(levels(df$group_colour)))
                
                
                ## plot 5 customize color for each model type
                (palette <- c(model_palette, box_colour))
                p <- p_one_plot + geom_boxplot(lwd=0.3, outlier.shape = NA)+
                    geom_jitter(width = jitter_w, height = jitter_h,size = 0.8, alpha = 0.5) +
                    geom_line(stat="summary", fun.y="median", aes(group=as.character(Model_types), colour = Model_types), 
                              size =0.7,alpha = 0.5) +
                    scale_colour_manual(values = palette) # customize color for each model type
                (f_plot <- paste0(f_out_dir,disease,'_', brain_cells,'_','box_indi_models', gsub(' ', '_', mouse_type), '.', plot_type))
                ggsave(filename = f_plot, plot = p, width = 8, height =9 , units ='in')
                
            } ## end loop to plot for each mouse models
            
        }  ## end  loop to creat individual plots for each mouse type (wilcox test redone)
        
        
        
        
        
        
    }## end loop disease
    
    
}







