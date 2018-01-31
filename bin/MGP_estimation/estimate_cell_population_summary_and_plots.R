#2016-09-17
#' updated: 2017-02-20: made into wrappers
##
#' after all the estimate cell populations
#' get a summary table for all and plots, saved in 'ND_results/cell_population/'
# disease_ls = c('AD', 'HD', 'PD')
# phase_ls =c('early', 'late')
# source('estimate_cell_population_summary_and_plots.R')

library('HelperFunctions')
library(dplyr)
library(ggplot2)


###########################
## helper functions
###########################



addSpace <- function(x, maxchar){
    x <- paste0(paste0(rep(' ', maxchar-nchar(x)), collapse=''), x)
    return(x)
}






addCellType <- function(df_disease_all){
    ###
    ## add cell type category, and reorder cell types (and rename)
    ###
    full_name <- c('GabaSSTReln', 'StriatumCholin','Pyramidal_Thy1', 'ForebrainCholin', 'Oligo', 'Spiny', 'Astrocyte')
    new_name <- c('GABAergic','Cholinergic', 'Pyramidal','Cholinergic', 'Oligodendrocytes','MSNs','Astrocytes')
    
    
    df_disease_all[, 'cell_type'] <- plyr::mapvalues(df_disease_all[, 'cell_type'], from = full_name, 
                                                     to = new_name, warn_missing = F)
    full_name <- c("Astrocytes", "Microglia",'Oligodendrocytes', 
                   "DentateGranule","GABAergic", "Pyramidal","MSNs","Cholinergic",
                   'Microglia_deactivation', 'Microglia_activation')
    new_name <- c(rep('glia', 3), rep('neurons', 5), rep('microglia_act', 2))
    
    df_disease_all$brain_cells <- df_disease_all[, 'cell_type']
    df_disease_all[, 'brain_cells'] <- plyr::mapvalues(df_disease_all[, 'brain_cells'], from = full_name, 
                                                       to = new_name, warn_missing = F)
    return(df_disease_all)
}








processCellPopDf <- function(disease, disease_dir,df_info, original_genotype =F,
                             write_df =F,
                             f_out_dir = NULL){
    #' get all the cell estimate for the disease from each individual studies (estimation is made by each study
    #' separately), add test significance
    #' df_info: e.g. dataset_info_all.tsv to add the disease phase info for each sample (by study and timepoint)
    #' 
    ## get the genotypes (including combined, but will be replaced later from the orginal genotype files)
    folder = paste0(disease_dir, '/MGP_estimation/')
    (folder = max(list.dirs(folder, recursive=F)))
    (f_ls = grep('_cell_proportion_estimation_scaled.tsv',list.files(folder, recursive=T), value=T))
    (f_estimate_ls = grep('estimatefile.tsv',list.files(folder, recursive=T), value=T))
    (keyword_ls = unlist(lapply(f_ls, function(x) unlist(strsplit(x, split = '/'))[1])))
    (f_ls = paste0(folder, '/',f_ls))
    (f_estimate_ls =paste0(folder, '/',f_estimate_ls))
    
    if(original_genotype){
        ## get from the orginal genotype (multiple disease genotypes) and add to the f_estimate_ls
        folder = paste0(disease_dir, '/results/Cell_population_estimates_Original_genotype/')
        (folder = max(list.dirs(folder, recursive=F)))
        (f_ls_g = noWarnings(grep('_cell_proportion_estimation_scaled.tsv',list.files(folder, recursive=T), value=T)))
        (f_estimate_ls_g = noWarnings(grep('estimatefile.tsv',list.files(folder, recursive=T), value=T)))
        if(length(f_ls_g) >0){
            (keyword_ls_g = unlist(lapply(f_ls_g, function(x) unlist(strsplit(x, split = '/'))[1])))
            (f_ls_g = paste0(folder, '/',f_ls_g))
            (f_estimate_ls_g = paste0(folder, '/',f_estimate_ls_g))
            
            index = which(keyword_ls %in% keyword_ls_g)
            print(paste0('replacing ', paste0(keyword_ls[index], collapse=',')))
            keyword_ls[index] = keyword_ls_g
            f_ls[index] = f_ls_g
            f_estimate_ls[index] = f_estimate_ls_g
        }
    }
    
    df_all=NULL
    df_disease_all = NULL
    ## for group/study in the disease
    for(i in 1:length(f_ls)){
        (f= f_ls[i])
        print(keyword_ls[i])
        (f_e= f_estimate_ls[i])
        # print(i)
        # print(f)
        # print(f_e)
        
        ## get the samples names
        df_e <- read.delim(f_e, comment.char = '#')
        df_e$Sample = rownames(df_e)
        df_e <- reshape2::melt(df_e, id.vars = c('Sample', 'groups'), variable.name="cell_type", value.name = 'PC1')
        # sub the .1 .2
        df_e$Sample <- gsub('GSE64398\\.1','GSE64398', df_e$Sample )
        df_e$Sample <- gsub('GSE63617\\.1|GSE63617\\.2','GSE63617', df_e$Sample )
        
        df = read.delim(f, comment= '#')
        df = noWarnings(left_join(df, df_e))
        
        df$keyword = keyword_ls[i]
        if(is.null(df_all)){
            df_all = df
        }else{
            df_all =rbind(df_all, df)
        }
    }
    df_all$disease = disease
    
    ## add the phase info
    
    print(paste0('Input dataset info: ', df_info))
    info_df <- read.delim(df_info, comment.char='#')
    info_df$group <- info_df$Genotype
    info_df$keyword <- paste0(info_df$Dataset, '_', info_df$Timepoint)
    ## also duplicate for WT
    info_df2 <- info_df
    info_df2$group <- 'WT'
    info_df <- rbind(info_df, info_df2)
    
    df_disease = noWarnings(left_join(df_all, info_df[, c(intersect(colnames(info_df), colnames(df_all)), 'Phase')]))
    
    if(is.null(df_disease_all)){
        df_disease_all = df_disease
    }else{
        df_disease_all =rbind(df_disease_all, df_disease)
    }
    df_disease_all <- rmDup(df_disease_all) ## gseGSE64398.1_4_months WT will be duplicated for early and late
    
    
    ## add the plot title
    df_disease_all$plot_title <- gsub('_months', 'm)', df_disease_all$keyword)
    df_disease_all$plot_title <- gsub('_weeks', 'w)', df_disease_all$plot_title)
    df_disease_all$plot_title <- gsub('\\..*[1|2|3]_', '(', df_disease_all$plot_title)
    df_disease_all$plot_title <- gsub('_', '(', df_disease_all$plot_title)
    
    ## assign the significance
    #' <0.01: ***
    #' <0.05 : **
    #' <0.1:*
    df_disease_all$significance =''
    df_disease_all$significance[which(df_disease_all$adj_p_wt<=0.1)] = '*'
    df_disease_all$significance[which(df_disease_all$adj_p_wt<=0.05)] = '**'
    df_disease_all$significance[which(df_disease_all$adj_p_wt<=0.01)] = '***'
    df_disease_all$significance <- as.factor(df_disease_all$significance)
    
    
    msg <- paste0('#Plot date: ', Sys.Date(), '\n#Input folder: ', folder, '\n#Input dataset info: ', df_info)
    
    
    ## add cell types to the table
    df_disease_all <- addCellType(df_disease_all)
    
    ## save disease all table
    
    if(write_df){
        dir.create(f_out_dir, recursive=T, showWarnings=F)
        f=paste0(f_out_dir, disease, '_cell_population_estimate_', Sys.Date(), '.tsv')
        writeTable(df_disease_all, f_out=f,msg = msg)
        print(paste0('output table: ', f))
    }
    
    ## prepare plot titles
    df_disease_all$plot_title <- gsub('\\(', '\n\\(', df_disease_all$plot_title)  ## make the timepoint in the second line
    
    
    
    return(list(df_disease_all, msg))
}





#####
#plot functions
#####

plotIndiStudyCellEstimates <- function(disease, phase, df_disease_all, f_out_dir, plot_indi=F,
                                       font_size=30,
                                       x_angle=90){
    
    #' f_out_dir: plot out dir all studies in the same disease phase
    print(paste0('############ Plot ',disease, ' ' , phase))
    

    
    ## get the disease, and phase
    index <- Reduce(intersect, list(which(df_disease_all$disease == disease),
                                    which(df_disease_all$Phase== phase)))
    df_cell <- df_disease_all[index, ]%>%droplevels()
    
    
    #####
    ## make all genotypes the same length(add space before the names)
    #####
    ## replace the long names
    new_name <- c('TPM','Pink1_KO\nSNCA_A53T','R6/2', 'HdhQ111/Q111','HdhQ92/Q92','CHL2Q150/Q150', 'Lrrk2_KO',
                  'MPTP(AchE-S)', 'App_KO', 'Aplp2_KO','N171(98Q)','R6/2_300Q','LRRK2_R1441G','NdC_KO')
    full_name <- c('het_TPM', 'PINK1_KO_SNCA_A53T', 'R62','HdhQ111','HdhQ92_Q92', 'CHL2Q150_Q150', 'LRRK2_KO',
                   'S_MPTP','APP_KO','APLP2_KO', 'N171_98Q', 'R62_300Q', 'hLRRK2_R1441G','N_dC_KO')
    
    df_cell[, 'group'] <- plyr::mapvalues(df_cell[, 'group'], from = full_name, 
                                          to = new_name, warn_missing = F)
    
    ## relevel
    new_levels <- c('WT', setdiff(unique(df_cell$group), 'WT'))
    df_cell$group <- factor(df_cell$group, levels = new_levels)
    ## WT has 1 colour, others have the other
    df_cell$group_colour <- 'WT'
    df_cell$group_colour[which(df_cell$group !='WT')] <- 'model'
    
    
    #####
    ## for thesis, plot all studies, and all cell types in glia and neurons
    #####
    ## make sure 
    df_cell$brain_cells <- as.factor(df_cell$brain_cells)
    
    for(brain_cells in levels(df_cell$brain_cells)){
        df_brain <- df_cell[which(df_cell$brain_cells == brain_cells), ]%>% droplevels()
        print(paste0(disease, ' ', phase, ' ', brain_cells))
        
        one_plot_font_size = 9
        
        
        #### get the text label for significance
        #' see https://trinkerrstuff.wordpress.com/2012/09/01/add-text-annotations-to-ggplot2-faceted-plot/
        sig_label_df <- rmDup(df_brain[, c('cell_type', 'plot_title','significance','group')])
        sig_label_df$y = 1.05
        sig_label_df$x = sig_label_df$group
        
        p_one_plot <- ggplot(df_brain,aes_string(x = 'group', y = 'PC1_scaled', colour='group_colour')) + 
            facet_grid(cell_type~plot_title,drop = TRUE,scale="free", space='free') + 
            #geom_boxplot(lwd=0.3, outlier.size = 0.9, outlier.shape = 2) +
            geom_boxplot(lwd=0.3, outlier.shape = NA) +
            theme_bw() + 
            geom_point(aes_string(colour='group_colour')) +
            geom_jitter(width = 0.2, height = 0.2,size = 0.8, alpha = 0.7)+
            ylab('Marker gene profiles')+
            ylim(-0.1, 1.1) +
            theme(text=element_text(family = 'Arial'),
                  legend.position="none",
                  #axis.text.y=element_blank(), ## the number on the y
                  #axis.ticks.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_text(size = one_plot_font_size),
                  title = element_text(size = one_plot_font_size, colour = 'black'),
                  axis.text.x = element_text(size = one_plot_font_size, angle = x_angle, hjust = 1, vjust=0.5), ## vjust 0.5 put x labels in the middle
                  axis.text.y=element_text(size = one_plot_font_size),
                  strip.text.x = element_text(size = one_plot_font_size),
                  strip.text.y = element_text(size = one_plot_font_size)) + 
            geom_text(aes(x,y, label=significance,group=NULL), data=sig_label_df, size = 6, color = "black") ## add the significance dots
        
        (f_plot <- paste0(f_out_dir,disease,'_',phase, '_', brain_cells,'.png'))
        ggsave(filename = f_plot, plot = p_one_plot, width = 9, height =6 , units ='in')
        print(f_plot)
        
    }
    
    
    
    
    
    if(plot_indi){
        (plot_indi_out = paste0(f_out_dir, '/', disease,'/',phase, '/'))
        dir.create(plot_indi_out, recursive = T, showWarnings = F)
        #####
        ## get the cell type, plot for individual studies (show the significant stars)
        #####
        for (cell_type in levels(df_cell$cell_type)){ # each cell type
            index <- which(df_cell$cell_type == cell_type)
            df <- df_cell[index, ]%>%droplevels()
            
            print(paste0('Plot ', disease, ', ', phase, ', ', cell_type))
            
            
            ## order the dataset by p value
            df_order <- orderCol(na.omit(rmDup(df[, c('group','keyword','adj_p_wt')])), cols = 'adj_p_wt')
            (g_order <- as.character(unique(df_order$keyword)))
            
            for (k in 1:length(g_order)){
                print(k)
                df_t <- filterContain(df, 'keyword', g_order[k])
                df_t$adj_p_wt <- round(df_t$adj_p_wt, digits = 3)
                df_t$adj_p_wt[which(is.na(df_t$adj_p_wt))] =''
                gtitle <- df_t$plot_title[1]
                print(gtitle)
                (f_plot <- paste0(plot_indi_out, cell_type,'_', k, '_',  df_t$keyword[1], '.png'))
                
                ## make sure the levels matched
                (factor_level <- c(setdiff(unique(df_t$group), 'WT'), 'WT'))
                df_t$group <- factor(df_t$group, levels=factor_level)
                label_df <- rmDup(df_t[, c('group', 'significance', 'adj_p_wt')])
                label_df <- noWarnings(left_join(data.frame(group = factor_level), label_df))
                
                
                
                #plt_ratio=14/nrow(label_df)
                
                # box of specified column, filled by group
                p <- ggplot(df_t, aes_string(x = 'group', y = 'PC1_scaled' ))+ 
                    geom_boxplot() +
                    theme_bw() + geom_point() +  ggtitle(gtitle) +
                    ylim(-0.1, 1.1) +
                    ## annotation at y=1.05 position
                    ## annotation will make facet not working
                    annotate("text", x=seq(1:nrow(label_df)), y= 1.05, 
                             label= as.character(label_df$significance),size =20,
                             colour ='red',angle = 90) +
                    #coord_fixed(ratio = plt_ratio) +
                    theme(text=element_text(family = 'Arial'),
                          legend.position="none",
                          #axis.text.y=element_blank(), ## the number on the y
                          #axis.ticks.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          title = element_text(size = font_size, colour = 'black'),
                          axis.text.x = element_text(size = font_size, angle = x_angle, hjust = 1,vjust=0.5), ## vjust 0.5 put x labels in the middle
                          axis.text.y=element_blank())## dont show number on the y)
                
                
                ##show number on the y)
                #         if(k!=1){
                #           p = p+theme(axis.text.y=element_blank())
                #             p_w= nrow(label_df) * 3 +0.4}  ## in cm ## for 2}
                
                ####
                ## plot size without y lable
                ####
                p_w= nrow(label_df) * 3  ## in cm
                p_h = 15*1.7
                
                ## save each indi plot
                ggsave(filename = f_plot, plot = p, width =p_w , height =p_h , units ='cm')
            }
        } # end each cell type
    }## end plot indi 
    

}

###########################
## mainFunction
###########################
mainPlotIndiStudyCellEstimates <- function(disease, phase_ls, disease_dir,df_info, original_genotype =F,
                                           write_df =T,
                                           f_out_dir = f_out_dir,
                                           plot_indi=F,
                                           font_size=30,
                                           x_angle=90){
    
    # to process data frame
    x <- processCellPopDf(disease, disease_dir,df_info, original_genotype =original_genotype,
                          write_df =write_df,
                          f_out_dir = f_out_dir)
    df_disease_all <- x[[1]]
    msg <- x[[2]]
    
    
    ## to plot
    for(phase in phase_ls){

        plotIndiStudyCellEstimates(disease, phase, df_disease_all, f_out_dir, plot_indi=plot_indi,
                                   font_size=font_size,
                                   x_angle=x_angle)
    }
    
    return(df_disease_all)
}
