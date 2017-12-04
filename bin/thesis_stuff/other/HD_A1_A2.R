# df <- all_ranks
# 
# up <- df[which(df$disease == 'HD' & df$cell_type =='Astrocyte' & df$Estimate >0 & df$P_adj <0.1), ]
# down <- df[which(df$disease == 'HD' & df$cell_type =='Astrocyte' & df$Estimate <0& df$P_adj <0.1), ]
# 
# up <- df[which(df$disease == 'HD' & df$cell_type =='Astrocyte' & df$Estimate >0), ]
# down <- df[which(df$disease == 'HD' & df$cell_type =='Astrocyte' & df$Estimate <0), ]
# 
# 
# table(up$phase)
# table(down$phase)

#' 2017-01-31
## looking at HD astrocyte marker rotation estimations
## and the logFC of these markers in LPS astrocytes (supposedly corresponding to A1 astrocytes)
## if the markers (pos and neg rotation) can seperate the A1 from WT (measured by logFC)
#' the seperation should be clear between pos and neg -- than the A1 change may affect HD astrocye estimation

f='C:/Users/User/Documents/lab_results/ND_project_combined/ND_results/cell_population/all_sample_estimation/PD/early/OrganismPart_Astrocyte_LPS/results/GSE35338_1_day_limma_toptable_GenotypeLPS-vs-baseline_WT_.tsv'
df <- read.delim(f,comment.char = '#')
f='C:/Users/User/Documents/lab_results/ND_project_combined/ND_results/cell_population/all_sample_estimation_test/HD/late/2017-01-31/Astrocyte rotTable.tsv'
hd <- read.delim(f,comment.char = '#',header = FALSE)
colnames(hd) <- c('GeneSymbols', 'rotation')
astro <- left_join(hd, df)
astro$abs_FC <- abs(astro$LogFC)

astro <- dplyr::arrange(astro, dplyr::desc(abs_FC))
astro <- astro[!duplicated(astro$GeneSymbols), ]


# f='C:/Users/User/Documents/lab_results/ND_project_combined/ND_results/cell_population/all_sample_estimation/PD/early/OrganismPart_Astrocyte_LPS/results/'
# df <- read.delim(f,comment.char = '#')
# table(astro$rotation_dir)

astro$rotation_dir =''
astro$rotation_dir[which(astro$rotation >0)] ='positive'
astro$rotation_dir[which(astro$rotation <=0)] ='negative'
astro$rotation_dir <- as.factor(astro$rotation_dir)
p_one_plot <- ggplot(astro,aes_string(y = 'LogFC', x = 'rotation_dir', colour='rotation_dir'))
    p_one_plot +
        geom_violin(lwd=0.3) +
        #geom_boxplot(lwd=0.3, outlier.shape = NA)+
    geom_jitter(width = 0.1, height = 0.1,size = 1.5, alpha = 0.5) +
    theme(legend.position="none")

    
    