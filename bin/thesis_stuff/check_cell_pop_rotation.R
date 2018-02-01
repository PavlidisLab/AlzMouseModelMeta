## 2017-01-09
## summary of cell markers and estimation of rotations and variations
## save results in, tables and rdata
#'outdir <- paste0(home_dir, '/git/results/ND_results/DE_genes_markers/rotations/')
#'('../../results/cell_marker_summary.Rdata')
#'
#+++++++++++++++++++
## get which DE genes are also cell type markers before and after correction
#+++++++++++++++++++
## thesis tables /results/ND_results/DE_genes_markers/'

# rm(list=setdiff(ls(),'home_dir'))
#  
library(HelperFunctions)
library(dplyr)
library(pheatmap)


source('mixed_models/plot_up_down_genes_thesis_figure_corrected_value.R')

phase_ls <- c('early', 'late')
dir.create(plotdir,recursive = T,showWarnings = F)




df_all <- NULL
df_ca <- NULL
msg <- paste0('# ', Sys.Date())


for(disease in disease_ls){## loop disease
    print(disease)
    msg <- paste0(msg, "\n#####\n# DISEASE: ", disease,"\n#####")
    (f_dir <- paste0(disease_dir, '/MGP_estimation/')) ## results of MGP estimation as input
    
    
    for(phase in phase_ls){ ## loop phase
        
        
        ## define heatmap input
        ## plot heatmaps
        legend_ls <- list(c('Genotype', 'Study'))
        regulation_ls <- c('up', 'down')
        mm_rdata_keyword <- c('_include_NA_low_exp_rm')  ## for only study corrected matrix
        model <- c('random_intercept')
        model_keyword <- c('_include_NA_low_exp_rm')  
        result_rank <- c('mixed_model_jackknife')  ## which rank list to get genes from
        opposite_phase = F
        top_genes = F
        rdata_keyword = 'mixed_model_results_exp_corrected.Rdata'  # choose which rdata to get from
        ## where the mixed model result files
        mm_rdata_dir <- paste0(disease_dir, '/mixed_model/', model,mm_rdata_keyword,'/') ## where the expression results
        mm_dir <- paste0(disease_dir,'/',result_rank,'/', model,mm_rdata_keyword,'/') ## where the ranked results
        
        
        ## input cell marker files
        print(phase)
        (file_dir <- max(list.dirs(paste0(f_dir, '/',phase))))
        cat(paste0('\n Input marker gene profiles, rotTables: ', file_dir, '\n'))
        ## list of rotation tables
        f_ls <- list.files(file_dir,pattern = 'rotTable.tsv')
        (f_ls <- grep('activation', f_ls, value = T, invert = T))
        
        msg <- paste0(msg, "\n# PHASE: ", phase, "\n# File dir: ", file_dir)
        for(i in 1:length(f_ls)){ ## for each cell type
            (f <- f_ls[i])
            print(f)
            (cell_t <- gsub(' rotTable.tsv', '', f))
            f <- paste0(file_dir, '/',f)
            ## read the first line (variation explained)
            (var_explained <- gsub('# Variation explained: ', '', readLines(f, 1)))
            var_explained <- as.numeric(unlist(strsplit(var_explained, split = ' ')))
            df <- read.delim(f, comment.char = '#', header = F)
            colnames(df) <- c('geneSymbol', 'rotation')
            
            ## get the markers and rotaion and variation explained
            if(length(var_explained) == nrow(df)){
                df <- cbind(cell_type = cell_t, df, var_explained=var_explained,disease=disease, phase=phase, marker_order =1:nrow(df)) 
            }else{
                df <- cbind(cell_type = cell_t, df, var_explained = NA,disease=disease, phase=phase, marker_order =1:nrow(df))
            }
            
            
            ##get the grp rotations for wt and disease
            f_g <- paste0(file_dir, '/',cell_t,' groupRots')
            df_g <- read.delim(f_g)
            df_g$geneSymbol <- rownames(df_g)
            df_g$D_WT_rot_diff <- df_g$Disease-df_g$WT
            
            df <- noWarnings(left_join(df, df_g))
            
            if(is.null(df_all)){
                df_all <- df
            }else{
                df_all <- rbind(df_all, df)
            }
            
            ## get the summary count
            df_count <- as.data.frame(table(df[, c(1, 5,6)]))
            
            df_count <- cbind(df_count, 
                              pos_rotation = length(which(df$rotation >=0)),
                              neg_rotation = length(which(df$rotation <0)),
                              pos_diff = length(which(df$D_WT_rot_diff >=0)),
                              neg_diff = length(which(df$D_WT_rot_diff <0)))
            
            if(is.null(df_ca)){
                df_ca <- df_count
            }else{
                df_ca <- rbind(df_ca, df_count)
            }
            
            ## plot heatmap
            gene_list <- as.character(df$geneSymbol)
            
            plotThesisHeatGenes(mm_dir,mm_rdata_dir,
                                plt_out=plotdir,
                                phase_ls = phase, 
                                threshold,
                                rdata_keyword = rdata_keyword,
                                rm_gene_list = rm_gene_list,
                                top_genes =top_genes, 
                                opposite_phase =opposite_phase, 
                                gene_list =gene_list,
                                row_width =25,
                                legend_ls =legend_ls,
                                plot_indi = F,
                                plot_suffix = cell_t,
                                plotdir=T)

            plotThesisHeatGenes(mm_dir,mm_rdata_dir,
                                plt_out=paste0(plotdir, '/clustered/'),
                                phase_ls = phase, 
                                threshold,
                                rdata_keyword = rdata_keyword,
                                rm_gene_list = rm_gene_list,
                                top_genes =top_genes, 
                                opposite_phase =opposite_phase, 
                                gene_list =gene_list,
                                row_width =25,
                                legend_ls =legend_ls,
                                plot_indi = F,
                                plot_suffix = paste0(cell_t, '_clustered'),
                                plotdir=T,
                                cluster_rows =T)
            
            
        }## loop cell type
        
    }## loop phase
    
}## loop disease






## write the tables

writeTable(df_all, f_out = paste0(outdir, 'Rotations_', Sys.Date(), '.tsv'), msg = msg)
writeTable(df_ca, f_out = paste0(outdir, 'marker_counts_', Sys.Date(), '.tsv'), msg = msg)

## save r object
cell_marker_rotations <- df_all
cell_marker_freq <- df_ca
cell_marker_msg <- msg
save(cell_marker_freq, cell_marker_rotations, cell_marker_msg, 
     file = f_rdata_out )

cat(paste0('\n#############\nR object out: ', f_out, '\n Tables out: ', outdir, '\n Plots : ', plotdir, 
           '\n need to run getCellMarkers(disease, phase, threshold) and load all_rank, all_ranks_adj to 
           create thesis_DE_markers.tsv file'))


