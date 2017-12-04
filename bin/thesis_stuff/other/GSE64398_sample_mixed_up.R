home_dir="/home/bzhuang/ND_project_combined/"
home_dir <- 'C:/Users/User/Documents/lab_results/ND_project_combined/'

rm(list=setdiff(ls(),'home_dir'))

setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source('helper_functions.R')


# ## run function:
# 
# ## with all 3 tissues
# robject <- paste0(home_dir,"AD_mouse_model_project/data_and_QC/GEO_data/explore_data/explore_mode/GSE64398/no_subset/results/GSE64398_objects.Rdata")
# load(robject)
# 
# p_out <- paste0('../../results/ND_results/quality_check/GSE64398_sample_mixed_up_3tissue.png')
# plot_PCA_GSE64398(array_dat, array_design, p_out)
# 
# 
# ## with only cortex and hippocampus
# robject <- paste0(home_dir,"/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/explore_mode/GSE64398/subset/OrganismPart_Hippocampus,Cortex/results/GSE64398_objects.Rdata")
# load(robject)
# p_out <- paste0('../../results/ND_results/quality_check/GSE64398_sample_mixed_up_2tissue.png')
# plot_PCA_GSE64398(array_dat, array_design, p_out)




plot_PCA_GSE64398 <- function(array_dat, array_design, p_out){
    
    ###############################
    ###### PCAs
    ###############################
    #prcomp() use the singular value decomposition (SVD).
    ## remove constant/zero column to unit variance
    row_variance <- apply(array_dat, 1, var)
    index <- which(row_variance ==0)
    if(length(index) >0){
        df <- array_dat[-index, ]
    }else{
        df <- na.omit(array_dat)
    }
    
    
    array_pca <- prcomp(t(na.omit(df)), scale = T)  ## must have a matrix with probes as columns, and samples as rows
    #quick plot for each loadings # normally PC1 should not > 0.8, a high PC1 may indicate an strong outliers (which )
    # seperate itself from the rest of the group
    proportion_of_variance <- summary(array_pca)$importance[2, ]
    #png(filename= paste0(plot_pre, 'PCA.png'), width = width, height = height)
    #barplot(proportion_of_variance, main = dataset, ylab = "proportion of variance")
    #dev.off()
    
    
    (nobs.factor <- sqrt(nrow(array_pca$x) - 1))
    (d <- array_pca$sdev)
    (u <- sweep(array_pca$x, 2, 1 / (d * nobs.factor), FUN = '*'))
    (v <- array_pca$rotation)
    
    pc_n = 1:4      # choose the first few pcs (must <= the number of samples)
    obs.scale = 1
    df.u <- as.data.frame(sweep(u[, pc_n], 2, d[pc_n]^obs.scale, FUN='*')) # datafame of the first 2 PCs for the samples
    
    #append rotations to the array_design
    rownames(array_design) <- array_design$Sample
    array_prin_comp <- mapBindAllColumns(array_design, df.u)
    # set the color group (Genotype as default)
    (color_col <- intersect(x_col, colnames(array_prin_comp)))
    
    
    # keep_names=c('TAU_18m_Hipp_rep1','TAU_18m_Hipp_rep2','TAU_2m_Hipp_rep1','TAU_18m_Hipp_rep3',
    #              'TPM_2m_Hipp_rep1','TPM_2m_Hipp_rep2','TPM_2m_Hipp_rep3','TPM_2m_Hipp_rep4')
    
    
    ## this is the sample mixed up
    keep_names=c('TPM_2m_Hipp_rep4')
    
    array_prin_comp$Sample_new <- as.character(array_prin_comp$Sample)
    n_index <- setdiff(1:nrow(array_prin_comp), which(array_prin_comp$Sample_new %in% keep_names))
    array_prin_comp$Sample_new[n_index] <- ''
    
    
    
    # relevel the tissues
    array_prin_comp$OrganismPart <- relevel(array_prin_comp$OrganismPart, ref = 'Hippocampus')
    colnames(array_prin_comp)[which(colnames(array_prin_comp) == 'OrganismPart')]='Tissue'
    
    i="Tissue"
    
    
    ## only PC1 and PC2
    ggplot(array_prin_comp, aes_string("PC2","PC1", label = "Sample_new", color = i, hjust = 1.1)) +
        geom_point() +
        theme_bw() +
        xlab(paste0("PC2 (", proportion_of_variance[2] *100, "%)"))+  #add proportion of variance for PC2
        ylab(paste0("PC1 (", proportion_of_variance[1] *100, "%)"))+ #add proportion of variance for PC1
        geom_text(size=3)
    #ggtitle(paste0("First two principal components: ", dataset, prefix)) 
    
    ggsave(filename = p_out, width = 8, height=5)
}



## run function:

## with all 3 tissues
robject <- paste0(home_dir,"AD_mouse_model_project/data_and_QC/GEO_data/explore_data/explore_mode/GSE64398/no_subset/results/GSE64398_objects.Rdata")
load(robject)

p_out <- paste0('../../results/ND_results/quality_check/GSE64398_sample_mixed_up_3tissue.png')
plot_PCA_GSE64398(array_dat, array_design, p_out)


## with only cortex and hippocampus
robject <- paste0(home_dir,"/AD_mouse_model_project/data_and_QC/GEO_data/explore_data/explore_mode/GSE64398/subset/OrganismPart_Hippocampus,Cortex/results/GSE64398_objects.Rdata")
load(robject)
p_out <- paste0('../../results/ND_results/quality_check/GSE64398_sample_mixed_up_2tissue.png')
plot_PCA_GSE64398(array_dat, array_design, p_out)