###### 2017-01-31

## looking at markers in GSE35338 (LPS and MCAO of reactive astrocytes)
## use striatum markers and A1 A2 markers from paper Liddelow et al 2017
## copy limma files to '/ND_results/cell_population/all_sample_estimation/PD/early/'

rm(list=setdiff(ls(),'home_dir'))
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
source("estimate_cell_population.R")

#for A1 and A2 markers
##__
## get A1 A2 markers
df <- read.delim('../../doc/A1_A2_non_activated_astrocyte_markers.tsv', comment.char = '#')
df <- filterContain(df, column ='regulation', value = 'up')

#setdiff(df$geneSymbol, all_ranks$geneSymbol)

mouseMarkerGenes$Striatum$Astrocyte= setdiff(mouseMarkerGenes$Striatum$Astrocyte, c("Slc1a2","Slc1a3"))
print(paste0('removed "Slc1a2","Slc1a3" from markers'))

(A1_markers <- as.character(df$geneSymbol[which(df$cell_type =='A1_reactive')]))
(A2_markers <- as.character(df$geneSymbol[which(df$cell_type =='A2')]))

intersect(mouseMarkerGenes$Striatum$Astrocyte, union(A1_markers, A2_markers))
intersect(mouseMarkerGenes$Striatum$Microglia_activation, A1_markers)
intersect(mouseMarkerGenes$Striatum$Microglia_deactivation, A2_markers)


intersect(mouseMarkerGenes$Hippocampus$Astrocyte, union(A1_markers, A2_markers))
intersect(mouseMarkerGenes$Hippocampus$Microglia_activation, A1_markers)
intersect(mouseMarkerGenes$Hippocampus$Microglia_deactivation, A2_markers)


mouseMarkerGenes$Striatum$A1 = setdiff(A1_markers, A2_markers)
mouseMarkerGenes$Striatum$A2 = setdiff(A2_markers, A1_markers)
mouseMarkerGenes$Hippocampus$A1 = setdiff(A1_markers, A2_markers)
mouseMarkerGenes$Hippocampus$A2 = setdiff(A2_markers, A1_markers)

##----



disease_ls =c('PD')
phase_ls = c('early')


# 
# disease_ls =c('HD')
# phase_ls = c('early', 'late')

all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Genotype'  # design_group = 'Original_genotype'
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT

threshold = -10 # filter probes exp less than threshold (input is study corrected value, and low expression filtered)
mixed_model =F
removeNegatives = F  ## ogan suggested use F, 


for (disease in disease_ls){
    if(disease == 'AD'){
        genes =mouseMarkerGenes$Hippocampus
    }else{
        genes =    mouseMarkerGenes$Striatum 
    } 
    
    
    for(phase in phase_ls){
        source('config/config_wrappers.R')
        (file_ls_all <- paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/',disease, '/', phase,'/Timepoint_1_day'))
        
        # file_ls_all <- list.dirs(file_ls, recursive = F, 'Timepoint_1_day')
        
        for(file_ls in file_ls_all){
            (r_ob <- list.files(file_ls, full.names=T, recursive=T, pattern= "limma_objects.Rdata"))
            (outDir  <- file_ls)
            dir.create(outDir,showWarnings=F, recursive=T)
            
            cellPopEstimatePlot(r_ob, outDir, design_group = design_group, mixed_model = mixed_model,
                                sigTest = sigTest, 
                                wt_only = wt_only,
                                removeNegatives = removeNegatives, threshold = threshold, genes = genes)
        }
        }

}
