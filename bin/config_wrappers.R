#' 2016-04-25
## load the config files for wrapper

## download theunfiletered profiles of GSE1556 from https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=27
## and save in the gemma_data folder 
disease ='AD'
disease_ls <- c('AD')

home_dir <- file.path(here::here(),'files_results')


dir.create(home_dir,showWarnings = F)

disease_dir <- file.path(home_dir, 'AD_mouse_model_project/')

platform_folder_biological_pathway_only <- file.path(home_dir,'/platforms/biological_pathway_only/')

platform_folder_biological_all_GO <- file.path(home_dir, '/platforms/')

(gemma_data <- file.path(home_dir, '/AD_mouse_model_project/data_and_QC/gemma_data/unfiltered_exp_data/'))


#####
## emineJ files
xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'




######
## experimental design info for all datasets

f_design <- file.path(here::here(),'configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv')

## a gene list for sanity heatmap check
f_sanity_check_genes <-file.path(here::here(),'configs/AD_mouse_dataset_doc/gene_list_for_sanity_check.tsv')

## must have home_dir defined
## define disease main dir



# ##***********************
# ## which dataset need batch correction
# ##***********************
batch_correction_config <- c('GSE64398.1', 'GSE50521','GSE36237') 

##***********************
## define limma dir based on brain region
##***********************
brain_region <- list(AD = 'hippocampus',
                     HD = 'striatum',
                     PD = 'striatum')

limma_dir <- file.path(disease_dir, '/limma_DE/',brain_region[disease], '/')
limma_dir_FG <- file.path(disease_dir, '/limma_DE_FG/',brain_region[disease], '/')


# ##***********************
# ## get the know modifier info
# ##***********************
# switch(disease,
#        AD = eval(parse(text = paste0("known_modifiers_config <- '", home_dir, "/git/doc/AD_mouse_dataset_doc/neurogem_modifiers_mouse_AD.tsv'"))),
#        HD = eval(parse(text = paste0("known_modifiers_config <- '", home_dir, "/git/doc/HD_mouse_dataset_doc/ref_genetic_modifiers.tsv'"))),
#        #PD = eval(parse(text = "known_modifiers_config <-'/home/bzhuang/git/doc/PD_mouse_dataset_doc/ref_genetic_modifiers.tsv'")))
#        PD = eval(parse(text = "known_modifiers_config <-''")))
# 
# ##***********************
# ## list of excluded genotypes in meta analysis
# ##***********************
# switch(disease,
#        AD = eval(parse(text = "exclude_genotype_config <- c('N_dC_KO', 'APLP2_KO')")),
#        #AD = eval(parse(text = "exclude_genotype_config <- c('N_dC_KO', 'APLP2_KO', 'ho_TASTPM','TAS10','TAU','het_TPM')")),
#        HD = eval(parse(text = "exclude_genotype_config <- c('Grm5_KO', 'HdhQ111_Grm5_KO', 'R62_p62_KO', 'p62_KO')")),
#        PD = eval(parse(text = "exclude_genotype_config <- c('Hsp70', 'Ache_S', 'Ache_R', 'R_MPTP', 'hLRRK2_WT')")))
# 
# 
# if(grepl('ND_project_combined', home_dir)){
#     switch(disease,
#            AD = eval(parse(text = "exclude_genotype_config <- NULL")),
#            HD = eval(parse(text = "exclude_genotype_config <- c('Grm5_KO', 'HdhQ111_Grm5_KO', 'R62_p62_KO', 'p62_KO')")),
#            PD = eval(parse(text = "exclude_genotype_config <- c('Ache_R', 'R_MPTP', 'hLRRK2_WT')")))
#     
# }
# 
# 
# # R_MPTP only has 1 sample GSE31458
# # 'hLRRK2_WT' is human LRRK2 wildtype (GSE52584)
# 
# ##***********************
# ## which dataset need batch correction
# ##***********************
# switch(disease,
#        AD = eval(parse(text = "batch_correction_config <- c('GSE36237', 'GSE64398.1', 'GSE50521')")),
#        HD = eval(parse(text = "batch_correction_config <- c('GSE9038','GSE9803','GSE9857','GSE10202', 'GSE32417')")),
#        PD = eval(parse(text = "batch_correction_config <- c('GSE60080')"))
#        #PD = eval(parse(text = "batch_correction_config <- c('GSE24838', 'GSE20547.2', 'GSE60080')"))
#        )
# 
# # #dataset_ls <- c('GSE36237','GSE13691.1', 'GSE13691.2')  # datasets for batch correction for AD for more datasets
# 
# ##***********************
# ## for cell population estimate, the marker genes are based on tissue:
# ##***********************
# switch(disease,
#        AD = eval(parse(text = "cell_pop_marker_genes_config <- './brainCellTypeSpecificGenes-master/analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Hippocampus/'")),
#        HD = eval(parse(text = "cell_pop_marker_genes_config <- './brainCellTypeSpecificGenes-master/analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Striatum/'")),
#        PD = eval(parse(text = "cell_pop_marker_genes_config <- './brainCellTypeSpecificGenes-master/analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Striatum/'")))
