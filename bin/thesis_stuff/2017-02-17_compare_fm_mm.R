# 2017-02-17
home_dir_ls <- c('C:/Users/User/Documents/lab_results/ND_project_combined/',
                 '/home/bzhuang/ND_project_combined/')

home_dir <- home_dir_ls[which(dir.exists(home_dir_ls))]
rm(list=setdiff(ls(),'home_dir'))
 
library(HelperFunctions)
library(dplyr)
library(ggplot2)

 

load('../../results/fixed_model_ranks_2017-02-17.Rdata')
load('../../results/ranks_paths_2017-02-02.Rdata')

######
# check correlation
######


corFunction <- function(df, df2, disease, phase){
    print(paste0(disease, ': ', phase))
    df <- df[which(df$disease == disease & df$phase == phase), ]
    df2 <- df2[which(df2$disease == disease & df2$phase == phase), ]
    df2 <- df2[, c('geneSymbol', 'up_jack')]
    colnames(df2) <- c('geneSymbol', 'fm_up_jack')
    df.a <-  noWarnings(na.omit(left_join(df[, c('geneSymbol', 'up_jack')],df2[, c('geneSymbol', 'fm_up_jack')])))
    
    x <- cor.test(df.a[, 2],df.a[, 3], method = 'spearman', exact = F)
    print(paste0('cor: ', x$estimate, ', pvalue: ', x$p.value))
}


disease_ls = c('AD', 'HD')
phase_ls = c('early', 'late')

for(disease in disease_ls){
    for(phase in phase_ls){
        corFunction(all_ranks,fm_all_ranks, disease, phase)
    }
}

for(disease in disease_ls){
    for(phase in phase_ls){
        corFunction(all_ranks_adj,fm_all_ranks_adj, disease, phase)
    }
}



for(disease in disease_ls){
    for(phase in phase_ls){
        corFunction(all_ranks,all_ranks_adj, disease, phase)
    }
}
