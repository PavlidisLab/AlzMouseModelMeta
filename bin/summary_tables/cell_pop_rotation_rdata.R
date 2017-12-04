#2016-010-27
##
#' after all the estimate cell populations
#' get a summary table for all and plots, saved in 'ND_results/cell_population/'
# disease_ls = c('AD', 'HD')
# phase_ls =c('early', 'late')
# source('estimate_cell_population_summary_and_plots.R')

library('HelperFunctions')
library(dplyr)
library(ggplot2)
print('NOT IN USE- 2017-02-02')


for(disease in disease_ls){ #loop 0 for each disease
    source('config/config_wrappers.R')
    
    #+++++++++++++++++++++++
    # get all the rotation tables
    #+++++++++++++++++++++++
    ## get the genotypes (including combined, but will be replaced later from the orginal genotype files)
    folder = paste0(disease_dir, '/results/Cell_population_estimates_Genotype/')
    (folder = max(list.dirs(folder, recursive=F)))

    (f_ls = grep('rotTable.tsv',list.files(folder, recursive=T), value=T))
    (keyword_ls = unlist(lapply(f_ls, function(x) unlist(strsplit(x, split = '/'))[1])))
    (f_ls = paste0(folder, '/',f_ls))

    
    ## get from the orginal genotype (multiple disease genotypes)
    folder = paste0(disease_dir, '/results/Cell_population_estimates_Original_genotype/')
    (folder = max(list.dirs(folder, recursive=F)))
    (f_ls_g = noWarnings(grep('rotTable.tsv',list.files(folder, recursive=T), value=T)))
    
    if(length(f_ls_g) >0){
        (keyword_ls_g = unlist(lapply(f_ls_g, function(x) unlist(strsplit(x, split = '/'))[1])))
        (f_ls_g = paste0(folder, '/',f_ls_g))
        
        (index = which(keyword_ls %in% keyword_ls_g))
        print(paste0('replacing ', paste0(keyword_ls[index], collapse=',')))
        keyword_ls[index] = keyword_ls_g
        f_ls[index] = f_ls_g
    }
    
    
    
    ## get the phase info
    df_info <- paste0(disease_dir, 'config_files/dataset_info_all.tsv') 
    info_df <- read.delim(df_info, comment.char='#')
    info_df$group <- info_df$Genotype
    info_df$keyword <- paste0(info_df$Dataset, '_', info_df$Timepoint)
    ## also duplicate for WT
    info_df2 <- info_df
    info_df2$group <- 'WT'
    info_df <- rbind(info_df, info_df2)
    

    
    for(phase in phase_ls){     ## loop1 for each phase
        (phase_keyword <- unique(as.character(info_df$keyword[which(info_df$Phase == phase)])))
        
        (phase_index <- which(keyword_ls %in%phase_keyword))
        f_ls_phase <- f_ls[phase_index]
        (keyword_ls_phase<- keyword_ls[phase_index])
        

        
        (cell_type_ls <- unique(unlist(lapply(f_ls_phase, function(x) tail(unlist(strsplit(x, split = ' |/')), n=2)[1]))))
        for(cell_type in cell_type_ls){        ## loop2 for each cell type
            ## for group/study in the disease
            df_all=NULL
            for(i in grep(cell_type, f_ls_phase)){ ## loop3 for each dataset in the phase of a cell type
                (f= f_ls_phase[i])
                print(cell_type)
                print(keyword_ls_phase[i])
                print(f)
                
                df = read.delim(f, comment= '#', header = F)
                colnames(df) <- c(cell_type, keyword_ls_phase[i])
                
                if(is.null(df_all)){
                    df_all = df
                }else{
                    df_all = noWarnings(full_join(df_all,df))
                }
            }## loop3 for each dataset in the phase
            df_all$disease = disease
            df_all$phase = phase
            
            assign(paste0(disease, '_',phase,'_', cell_type, '_rotations'), df_all)
        }## loop2 for each cell type
    } ## loop1 for each phase
}#loop 0 for each disease



## save a r data
f_out <- paste0(home_dir, '/ND_results/cell_population/rotation/')
dir.create(f_out, showWarnings = F)
f_out <- paste0(f_out, 'rotation_', Sys.Date(), '.Rdata')

(cmd <- paste0('save(', paste0(grep('_rotation', ls(), value = T), collapse = ','), ', file = "', f_out, '")'))
eval(parse(text = cmd))
    
    