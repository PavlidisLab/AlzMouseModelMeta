#2016-11-14
# make the all_ranks and all_ranks_adj

library(dplyr)
library(HelperFunctions)
library('Hmisc')



## must load adj jackknife first
f_ls=c(paste0(mm_adj_dir, 'jackknife_early_mm_results.Rdata'),paste0(mm_adj_dir, 'jackknife_late_mm_results.Rdata'))
for(f in f_ls){
    load(f)
}
AD_early_cell_adj <- AD_early_mm_variable
HD_early_cell_adj <- HD_early_mm_variable
AD_late_cell_adj <- AD_late_mm_variable
HD_late_cell_adj <- HD_late_mm_variable


# Load results before correction
(f_ls=c(paste0(mm_dir, 'jackknife_early_mm_results.Rdata'),paste0(mm_dir, 'jackknife_late_mm_results.Rdata')))
for(f in f_ls){
    load(f)
}

## load ermineJ data
# 
# load("C:/Users/User/Documents/lab_results/ND_project_combined/git/results/ermineJ11-15.Rdata")


#+++++++++++++++++++ 
# add cell type info
#+++++++++++++++++++ 
## check cell type markers in the top genes
cell_df <- NULL
(cell_f <-max(grep('2016.*Striatum', list.dirs('../../doc/cell_type_markers/',recursive = T, full.names = T), value = T)))
(f_ls <- grep('activation', list.files(cell_f, full.names = T), value = T, invert = T))
for(f in f_ls){
    tmp <- read.delim(f, comment.char = '#')
    tmp$striatum ='yes'
    tmp$disease ='HD'
    if(is.null(cell_df)){
        cell_df <- tmp
    }else{
        cell_df <- rbind(cell_df, tmp)
    }
}

cell_df2 <- NULL
cell_f <- max(grep('2016.*Hippocampus', list.dirs('../../doc/cell_type_markers/',recursive = T, full.names = T), value = T))
(f_ls <- grep('activation', list.files(cell_f, full.names = T), value = T, invert = T))
for(f in f_ls){
    tmp <- read.delim(f, comment.char = '#')
    tmp$hippocampus ='yes'
    tmp$disease ='AD'
    if(is.null(cell_df2)){
        cell_df2 <- tmp
    }else{
        cell_df2 <- rbind(cell_df2, tmp)
    }    
}

cell_df <- noWarnings(full_join(cell_df[, c(1,2,4,5)], cell_df2[, c(1,2,4,5)]))
colnames(cell_df)[1]='geneSymbol'
rm(cell_df2)


#+++++++++++++++++++ 
# add cell type info
#+++++++++++++++++++ 
###### preprocess, add cell type markers

aggregateDF <- function(keyword, cell_df){
    all_ranks = NULL
    need_col = c("disease","phase", "geneSymbol", "down_final_jackknife_rank", "up_final_jackknife_rank",
                 "Name", "pvalue", "pvalue_adj", 
                 "Estimate","down_max_jackknife_rank", 
                 "up_max_jackknife_rank", "Std..Error", "Feature.Type", "Chr", "Strand", "Start", "End", "Human_Disease", 
                 "up_pval", "down_pval", "up_padj", "down_padj", "up_rank", "down_rank", "Input.Type")
    
    for(disease in c('AD', 'HD')){
        for(phase in c('early', 'late')){
            assign('df', eval(parse(text = paste0(disease, '_', phase, keyword))))
            df$disease = disease
            df$phase = phase
            df <- df[,need_col ]%>%droplevels()
            colnames(df)[c(4,5)]= c('down_jack','up_jack')
            
            if(is.null(all_ranks)){
                all_ranks <- df
            }else{
                all_ranks <- rbind(all_ranks,df)
            }
        }
    }
    
    
    #############
    all_ranks <- left_join(all_ranks, cell_df)
    need_col2 = c("disease","phase", "geneSymbol", "down_jack", "up_jack",
                  "Name", "pvalue", "pvalue_adj", 
                  "Estimate","cell_type","down_max_jackknife_rank", 
                  "up_max_jackknife_rank", "Std..Error", "Feature.Type", "Chr", "Strand", "Start", "End", "Human_Disease", 
                  "up_pval", "down_pval", "up_padj", "down_padj", "up_rank", "down_rank", "Input.Type")
    all_ranks <- rmDup(all_ranks[, need_col2])
    return(all_ranks)
    
}


### aggregate all the ranks for both diseases only for the cell adj results
keyword = '_mm_variable'
all_ranks <- aggregateDF(keyword, cell_df)

keyword = '_cell_adj'
all_ranks_adj <- aggregateDF(keyword, cell_df)

# df <- sig_path_genes[, c('geneSymbol',"disease", "phase", "count_bio", "sig_paths_bio","count_all", "sig_paths_all")] %>%droplevels()
# df <- left_join(all_ranks_adj, df)
# df <- df[, c(1:10, 28,27, 30, 29, 19, 11:18)]
# all_ranks_adj <- df



## get the combined ranks
df <- all_ranks_adj
z <- NULL
for (phase in c('early', 'late')){
    y=NULL
    for(disease in c('AD', 'HD')){
        x <- df[which(df$disease == disease & df$phase == phase), ]
        x <- x[,c(2:5)]
        colnames(x)[3:4] <- paste0(disease, colnames(x)[3:4])
        if(is.null(y)){
            y <- x
        }else{
            y <- na.omit(left_join(y, x))
        }
    }
    ## combine ranks for a phase
    colnames(y)
    y$combined_up <- rank((y$ADup_jack+y$HDup_jack)/2)
    y$combined_down <- rank((y$ADdown_jack+y$HDdown_jack)/2)
    if(is.null(z)){
        z <- y
    }else{
        z <- rbind(z, y)
    }
}

z <- z[, c(1,2,7,8)]
df <- left_join(df, z)
all_ranks_adj <- df


## reorder some columns
all_ranks$P_adj <- round(all_ranks$pvalue_adj, digits = 3)
all_ranks_adj$P_adj <- round(all_ranks_adj$pvalue_adj, digits = 3)
all_ranks <- all_ranks[, colnames(all_ranks)[c(1:6, 27, 9:19, 7:8, 20:26)]]
all_ranks_adj <- all_ranks_adj[, colnames(all_ranks_adj)[c(1:6, 29, 9:19, 7:8, 20:28)]]
