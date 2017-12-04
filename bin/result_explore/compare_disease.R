#' 2016-08-11
#' compare diseases (results from meta, jack and mixed model)
#' 1. get the average rank of each comparison (same phase, same regulation), and save the results
#' 2. look for top overlapped genes
# 
# home_dir= '/home/bzhuang/ND_project_combined/'
#  
# rm(list=setdiff(ls(),'home_dir'))
# 
# 
# disease_ls= c('AD','HD','PD')
# model_ls = c('random_intercept', 'random_slope')
# model_keyword_ls = c('_include_NA')
# out_dir <- paste0(home_dir, '/ND_results/Disease_comparison/')



###############
## define values
###############
source("rank_product/rankprodbounds.R")
library(HelperFunctions)
library(dplyr)

dir.create(out_dir,recursive = T, showWarnings = F)
value_col <- c('up_pval', 'down_pval','adj_combined_max_p', 'Fisher')
all_col <- c(value_col, 'geneSymbol','geneNames')


###############
## grep all up and down regulation all phase of jackknife,  meta, and mixed model (specify model name and model keyword)
###############

for(disease in disease_ls){
    source('config/config_wrappers.R')
    ## get the meta, jack values and mixed models
    ## get meta jack for combined geneotypes
    
    ## meta and jack
    (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                              paste0(disease_dir, 'meta_analysis/meta_jack/')))
    for(jack_meta_folder in jack_meta_folder_ls){
        mj_f_ls <- grep('regulation', list.files(jack_meta_folder, full.names = T), value = T)
    }
    
    ## mixed models
    for(model in model_ls){
        for(model_keyword in model_keyword_ls){
            mm_dir = paste0(disease_dir, 'mixed_model/',model, model_keyword,'/') ## mixed model dir(parent dir)
            ## get results
            (file_ls <- grep('regulation.*.tsv', list.files(mm_dir, recursive = T, full.names = T), value = T))
            
            ##combine all files
            mj_f_ls <- c(mj_f_ls, file_ls)
        }
    }

    
    ## remove archival
    mj_f_ls <- sort(grep('archive|archival', mj_f_ls, invert = T, value = T))
    
    eval(parse(text = paste0('disease_', disease, '_ls = mj_f_ls')))
}

###############
## read all files
###############
input_ls <- vector()
result_ls <- vector()
top_ls <- vector()

## including overlap genes from each ND

for(i in 1:length(disease_AD_ls)){  ## loop for all files in disease_*D_ls
    df_all <- NULL
    f_all <- NULL
    print(i)
    for (disease in disease_ls){

        d_file <- eval(parse(text = paste0('disease_', disease, '_ls[i]')))
        print(d_file)
        
        df <- read.delim(d_file, comment.char = '#')
        (target_col <- intersect(value_col,colnames(df)))
        df <- df[, intersect(all_col,colnames(df))]

        colnames(df)[which(colnames(df) == target_col)] <- paste0(disease, '_', target_col, '_value')
        if(is.null(df_all)){
            df_all <- df
        }else{
            if('geneNames' %in% colnames(df)){
                df <- df[, -which(colnames(df) == 'geneNames')] %>%droplevels()
            }
            df_all <- noWarnings(full_join(df_all, df, by = c('geneSymbol')))
                
#                 df_all <- noWarnings(full_join(df_all, df, by = c('geneSymbol', 'geneNames')))
#             }else{
#                 df_all <- noWarnings(full_join(df_all, df, by = c('geneSymbol')))
#             }
                
            
        } 
        f_all <- c(f_all, d_file)
    }
    
    #####
    ## remove NA genes and rank genes for each disease, get the mean rank and rank product
    #####
    df_rp <- na.omit(df_all)   ## remove genes with missing value from any of the NDs
    (n_col <- grep('_value', colnames(df_rp), value = T))
    for(j in n_col){
        df_rp$new <- rank(df_rp[, j], na.last = 'keep')
        colnames(df_rp)[colnames(df_rp) == 'new'] <- gsub('_value', '_rank', j)
    }
    
    ## calculate p value for rankproducts k: number of input disease
    (index <- grep("_rank", colnames(df_rp))) ## grep all ranks
    df_rp$rank_avg <- apply(df_rp[, index], 1, mean)
    
    rownames(df_rp) <- df_rp$geneSymbol
    
    ## get the overlap genes in top threshold
    df_top=NULL
    gene_pool=NULL
    seq_i <- seq(from = 0, to = 1000, by=100)
    for(j in 1:length(seq_i)-1){
        (low_threshold = seq_i[j])
        (top_threshold = seq_i[j+1])
        topgenes_overlap <- df_rp$geneSymbol
        for(rank_index in index){
            (topgenes <- df_rp[which(df_rp[, rank_index]<=top_threshold), 'geneSymbol'])
            (topgenes_overlap <- intersect(topgenes_overlap,topgenes))
        }
        
        ## record which threshold has overlapped genes
        #print(gene_pool)
        #print(topgenes_overlap)
        topgenes_overlap <- setdiff(topgenes_overlap, gene_pool)
        if(length(topgenes_overlap)>0){
            df_tmp <- df_rp[topgenes_overlap, ]
            df_tmp$top_threshold = top_threshold
            if(is.null(df_top)){
                df_top = df_tmp
            }else{
                df_top <- rbind(df_top,df_tmp)
            }
        }
        gene_pool <- union(gene_pool, topgenes_overlap)
        print(gene_pool)
    }
    
    
    
    df_rp$rank_avg_percent <- df_rp$rank_avg/max(df_rp$rank_avg, na.rm = T)  # will be used for ermineJ
    df_rp$rank_product <- apply(df_rp[, index], 1, prod)
    
    ## gene count
    (n_gene <- nrow(df_rp))
    df_rp$RP_p_upper <- rankprodbounds(rho = df_rp[, "rank_product"], n = n_gene, k= 3, Delta= 'upper')
    df_rp$RP_p_lower <- rankprodbounds(rho = df_rp[, "rank_product"], n = n_gene, k= 3, Delta= 'lower')
    df_rp$RP_p_geometric <- rankprodbounds(rho = df_rp[, "rank_product"], n = n_gene, k= 3, Delta= 'geometric')
    df_rp$PR_rank <- rank(df_rp[, 'RP_p_geometric'], na.last = 'keep')
    df_rp$avg_rank <- rank(df_rp[, 'rank_avg'], na.last = 'keep')
    df_rp <- df_rp[, setdiff(colnames(df_rp),n_col)]
    
    ##

    
    ## combine with all the genes
    df_all <- left_join(df_all, df_rp, by='geneSymbol')
    df_all <- df_all[, c('geneSymbol', 'rank_avg_percent', setdiff(colnames(df_all),c('geneSymbol', 'rank_avg_percent')))]
    df_top <- df_top[, c('geneSymbol', 'top_threshold', setdiff(colnames(df_top),c('geneSymbol','top_threshold')))]
    
    
    result_ls[i] <- list(df_all)
    input_ls[i] <- list(f_all)
    top_ls[i] <- list(df_top)
}


##########
## analysis
##########

out_top_dir <- paste0(out_dir, '/top_overlap_genes/')
dir.create(out_top_dir,recursive = T, showWarnings = F)
for (k in 1: length(result_ls)){
    print(input_ls[k])
    df <- result_ls[[k]]
    df <- df[with(df, order(avg_rank)), ]  # order by the average rank
    (f_out <- paste0(grep('meta_genes|jackknife|random|.tsv', unlist(strsplit(input_ls[k][[1]][1], split = '/')), value = T), collapse='_'))
    (f_top_out <- gsub('.tsv', '_top_overlap.tsv', f_out))
    (f_out <- paste0(out_dir, f_out))
    (f_top_out <- paste0(out_top_dir, f_top_out))
    msg <- paste0('# ', Sys.Date(),
                  '\n# Input diseases: ', paste0(disease_ls, collapse = ', '), 
                  '\n# ', paste0(input_ls[[k]], collapse = '\n# '))
    
    writeTable(df, f_out = f_out, msg = msg)
    print(f_out)
    
    df_top <- top_ls[[k]]
    writeTable(df_top, f_out = f_top_out, msg = msg)
    print(f_top_out)
}




# 
# k=1
# input_ls[k]
# early_d <- result_ls[[k]]
# early_d <- early_d[with(early_d, order(PR_rank)), ]
# writeTable(early_d[1:10, ], f_out = 'early_down_top.tsv')
# 
# k=10
# input_ls[k]
# early_u <- result_ls[[k]]
# early_u <- early_u[with(early_u, order(PR_rank)), ]
# writeTable(early_u[1:10, ], f_out = 'early_up_top.tsv')
# k=11
# input_ls[k]
# late_d <- result_ls[[k]]
# 
# k=12
# input_ls[k]
# late_u <- result_ls[[k]]
