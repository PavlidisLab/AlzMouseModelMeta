home_dir_ls <- c('C:/Users/User/Documents/lab_results/ND_project_combined/',
                 '/home/bzhuang/ND_project_combined/')

home_dir <- home_dir_ls[which(dir.exists(home_dir_ls))]
rm(list=setdiff(ls(),'home_dir'))
 
library(HelperFunctions)
library(dplyr)
library(ggplot2)

 

###########################
######### make all ranks and all ranks adj for mixed model results
###########################
# mm_dir <- paste0(home_dir, "/ND_results/mm_results/2017-02-02/")
# mm_adj_dir <- paste0(home_dir, "/ND_results/mm_results_cell_adj/2017-02-02/")
# # # make the all_ranks and all_ranks_adj
# source('./thesis_stuff/mk_all_ranks_all_ranks_adj.R')
# ## make GO_non_adj and GO_adj for mixed models:
# folder_pre = 'ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm'
# f_result_out =paste0('../../results/emineJ_', Sys.Date(),'.Rdata')
# go_plot_out <- '../../results/ND_results/'
# source('thesis_stuff/Go_enrichment_summary.R')
# save(all_ranks, all_ranks_adj, GO_adj, GO_non_adj, file = '../../results/ranks_paths_2017-02-02.Rdata')

###########################


###########################
######### make all ranks and all ranks adj for FIXED model results ########
###########################
# mm_dir <- paste0(home_dir, "/ND_results/fm_results/2017-02-17/")
# mm_adj_dir <- paste0(home_dir, "/ND_results/fm_results_cell_adj/2017-02-17/")
# # # make the all_ranks and all_ranks_adj
# source('./thesis_stuff/mk_all_ranks_all_ranks_adj.R')
# fm_all_ranks <- all_ranks
# fm_all_ranks_adj <- all_ranks_adj
## get
##
# ## make GO_non_adj and GO_adj for FIXED models: 2017-02-21
# folder_pre = 'ermineJ/mixed_model_jackknife/fixed_model_include_NA_low_exp_rm'
# go_plot_out = '../../results/ND_results/GO_figures_fixed_model/'
# f_result_out =paste0('../../results/emineJ_fixed_models_', Sys.Date(),'.Rdata')
# source('thesis_stuff/Go_enrichment_summary.R')
# fm_GO_adj <- GO_adj
# fm_GO_non_adj <- GO_non_adj
# save(fm_all_ranks, fm_all_ranks_adj,  fm_GO_adj, fm_GO_non_adj, file = '../../results/fixed_model_ranks_2017-02-17.Rdata')

# load('../../results/fixed_model_ranks_2017-02-17.Rdata')

###########################
# 2017-02-01 not in use anymore
###########################
# # # add the significantly enriched pathways # 11-15, genset gene set 500, iteration 200k
# load("../../results/ermineJ11-26_sig_genes.R")


###########################
# load saved Rdata
###########################
rm(list=setdiff(ls(),'home_dir'))
 

## load all ranks, adj, and enrichment results for adj
load('../../results/ranks_paths_2017-02-02.Rdata')




#' NOTES:
#' GO_non_adj : go enrichment for non-cell_pop-corrected jackknife results FDR <0.05 results
#' GO_adj : go enrichment for cell_pop-corrected jackknife results FDR <0.05 results
#' all_ranks: NON-corrected jackknife results
#' all_ranks_adj: cell_pop-corrected jackknife results




## load cell marker rotation estimations
#'summary of cell markers and estimation of rotations and variations
# source('thesis_stuff/check_cell_pop_rotation.R')
load(('../../results/cell_marker_summary.Rdata'))

## source helper functions
source('./thesis_stuff/check_genes_helpers.R')

#+++++++++++++++++++
# get the gender gene non-expression median (reason for the threshold) # must use the expression.R in the 
# /mixed_model/random_intercept_include_NA/
source('thesis_stuff/sex_gene_threshold.R')
#+++++++++++++++++++

#+++++++++++++++++++ thesis table
#' summary of all samples and genes counts for input data
#' genes are after filter: input data from mixed_model/random_intercept_include_NA_low_exp_rm/
source('summary_tables/summary_mm_input_data.R')
#+++++++++++++++++++

#+++++++++++++++++++ 
# 2017-03-01 # thesis plots # updated 2017-03-21
# plot the expression of a gene (or gene list) from each study
# with raw expression, study corrected or study and MGP corrected
#+++++++++++++++++++ 
source('thesis_stuff/plot_indi_gene_exp_per_study.R')


f_out_dir <- paste0('../../results/ND_results/gene_expr_before_after_MGP/', Sys.Date(), '/')
dir.create(f_out_dir, showWarnings = F,recursive = T)

#++
disease <- 'AD'
phase <- 'early'
#++
geno_f <- paste0('../../doc/',disease, '_mouse_dataset_doc/',disease,'_', phase, '_design_model_types.tsv')

x <- arrayData(disease,phase,geno_f=geno_f)
array_dat_raw <- x[[1]]
array_dat_MGP <- x[[2]]
array_dat_study <- x[[3]]
array_design <- x[[4]]

### compare before and after MGP ( 1 gene a plot) ## plot
gene_ls <- c('Sqle','Msmo1', 'Nsdhl')
df <- mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,return_df = T, 
               one_plot_font_size =9,y_axis_null =F)
# 


#++
disease <- 'AD'
phase <- 'late'
#++
geno_f <- paste0('../../doc/',disease, '_mouse_dataset_doc/',disease,'_', phase, '_design_model_types.tsv')
x <- arrayData(disease,phase,geno_f=geno_f)
array_dat_raw <- x[[1]]
array_dat_MGP <- x[[2]]
array_dat_study <- x[[3]]
array_design <- x[[4]]

### compare before and after MGP ( 1 gene a plot) ## plot
gene_ls <- c('Trem2','Msmo1')
gene_ls <- c('Sqle','Msmo1', 'Nsdhl')
df <- mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,return_df = T, 
               one_plot_font_size =9,y_axis_null =F)



#++
disease <- 'HD'
phase <- 'early'
#++
geno_f <- paste0('../../doc/',disease, '_mouse_dataset_doc/',disease,'_', phase, '_design_model_types.tsv')
x <- arrayData(disease,phase,kuhn =T,geno_f=geno_f)
array_dat_raw <- x[[1]]
array_dat_MGP <- x[[2]]
array_dat_study <- x[[3]]
array_design <- x[[4]]

### compare before and after MGP ( 1 gene a plot) ## plot
gene_ls <- c('Coa3', 'Ndufs3')
df <- mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,return_df = T, 
               one_plot_font_size =9,y_axis_null =F)



#++
disease <- 'HD'
phase <- 'late'
#++
geno_f <- paste0('../../doc/',disease, '_mouse_dataset_doc/',disease,'_', phase, '_design_model_types.tsv')
x <- arrayData(disease,phase,kuhn =T,geno_f=geno_f)
array_dat_raw <- x[[1]]
array_dat_MGP <- x[[2]]
array_dat_study <- x[[3]]
array_design <- x[[4]]

### compare before and after MGP ( 1 gene a plot) ## plot
gene_ls <- c('Ddit4l', 'Nrep')
df <- mainPlot(gene_ls, array_dat_raw, array_dat_MGP, disease, array_design,f_out_dir,return_df = T, 
               one_plot_font_size =9,y_axis_null =F)


#+++++++++++++++++++ 
# checking gene ranking changes before and after MGP
#' genes that are FDR <0.05, or under gene ranking threshold
#+++++++++++++++++++ 

fdr = 0.05
threshold =52
## preprocess
df1 <- all_ranks[, 1:9]
df2 <- all_ranks_adj[, 1:8]
colnames(df1)[c(4:5,7,8)] <- paste0('before_',colnames(df1)[c(4:5,7,8)])
df1$cell_type <- as.character(df1$cell_type)
df1$cell_type[which(is.na(df1$cell_type))] <- 'not_cell_type'
df1$before_regulation = ''
df1$before_regulation[which(df1$before_Estimate >0)] = 'up'
df1$before_regulation[which(df1$before_Estimate <0)] = 'down'
colnames(df2)[c(4:5,7,8)] <- paste0('after_',colnames(df2)[c(4:5,7,8)])
df <- dplyr::left_join(df1, df2)

## function
summary_DE_b_a <- function(df_t, notes, threshold=NULL, col_prefix=NULL){
    #col_prefix; 'after', 'before', 'both'
    
    if(!is.null(threshold)){
       
        if(col_prefix =='both'){
            index <- union(which(df_t[, paste0('before_up_jack')] <= threshold), 
                           which(df_t[, paste0('before_down_jack')] <= threshold))
            
            (index <- intersect(index, union(which(df_t[, paste0('after_up_jack')] <= threshold), 
                           which(df_t[, paste0('after_down_jack')] <= threshold))))
        }else{
            index <- union(which(df_t[, paste0(col_prefix, '_up_jack')] <= threshold), 
                           which(df_t[, paste0(col_prefix, '_down_jack')] <= threshold))
        }
            
            
            if(length(index) >0){
                df_t <- df_t[index, ]
            }else{
                stop(print('no genes under threshold'))
            }
        
    }
    
    df_summary <- data.frame(table(df_t[c(1,2,9,10)]))
    df_summary$note <- notes
    df_threshold <- paste0(col_prefix, '_', threshold)
    ## rm 0
    df_summary <- df_summary[which(df_summary$Freq >0), ]
    df_summary <- orderCol(df_summary,c('disease', 'phase', 'cell_type'))
    
    return(list(df_summary, df_t))
}




## DE before and DE after
df_de_ba <- df[which(df$before_P_adj <=fdr & df$after_P_adj <=fdr), ]

#**********
## DE before and DE after: changed direction fdr
#**********
tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate <0), ]
notes <- paste0('DE_before_DE_after_changed_dir_FDR', fdr)
x <- summary_DE_b_a(tmp, notes)
df_t <- x[[2]]
y <- x[[1]]
print(y)



#**********
## DE before and DE after: changed direction fdr and threshold (before) # no genes
#**********
# tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate <0), ]
# notes <- paste0('DE_before_DE_after_changed_dir_FDR', fdr)
# x <- summary_DE_b_a(tmp, notes, threshold = 52, col_prefix = 'before')
# df_t <- x[[2]]
# y <- x[[1]]
# print(y)


#**********
## DE before and DE after: changed direction fdr and threshold (after)
#**********
prefix <- 'after'
tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate <0), ]
notes <- paste0('DE_before_DE_after_changed_dir_FDR', fdr, '_threshold_',prefix, threshold)
x <- summary_DE_b_a(tmp, notes, threshold = threshold, col_prefix = prefix)
df_t <- x[[2]]
y <- x[[1]]
print(y)


#**********
## DE before and DE after: same direction
#**********
tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate >0), ]
notes <- paste0('DE_before_DE_after_same_dir_FDR', fdr)
x <- summary_DE_b_a(tmp, notes)
print(x)

#**********
## DE before and DE after: same direction and threshold after
#**********
prefix <- 'after'
tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate >0), ]
notes <- paste0('DE_before_DE_after_same_dir_FDR', fdr, '_threshold_',prefix, threshold)
x <- summary_DE_b_a(tmp, notes, threshold = threshold, col_prefix = prefix)
df_t <- x[[2]]
y <- x[[1]]
print(y)

#**********
## DE before and DE after: same direction and threshold before
#**********
prefix <- 'before'
tmp <- df_de_ba[which(df_de_ba$before_Estimate * df_de_ba$after_Estimate >0), ]
notes <- paste0('DE_before_DE_after_same_dir_FDR', fdr, '_threshold_',prefix, threshold)
x <- summary_DE_b_a(tmp, notes, threshold = threshold, col_prefix = prefix)
df_t <- x[[2]]
y <- x[[1]]
print(y)


#**********
## DE before and not DE after
#**********
tmp <- df[which(df$before_P_adj <=fdr & df$after_P_adj >fdr), ]
notes ='DE_before_not_DE_after'
x <- summary_DE_b_a(tmp, notes)
print(x)


#**********
## DE after and not DE before, threshold after
#**********
prefix <- 'after'
tmp <- df[which(df$before_P_adj >fdr & df$after_P_adj <=fdr), ]
notes <- paste0('not_DE_before_DE_after_FDR', fdr, '_threshold_',prefix, threshold)
x <- summary_DE_b_a(tmp, notes, threshold = threshold, col_prefix = prefix)
df_t <- x[[2]]
y <- x[[1]]
print(y)





#+++++++++++++++++++
# compare my top genes and the reported genes (thesis tables)
#+++++++++++++++++++
source('thesis_stuff/compare_pub_thesis.R')
threshold_ls =c(20,50)
phase_ls = c('early', 'late')
disease_ls =c('AD', 'HD')

out_dir <- (paste0('../../results/ND_results/compare_to_original_study/', Sys.Date(), '/'))

for(phase in phase_ls){
    for(threshold in threshold_ls){
        for(disease in disease_ls){

            ## get the reported genes from repo
            (gene_list_dir <- paste0('../../doc/original_publication_gene_list/',disease, '/'))
            ## get the most recent results from adj and non adj
            keyword=''
            (input_dir <- max(list.dirs(paste0('../../results/ND_results/top_gene_heatmaps/', disease, 
                                               '/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm',keyword, '/'), recursive = F)))
            f_exp_non <- paste0(input_dir,'/',
                                phase, '/', disease, '_', phase, '_inputdata.tsv')
            
            keyword='_adj_cell_pop'
            (input_dir <- max(list.dirs(paste0('../../results/ND_results/top_gene_heatmaps/', disease, 
                                               '/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm',keyword, '/'), recursive = F)))
            f_exp_adj <-paste0(input_dir,'/',
                                phase, '/', disease, '_', phase, '_inputdata.tsv')
 
            print(f_exp_non)
            x <- comparePublication(df=all_ranks, threshold, phase, disease, gene_list_dir, prefix ='_',
                                    out_dir = out_dir,write_out = T,
                                    f_exp = f_exp_non)
            print(f_exp_adj)
            y <- comparePublication(df = all_ranks_adj, threshold, phase, disease, gene_list_dir, prefix ='_adj', 
                                    out_dir = out_dir,
                                    write_out = T, f_exp = f_exp_adj)
        }
    }
}


#+++++++++++++++++++
# plot pvalue distribution
#+++++++++++++++++++
##############
## thesis plot - histo plot for pvalues before and after cell type correction
##############

## plot pvalues
df <- all_ranks[, c('disease', 'phase','pvalue')]%>%droplevels()
df$correction <- 'Before'
df2 <- all_ranks_adj[, c('disease', 'phase','pvalue')]%>%droplevels()
df2$correction <- 'After'
df <- rbind(df,df2)
df$phase <- paste0(df$disease, ' ', df$phase)
df <- as.data.frame(unclass(df))
df$correction <- factor(df$correction, levels= c('Before', 'After'))


plot_dir <- '../../results/ND_results/figures/pvalues/'
dir.create(plot_dir, recursive = T,showWarnings = F)
for(disease in c("AD", "HD")){
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_', disease, '_pvalue.png') 
    plotPvalue(df2, one_plot_font_size =14, f_out)
}

## save as svg for later modification
for(disease in c("AD", "HD")){
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_', disease, '_pvalue.svg') 
    plotPvalue(df2, one_plot_font_size =14, f_out)
}


#+++++++++++++++++++
##### get the enriched GO terms and cross disease overlap terms -thesis
#' 
#+++++++++++++++++++
outdir <- '../../results/ND_results/sig_go_terms/'
dir.create(outdir,recursive = T, showWarnings = F)
df <- GO_adj
df <- df[, c(1:3,19, 5:9, 11, 15:16)]
## get the ranking of the top genes in the gene sets

countTop <- function(x, threshold=100,return_count=T, return_gene=F){
    x <- as.character(x)
    x_char <- unlist(strsplit(x, ', '))
    x <- strsplit(x_char, '\\(|\\)')
    y=NULL
    for(i in 1:length(x)){
        y <- c(y, as.numeric(x[[i]][2]))
    }
    y_len <- length(which(y <= threshold))
    
    if(return_count) return(y_len)
    if(return_gene) return(paste0(x_char[1:y_len],collapse = ', '))
    
}


df$top_count_20 <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =20))
df$top_count_20_gene <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =20, return_count = F, return_gene = T))
# df$top_count_100 <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =100))
# df$top_count_200 <- apply(as.matrix(df$top_genes), 1, function(x) countTop(x, threshold =200))
df$top_20 <- paste0(df$NumGenes, ' (',df$top_count_20, ')')


## MP less than 0.1
#df_mp <- df[which(df$MP_adj <= 0.1), ]



## remove pathways if there's no top 20 genes involved
df <- df[which(df$top_count_20 >=2), ]%>%droplevels()


## get the repeated top 20 genes
overlap_term <- data.frame(top_count_20_gene=df$top_count_20_gene[duplicated(as.character(df$top_genes))])
overlap_term$overlap_top_20 = 'overlap_20'
overlap_term <- rmDup(overlap_term)

df_out <- noWarnings(left_join(df,overlap_term))



## overlapped terms across disease, 
phase_ls =c('early', 'late')
tmp <- df_out

overlap_df_all <- NULL
for(phase in phase_ls){
    overlap_term <- intersect(as.character(tmp$Name[which(tmp$phase == phase & tmp$disease=='AD')]), as.character(tmp$Name[which(tmp$phase == phase & tmp$disease=='HD')]))
    
    if(length(overlap_term) >0){
        overlap_df <- data.frame(Name = overlap_term,
                                 phase = phase)
        if(is.null(overlap_df_all)){
            overlap_df_all <- overlap_df
        }else{
            overlap_df_all <- rbind(overlap_df_all,overlap_df)
        }
    }
}
overlap_df_all$cross_disease_overlap = 'overlap'

df_out <- noWarnings(left_join(df_out,overlap_df_all))


## if top genes are overlapped, maybe select the lowest multifunctionality
tmp <- df_out[,c('Multifunctionality', 'Name','overlap_top_20', 'top_count_20_gene', 'NumGenes' )]
tmp$agg <- paste0(tmp$top_count_20_gene, tmp$overlap_top_20)

df2 <- aggregate(Multifunctionality ~ agg, tmp[,c(1:4,6)], function(x) min(x, na.rm = T))
df2$min_multi <- 'yes'
tmp <- noWarnings(left_join(tmp, df2))

## get the pathways with min gene set number
df3 <- aggregate(NumGenes ~ agg, tmp[,c(1:3,5:6)], function(x) min(x, na.rm = T))
df3$min_n_genes <- ''
tmp <- noWarnings(left_join(tmp, df3))

df_out <- left_join(df_out,tmp)

## final filter:: remove pathways has the same top 20 genes with higher multifunctionality
df_min <- df_out[which(!is.na(df_out$min_multi)),]
df_min <- orderCol(df_min, c('disease', 'phase', 'regulation', 'Multifunctionality','min_n_genes'))
## replace NA with ''
index <- which(is.na(df_min$min_n_genes))
df_min$min_n_genes[index] <- 'No'
index <- which(is.na(df_min$overlap_top_20))
df_min$overlap_top_20[index] <- ''
index <- which(is.na(df_min$cross_disease_overlap))
df_min$cross_disease_overlap[index] <- ''

df_min$regulation = paste0(df_min$phase, ': ', df_min$regulation, '-regulated')
df_min$tmp <- as.factor(paste0(df_min$disease,df_min$regulation))

index <- NULL
for(keep in levels(df_min$tmp)){
    index <- c(index, which(df_min$tmp == keep)[1])
}

## get the overlap across diseases
df_o <- df_min[which(df_min$cross_disease_overlap =='overlap'), ]
df1 <- df_o[which(df_o$disease == 'AD'), c('phase','ID', 'Name','top_count_20_gene')]
colnames(df1)[4]='AD'
df2 <- df_o[which(df_o$disease == 'HD'), c('phase','ID', 'Name','top_count_20_gene')]
colnames(df2)[4]='HD'
df_o <- left_join(df1, df2)
f <- paste0(outdir, Sys.Date(), "_cross_dis_sig_go_terms.tsv")
writeTable(df_o, f_out = f)

## save tables
df_min$disease <- as.character(df_min$disease)
df_min$disease[setdiff(1:nrow(df_min),index)] <- ''
df_min$regulation <- as.character(df_min$regulation)
df_min$regulation[setdiff(1:nrow(df_min),index)] <- ''


need_col <- c('disease','regulation','phase','ID','Name','P_adj','Multifunctionality','top_20',
              'top_count_20_gene','cross_disease_overlap','min_n_genes')
new_col <- c('Disease','Regulation','phase', 'GO Term ID','GO Term Description','FDR','Multifunctionality Score','Number of Genes',
         'Top Hits','cross_disease_overlap','min_n_genes')


df_min <- df_min[,new_col]
colnames(df_min) <- new_col



f <- paste0(outdir, Sys.Date(), "_sig_go_terms.tsv")
writeTable(df_min, f_out = f)

#+++++++++++++++++++
##### check the # of genes up or down regulated within padj<threshold
#' the correlation (all genes) before and after adj
#' write a table and do the thesis bar plot for each disease
#+++++++++++++++++++
#' HD late has large # of genes up or down regulated
#' 
#' 
#### check correlations of ranks before and after 
source('./thesis_stuff/check_genes_helpers.R')
df_all_count <- NULL

for (disease in c('AD', 'HD')){
    for(phase in c('early', 'late') ){
        for(regulation in c('up', 'down')){
            print('compare cell adj vs non cell adj correclations ')
            df <- compareAdjCor(disease, phase, regulation)
            if(is.null(df_all_count)){
                df_all_count <- df
            }else{
                df_all_count <- rbind(df_all_count,df)
            }
        }
    }
}


####################################
#### check correlations of ranks early and late
####################################
df_all_el <- NULL
for (disease in c('AD', 'HD')){
            print('compare cell early and late correclations ')
            df <- compareELCor(all_df = all_ranks_adj, disease, prefix = 'adj',gene_names=F)
            df2 <- compareELCor(all_df = all_ranks, disease, prefix = '',gene_names=F)
            if(is.null(df_all_el)){
                df_all_el <- rbind(df, df2)
            }else{
                df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
            }
}
## check top genes early late correction
top_t_ls = c(100, 200)
for (disease in c('AD', 'HD')){
    print('compare cell early and late correclations ')
    for(regulation in c('up', 'down')){
        for(top_t in top_t_ls){
            df2 <- compareELCor(all_df = all_ranks[which(all_ranks[paste0(regulation, '_jack')] <= top_t), ], disease, prefix = paste0('top_', regulation, top_t),gene_names=T)
            df <- compareELCor(all_df = all_ranks_adj[which(all_ranks_adj[paste0(regulation, '_jack')] <= top_t), ], disease, prefix = paste0('adj_top_', regulation, top_t),gene_names=T)
            df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
        }
    }
}

### compare between AD and HD
for(phase in c('early', 'late')){
    print('compare AD and HD ')
    df <- compareADHDCor(all_df = all_ranks_adj, disease,  phase, prefix = 'adj')
    df$gene_names=""
    df2 <- compareADHDCor(all_df = all_ranks, disease,  phase, prefix = '')
    df2$gene_names=""
    
    df_all_el <- Reduce('rbind', list(df_all_el,df,df2))
    
}




outdir <- '../../results/ND_results/early_late_corr/'
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, Sys.Date(), "_early_late_corr.tsv")
writeTable(df_all_el, f_out = f)

#+++++++++++++++++++
########## thesis plot [not in use]
####### summary of each mouse model is in how many studies
#+++++++++++++++++++
df <- read.delim('../../doc/AD_mouse_dataset_doc/dataset_info_mixed_model.tsv', comment.char = '#')
df$Dataset <- gsub('\\.1|\\.2', '',df$Dataset)
df2 <- read.delim('../../doc/HD_mouse_dataset_doc/dataset_info_mixed_model.tsv', comment.char = '#')
df2$Genotype <- gsub('R62_300Q', 'R62', df2$Genotype)
df <- rmDup(rbind(df[, c('Dataset', 'Genotype')],df2[, c('Dataset', 'Genotype')]))

(df <- as.data.frame(table(df$Genotype)))
df <- as.data.frame(table(df$Freq))

ggplot(df, aes(Freq)) + geom_bar()

one_plot_font_size=10
f_out <- paste0('../../results/ND_results/figures/number_of_studies_per_mouse_model_', Sys.Date(), '.png')

p <- ggplot(df, aes(x = Var1)) +
    geom_bar(stat="identity", aes_string(y='Freq'), position="dodge") +
    geom_text(aes_string(x='Var1', y='Freq', label='Freq'),nudge_y = 0.9, size =3.5) +
    xlab('Number of studies per mouse model') +
    ylab('Count') +
    theme_bw() +
    
    theme(text=element_text(family = 'Arial'),
          #legend.position="none",
          axis.title.x=element_text(size = one_plot_font_size),
          axis.title.y=element_text(size = one_plot_font_size),
          title = element_text(size = one_plot_font_size, colour = 'black'),
          axis.text.x = element_text(size = 10), ## vjust 0.5 put x labels in the middle
          axis.text.y=element_text(size = 10),
          strip.text.x = element_text(size = one_plot_font_size),
          strip.text.y = element_text(size = one_plot_font_size))
ggsave(filename = f_out, plot=p, width = 6.45, height =3.15, units = "in")  

################################
########
## thesis bar plot DE genes before and after cell type correction
########
source('./thesis_stuff/check_genes_helpers.R')
threshold <- 0.05
df <- upDownCount(all_ranks, threshold)
df <- df[,c('disease','phase', 'regulation', 'Freq', 'ratio')]
df$correction <- 'Before'
df_before <- df

df <- upDownCount(all_ranks_adj, threshold)


df <- df[,c('disease','phase', 'regulation', 'Freq', 'ratio')]
df$correction <- 'After'
df_after <- df

df <- rbind(df_before, df_after)
df$Count <- paste0(df$Freq, '(', round(df$ratio*100,2), '%)')
df <- left_join(df, df_all_count)
df <- as.data.frame(unclass(df))
df$correction <- factor(df$correction, levels= c('Before', 'After'))

# get no direction
df1 <- aggregate(Freq ~ disease + phase +correction , data = df, FUN = sum)
df2 <- aggregate(ratio ~ disease + phase +correction , data = df, FUN = sum)
df_no_reg_dir <- left_join(df1, df2)
df_no_reg_dir$Count <- paste0(df_no_reg_dir$Freq, '(', round(df_no_reg_dir$ratio*100,2), '%)')


plot_dir <- '../../results/ND_results/pvalues/'
dir.create(plot_dir, recursive = T,showWarnings = F)
writeTable(df, f_out = paste0(plot_dir, Sys.Date(), '_padj5_count.tsv'))
writeTable(df_no_reg_dir, f_out = paste0(plot_dir, Sys.Date(), '_padj5_count_no_dir.tsv'))
########
## thesis bar plot DE genes before and after cell type correction
########
for(disease in c('AD', 'HD')){
    ## x label is regulation (up or down, facet by phase)
    df2 <- filterContain(df, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar.png') 
    p <- plotBarDE(df2, one_plot_font_size=10, f_out, x_col = 'regulation',
              y_col = 'Freq' ,
              x_label = 'Regulation',
              y_label = 'Count',
              return_p = T)
}



for(disease in c('AD', 'HD')){
    df2 <- filterContain(df_no_reg_dir, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar_no_dir.png') 
    p <- plotBarDE(df2, one_plot_font_size=14, f_out, x_col = 'phase',
                   y_col = 'Freq' ,
                   x_label = 'Disease Phase',
                   y_label = 'DE Gene Count',
                   return_p = T, facet_phase =F, save_p = T)
}

##save as svg to move the lables around
for(disease in c('AD', 'HD')){
    df2 <- filterContain(df_no_reg_dir, column = 'disease', value = disease)
    f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar_no_dir.svg') 
    p <- plotBarDE(df2, one_plot_font_size=14, f_out, x_col = 'phase',
                   y_col = 'Freq' ,
                   x_label = 'Disease Phase',
                   y_label = 'DE Gene Count',
                   return_p = T, facet_phase =F, save_p = T)
}


# for(disease in c('AD', 'HD')){
#     df2 <- filterContain(df, column = 'disease', value = disease)
#     f_out <-paste0(plot_dir,Sys.Date(),'_',disease,'_pvalue_bar.svg') 
#     plotBarDE(df2, one_plot_font_size=10, f_out)
# }

    



#+++++++++++++++++++
# look at genes
#+++++++++++++++++++
gene = 'Prkcq'
gene = 'Gnai2'
gene = 'Gigyf2'
gene = 'Trpc1'
gene = 'Acly'
gene = 'Sqle'
gene = 'Gldc'

## AD genes
gene = c('C1qa', 'C1qb', 'C1qc')
gene = c('Picalm', 'Cr1','Clu')
gene = c('Ide') 
gene = c('Pld3', 'Unc5c', 'Akap9')

gene = c('Drd2', 'Adora2a', 'Cnr1', 'Penk', 'Rasgrp2', 'Myt1l', 'Ca12')  ## HD genes in Khun 2007
(gene = Hmisc::capitalize(tolower(c('Becn1', 'ATG5', 'ATG7', 'SQSTM1'))))  ## autophage genes

# another list from lilah, autophage genes
(genes = Hmisc::capitalize(tolower(c('RSL1D1',
                                     'RAB23',
                                     'Dido1', #(autophagy and apoptosis)
                                     'Fbxo45')))) #(this one relates to ubiquitination))

gene = c('Bdnf')

gene =c('Bcl11b', '')  ## HD genes interact with mHTT or HTT

gene =c('Crhbp')
gene =c('Ccl3', 'Clec7a', 'Treml2')




##HD
gene=c('Adam10')
gene=c('Cox14')
gene=c('Coa3', 'Ndufa12','Ndufs3', 'Mrpl41')  # my top complex I, IV and other mitochondra genes in early HD

gene=c('Sdha','Sdhb','Sdhc','Sdhd') ## complex II genes



## HD genetic modifier listed in HD primer
gene=Hmisc::capitalize(tolower(c('ADORA2A','ATG7','CNR1','GRIK2','GRIN2A','GRIN2B',
                                 'HAP1','PPARGC1A','MAP2K6','MAP3K5','NPY','NPY2R',
                                 'OGG1','PEX7','TP53','UCHL1')))
## replicated in other studies
gene=Hmisc::capitalize(tolower(c('ADORA2A','ATG7','GRIK2','GRIN2A','GRIN2B',
                                 'HAP1','PPARGC1A', 'TCERG1')))

gene=c('CNR1', 'ADORA2A','DRD2')
gene=Hmisc::capitalize(tolower(gene))



gene=c('E2f2', 'Mlh1', 'Sp1') # genetic modifiers

gene= c('Brn-2', 'Bcl11b')
gene =c('Nfya','Nfyc')  # regualte gene expression up in HD late
gene =c('Nr2b')


diseases <- c('AD', 'HD')
diseases <- c('HD')
phases <- c('early', 'late')

df1 <- all_ranks[which(all_ranks$geneSymbol %in% gene & all_ranks$disease %in% diseases & all_ranks$phase %in% phases),]
df2 <- all_ranks_adj[which(all_ranks_adj$geneSymbol %in% gene& all_ranks_adj$disease %in% diseases& all_ranks_adj$phase %in% phases),]
#' Clu is a astrocyte marker, down in late AD HD,
#' no change in Picalm, no Cr1 genes (complement component )


#+++++++++++++++++++
## get which DE genes are also cell type markers before and after correction
#+++++++++++++++++++
## thesis tables /results/ND_results/DE_genes_markers/'
source('./thesis_stuff/check_genes_helpers.R')

threshold = 50
disease_ls=c('AD', 'HD')
phase_ls =c('early','late')
outdir <- '../../results/ND_results/DE_genes_markers/'
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, Sys.Date(), "_DE_markers_top",threshold,".tsv")
if(file.exists(f)){file.remove(f)}
for(disease in disease_ls){
    for(phase in phase_ls){
        print(paste0(disease, '_', phase))
        df_m <- getCellMarkers(disease, phase, threshold)
        writeTable(df_m, f_out = f, append=T)
    }
}

fdr=0.05
disease_ls=c('AD', 'HD')
phase_ls =c('early','late')
outdir <- '../../results/ND_results/DE_genes_markers/'
dir.create(outdir,recursive = T, showWarnings = F)

f <- paste0(outdir, Sys.Date(), "_DE_markers_fdr",fdr,".tsv")
if(file.exists(f)){file.remove(f)}
for(disease in disease_ls){
    for(phase in phase_ls){
        print(paste0(disease, '_', phase))
        df_m <- getCellMarkers(disease, phase, fdr=0.05)
        writeTable(df_m, f_out = f, append=T)
    }
}





#+++++++++++++++++++
# overlap of genes early vs. late
#+++++++++++++++++++
df1= all_ranks_adj
disease1 <- 'AD'
phase1 <- 'early'
rank1 <- 'up_jack'

df2= all_ranks_adj
disease2 <- 'AD'
phase2 <- 'late'
rank2 <- 'up_jack'
threshold <- 20

overlapGenes(df1, disease1, phase1, rank1, df2, disease2, phase2, rank2, threshold)


#+++++++++++++++++++
# overlap of same phase, corrected or not
#+++++++++++++++++++
df1= all_ranks_adj
df2= all_ranks
threshold <- 20

for(disease in c('AD', 'HD')){
    cat('\n######\ncompare genes before and after cell type correction\n######\n')
    for(phase in c('early', 'late')){
        for (reg in c('up_jack', 'down_jack')){
            disease1 <- disease
            phase1 <- phase
            rank1 <- reg
            disease2 <- disease
            phase2 <- phase
            rank2 <- reg
            tmp <- overlapGenes(df1, disease1, phase1, rank1, df2, disease2, phase2, rank2, threshold)
            cat(paste0('\n\n',disease, ', ', phase, ', ', reg, ': overlap ', length(tmp),'\n'))
            cat(tmp)
        }
    }
}


#+++++++++++++++++++
# overlap of early and late comp
#+++++++++++++++++++

df1= all_ranks
df2= all_ranks
threshold <- 20

for(disease in c('AD', 'HD')){
    cat('\n######\ncompare genes before cell type correction: overlap genes of early and late\n######\n')
        for (reg in c('up_jack', 'down_jack')){
            disease1 <- disease
            phase1 <- 'early'
            rank1 <- reg
            disease2 <- disease
            phase2 <- 'late'
            rank2 <- reg
            tmp <- overlapGenes(df1, disease1, phase1, rank1, df2, disease2, phase2, rank2, threshold)
            cat(paste0('\n\n',disease, ', ', reg, ': overlap ', length(tmp),'\n'))
            cat(tmp)
        }

}


df1= all_ranks_adj
df2= all_ranks_adj
threshold <- 40

for(disease in c('AD', 'HD')){
    cat('\n######\ncompare genes after cell type correction: overlap genes of early and late\n######\n')
    for (reg in c('up_jack', 'down_jack')){
        disease1 <- disease
        phase1 <- 'early'
        rank1 <- reg
        disease2 <- disease
        phase2 <- 'late'
        rank2 <- reg
        tmp <- overlapGenes(df1, disease1, phase1, rank1, df2, disease2, phase2, rank2, threshold)
        cat(paste0('\n\n',disease, ', ', reg, ': overlap ', length(tmp),'\n'))
        cat(tmp)
    }
    
}

#+++++++++++++++++++
# overlap of AD and HD (same phase)
#+++++++++++++++++++
df1= all_ranks_adj
df2= all_ranks_adj
threshold <- 200

for(disease in c('AD', 'HD')){
    cat('\n######\ncompare genes after cell type correction: overlap genes of early and late\n######\n')
    for (reg in c('up_jack', 'down_jack')){
        for(phase1 in c('early', 'late')){
            disease1 <- 'AD'
            rank1 <- reg
            disease2 <- 'HD'
            rank2 <- reg
            tmp <- overlapGenes(df1, disease1, phase1, rank1, df2, disease2, phase1, rank2, threshold)
            cat(paste0('\n\n',disease1, ', ', disease2,', ', reg,', ',phase1, ': overlap ', length(tmp),'\n'))
            cat(tmp)
        }
        
    }
}

#+++++++++++++++++++
# overlap of AD and HD paths (adj cell)
#+++++++++++++++++++

fdr = 0.05
df=GO_adj

for (reg in c('up', 'down')){
    for(phase in c('early', 'late')){
df1 <- as.character(df$Name[which(df$disease == 'AD' & df$phase == phase &df$CorrectedPvalue <= fdr &df$regulation ==reg)])
df2 <- as.character(df$Name[which(df$disease == 'HD' & df$phase == phase &df$CorrectedPvalue <= fdr &df$regulation ==reg)])
tmp <- intersect(df1,df2)


cat(paste0('\n\n', reg,', ',phase, ': overlap ', length(tmp),'\n'))
cat(paste0(tmp, collapse = '\n'))
}}

# 
# ## print the top genes under the treshold
# tmp_all <- NULL
# check_list=c('AD_early', 'AD_late', 'HD_early', 'HD_late')
# for(j in 1:length(check_list)){
#     i = check_list[j]
#     assign('df', eval(parse(text=paste0(i, '_mm_variable'))))
#     
#     (x <- length(which(df$pvalue_adj<= threshold)))
#     
#     down_index = which(df$pvalue_adj<= threshold & df$Estimate <0)
#     up_index = which(df$pvalue_adj<= threshold & df$Estimate >0)
#     
#     down_jack = which(df$down_final_jackknife_rank<= length(down_index) & df$pvalue_adj<= threshold)
#     up_jack = which(df$up_final_jackknife_rank<= length(up_index) & df$pvalue_adj<= threshold)
#     
#     
#     y1 = orderCol(df[up_jack, ], 'up_final_jackknife_rank')
#     y2 = orderCol(df[down_jack, ], 'down_final_jackknife_rank')
#     y3 <- list(up = y1$geneSymbol[1: min(nrow(y1), 20)],down = y2$geneSymbol[1: min(nrow(y2), 20)])
#     print(i)
#     print(y3)
# }

###########################
## for HD, look at M2 and M7 genes
###########################
filterDP <- function(df, disease, phase){
    return(df[which(df$disease == disease & df$phase == phase), ]%>%droplevels())
}

df_tmp <- all_ranks_adj

df <- read.delim('../../doc/langfelder_module_M2.tsv', comment.char = '#')
tmp_e_m2 <- noWarnings(left_join(df, filterDP(df_tmp, 'HD', 'early')))

tmp_l_m2 <- noWarnings(left_join(df, filterDP(df_tmp, 'HD', 'late')))

df <- read.delim('../../doc/langfelder_module_M7.tsv', comment.char = '#')
tmp_e_m7 <- noWarnings(left_join(df, filterDP(df_tmp, 'HD', 'early')))

tmp_l_m7 <- noWarnings(left_join(df, filterDP(df_tmp, 'HD', 'late')))
















#+++++++++++++++++++
#' plot the before and after ranking differences, top 50 up and down before 
#' MGP correction are highlighted


#+++++++++++++++++++
f_out <- paste0('../../results/ND_results/figures/before_after_MGP_ranks_', Sys.Date(), '.png')

threshold = 50


df <- all_ranks[, c(1:5,8)]
colnames(df)[c(4:6)] <- paste0(colnames(df)[c(4:6)], '_before')

df_b <- df

df <- all_ranks_adj[, c(1:5, 8)]
colnames(df)[c(4:6)] <- paste0(colnames(df)[c(4:6)], '_after')

df <- na.omit(dplyr::left_join(df_b, df))


## mark the top ranked genes (up and down in separate colomns)
df$top_up_before <- ""
index <- which(df$up_jack_before <= threshold)
df$top_up_before[index] <- 'Top_ranked'
df$top_up_before <- as.factor(df$top_up_before)

df$top_down_before <- ""
index <- which(df$down_jack_before <= threshold)
df$top_down_before[index] <- 'Top_ranked'
df$top_down_before <- factor(df$top_down_before,levels=c('Top_ranked', '') )

## get correlation before and after MGP
disease_ls <- c('AD', 'HD')

phase_ls <- c('early', 'late')
for(disease in disease_ls){
    for(phase in phase_ls){
        df_t <- df[which(df$disease == disease & df$phase == phase), ]
        x <- cor(df_t$down_jack_before,df_t$down_jack_after, method = 'spearman')
        print(paste0(disease, ': ', phase, ': ', x))
    }
}


## plot

reg='up'

x=paste0(reg, '_jack_after')
y=paste0(reg, '_jack_before')
rank_col=paste0('top_',reg,'_before')
df_top_up <- df[which(df$top_up_before == 'Top_ranked') , ]
df_top_down <- df[which(df$top_down_before == 'Top_ranked') , ]

getRankThreshold <- function(phase_ls, disease_ls, df){
    # get the rank to separate FC>0 and FC <0 for each disease phase
    df_all <- NULL
    for(disease in disease_ls){
        for(phase in phase_ls){
            df_tmp <- df[which(df$phase == phase & df$disease == disease), ]
            
            # what's the rank separate the FC>0 beforeMGP
            (before_threshold <- max(df_tmp$up_jack_before[which(df_tmp$Estimate_before >0)]))
            # what's the rank separate the FC>0 afterMGP
            (after_threshold <- max(df_tmp$up_jack_after[which(df_tmp$Estimate_after >0)]))
            
            df_tmp <- data.frame(phase = phase, disease = disease, before_threshold=before_threshold,
                                 after_threshold=after_threshold)
            
            if(is.null(df_all)){
                df_all <- df_tmp
            }else{
                df_all <- rbind(df_all,df_tmp)
            }
        }
    }
    return(df_all)
}

df_t <- getRankThreshold(phase_ls=c('early', 'late'), disease_ls=c('AD', 'HD'), df)


#### plot
one_plot_font_size=14
top_dot_size =2

p <- ggplot(df, aes_string(x = x, y = y))+
    geom_point(alpha = 0.05, size =0.5)+
    theme_bw()

p <- p+ 
    geom_point(data=df_top_up, aes_string(x = x, y = y), colour='red', size=top_dot_size, alpha=0.5)+ # top up ranked
    geom_point(data=df_top_down, aes_string(x = x, y = y), colour='blue', size=top_dot_size, alpha=0.5)+
    facet_grid(phase ~ disease, scales = 'free') +
    ylab('Rank (before MGP correction)') +
    xlab('Rank (after MGP correction)') +
    scale_x_continuous(breaks = c(0, 5000,10000)) +
    scale_y_continuous(breaks = c(0, 5000,10000))
    # scale_x_log10() +
    # scale_y_log10()

p <- p+ theme(
         axis.title.x=element_text(size = one_plot_font_size),
         axis.title.y=element_text(size = one_plot_font_size),
         axis.text.x = element_text(size = 12), 
         axis.text.y=element_text(size = 12),
         strip.text.x = element_text(size = one_plot_font_size),
         strip.text.y = element_text(size = one_plot_font_size))+
    geom_vline(aes(xintercept = before_threshold), df_t)+
geom_hline(aes(yintercept = after_threshold), df_t)
p
ggsave(filename = f_out, plot=p, width = 6.45, height =4.57, units = "in")  
