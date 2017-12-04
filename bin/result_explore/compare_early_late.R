#'2016-07-25
#'compare up or down regulated genes between early and late 
#' (meta, jack, and mixed models of the same disease)
#' are they similar?
#'

home_dir= '/home/bzhuang/ND_project_combined/'
 
rm(list=setdiff(ls(),'home_dir'))
disease = 'PD'
source('config/config_wrappers.R')


#### define folders
####################
## get the meta, jack values and mixed models
## get meta jack for combined geneotypes
####################

## meta and jack
(jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                          paste0(disease_dir, 'meta_analysis/meta_jack/')))
for(jack_meta_folder in jack_meta_folder_ls){
    mj_f_ls <- grep('regulation', list.files(jack_meta_folder, full.names = T), value = T)
}
## mixed models
mm_dir = paste0(disease_dir, 'mixed_model/') ## mixed model dir(parent dir)
## get results
(file_ls <- grep('regulation', list.files(mm_dir, recursive = T, full.names = T), value = T))

##combine all files
mj_f_ls <- c(mj_f_ls, file_ls)

## remove archival
mj_f_ls <- grep('archive|archival', mj_f_ls, invert = T, value = T)

## keyword:
regulation ='up'
threshold =100

f=early_ls[i]
phase ='early'
value_col <- c('geneSymbol','up_pval', 'down_pval','adj_combined_max_p', 'Fisher')

(early_ls <- sort(grep(paste0('early.*', regulation), mj_f_ls, value = T)))
(late_ls <- sort(grep(paste0('late.*', regulation), mj_f_ls, value = T)))

processDf <- function(f, phase, value_col){
    df <- read.delim(f, comment.char = '#')
    (need_col <- intersect(value_col, colnames(df)))
    df <- df[, need_col] 
    colnames(df) <- c('geneSymbol',phase)
    df <- eval(parse(text = paste0("df[with(df, order(", phase, ")), ]"))) ## reorder by pvalue
    return(df)
}


for(i in 1:length(early_ls)){
    df_early <- processDf(early_ls[i], phase ='early', value_col)
    df_late <- processDf(late_ls[i], phase ='late', value_col)
    df <- na.omit(noWarnings(left_join(df_early, df_late)))
    cor(df$early, df$late, method = 'spearman')
    genes <- intersect(df_early[1:threshold, 'geneSymbol'], df_late[1:threshold, 'geneSymbol'])
    df_early[which(df_early$geneSymbol %in% genes), ]
    df_late[which(df_late$geneSymbol %in% genes), ]
    
    }




mj_df <- vector()
msg_df <- vector()
## read dfs
for (i in 1:length(mj_f_ls)){
    f <- mj_f_ls[i]
    df <- read.delim(f, comment.char = '#', stringsAsFactors = F)
    msg <- paste0('###\n',f, '\n###\n',paste0(as.character(colnames(df)), collapse = '\t'), '\n')
    mj_df[i] <- list(df)
    msg_df[i] <- msg
}

