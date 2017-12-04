#' 2016-07-20
#' explore genes, get results (meta and mixed model) and plots for the gene
library(dplyr)
library(HelperFunctions)

## loop for each diseae

#'@example
#'
# home_dir= '/home/bzhuang/'
# # compare mixed model results to metaanalysis results
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease = 'HD'
# source('config/config_wrappers.R')
# 
# 
# (combined_dir = paste0('/home/bzhuang/ND_project_combined/', disease, '_mouse_model_project/'))

# 
# source('./result_explore/explore_gene.R')

# gene_ls <- 'Gpx6'
# for(gene in gene_ls){
#     (outdir <- paste0(disease_dir,'gene_results/', gene,'/'))
#     unlink(outdir, recursive = T) ## delete the old file
#     msg <- getMsg(gene, mj_df, msg_df,outdir)
#     cpImg(gene, outdir, disease_dir, combined_dir=combined_dir)
# }


####################
## get the meta, jack values and mixed models
## get meta jack for combined geneotypes
####################

## meta and jack
(jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                          paste0(disease_dir, 'meta_analysis/meta_jack/'),
                          paste0(combined_dir, 'meta_analysis/low_exp_rm/'),
                          paste0(combined_dir, 'meta_analysis/meta_jack/')))
for(jack_meta_folder in jack_meta_folder_ls){
    mj_f_ls <- grep('regulation', list.files(jack_meta_folder, full.names = T), value = T)
}
## mixed models
mm_dir = paste0(disease_dir, 'mixed_model/') ## mixed model dir(parent dir)
## get results
(file_ls <- grep('mixed_model_results.tsv', list.files(mm_dir, recursive = T, full.names = T), value = T))

##combine all files
mj_f_ls <- c(mj_f_ls, file_ls)

# remove archive, archival files
mj_f_ls <- grep('archive|archival', mj_f_ls, value = T, invert = T)

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



####################
## functions
####################

getMsg <- function(gene, mj_df, msg_df,outdir){
    msg_all <- paste0(Sys.Date(), ': ', gene)
    for(i in 1:length(mj_df)){
        x <- filterContain(mj_df[[i]], column = 'geneSymbol', value = gene)
        msg <- msg_df[i]
        if(nrow(x) ==0){
            msg <- paste0('\n',msg, 'NA\n')
        }else{
            for (j in 1:nrow(x)){
                msg <- paste0('\n',msg,paste0(x[j,], collapse = '\t'), '\n')
            }
        }
        msg_all <- paste0(msg_all, msg, '\n')
    }
    dir.create(outdir, recursive = T, showWarnings = F)
    writeTable(df=NULL, f_out= paste0(outdir, '/results_', gene,'.tsv'), msg=msg_all)
    return(msg_all)
}



#' copy figures from significant genes
cpImg <- function(gene, outdir, disease_dir, combined_dir=''){
    print(gene)
    #' combined_dir: if genotypes are combined (e.g. AD, and PD)
    #' outdir: copy the images to the outdir
    mm_dir = paste0(disease_dir, 'mixed_model/') ## mixed model dir(parent dir)
    sig_dir = paste0(disease_dir, 'results/meta_analysis_significant_genes/')
    (sig_ls <- grep(gene, list.files(sig_dir, recursive = T, full.names = T), value = T))
    if(combined_dir !=''){
        sig_combined_dir = paste0(combined_dir, 'results/meta_analysis_significant_genes/')
        (sig_c <- grep(gene, list.files(sig_combined_dir, recursive = T, full.names = T), value = T))
        (sig_c <- gsub(combined_dir, 'combined_', sig_c)) ## label combined dir
    }else{
      sig_c=NULL
    }
    
    img_ls <- c(sig_ls, sig_c, 
                grep(gene, list.files(mm_dir, recursive = T, full.names = T), value = T))
    (to_img_dir <- gsub(paste0(disease_dir,'|results/'), '', img_ls))
    (to_img_dir <- gsub('//|/', '_', to_img_dir))
    (to_img_dir <- paste0(outdir, to_img_dir))
    

    
    file.copy(from = img_ls, to = to_img_dir)
}





## expression




