#' 2016-07-11
#' summary, run this script
#' to get 3 files, summary of disease, study and platform inputs
#' 
#' 
#home_dir= 'C:/Users/User/Documents/lab_results/ND_project_combined/'
home_dir ='/home/bzhuang/ND_project_combined/'
rm(list=setdiff(ls(),'home_dir'))
 
disease_ls =c('AD', 'HD', 'PD')

keyword <- 'all_genotypes'
pf_info <- paste0(home_dir, "/git/doc/ND_files/platform_names.tsv")
f_out <- paste0(home_dir, "/git/results/ND_results/ND_summary/")


library('dplyr')
library(HelperFunctions)

######

df_pfinfo <- read.delim(pf_info, comment.char = '#')


#########################################################
## functions
#########################################################

summarizeOneDisease <- function(disease, df_pfinfo){
    print(disease)
    #' input a limma summary table: input_data_summary.tsv
    #' calculate 1. platform info, info per study, and info per disease
    #' output 3 dfs
    
    source('config/config_wrappers.R')
    result_dir <- paste0(disease_dir, 'results/limma_DE_summary/')
    ## get the most result result folder and sub folder with the key word
    (f_dir <- max(grep('2016', list.dirs(result_dir,recursive = F), value = T)))
    (f_dir <- grep(keyword, list.dirs(f_dir), value = T))
    

    (input_f <- grep('input_data', list.files(f_dir, full.names = T), value = T))
    print(input_f)
    if(length(input_f) !=1){
        stop('no input_data_summary.tsv file found in ',f_dir)
    }
    
    df_input <- read.delim(input_f, comment.char = '#')
    ## add study
    df_input$study <- gsub('\\.[0-9]','', df_input$dataset)
    
    ## remove rows with phase =NA
    index <- which(is.na(df_input$phase))
    if(length(index) >=1){
        df_input <- df_input[-index, ]%>% droplevels()
    }
    
    
    #######
    # get platforms (count of genes and probes)
    #######
    df_pf <- rmDup(df_input[, c("platform","probes","genes")]%>% droplevels())
    
    ## get the full name
    df_pf <- noWarnings(left_join(df_pf, df_pfinfo))
    
    
    #######
    # get a disease overview for each dataset:
    # number of case and control samples, number of disease genotypes
    #######
    df_input$agg_col <- paste(df_input$dataset,df_input$Timepoint,df_input$phase, sep = '_')
    
    ## some wt are shared between genotypes if it's the same timepoint
    df_control <- rmDup(noWarnings(left_join(df_input[, c('dataset', 'Timepoint', 'phase', 'agg_col')], 
                                             aggregate(control_n ~ agg_col, max, data = df_input))))
    ## case number count is per genotype
    df_case <- aggregate(case_n ~ agg_col, sum, data = df_input)
    
    ## get number of genotypes
    ## case number count is per genotype
    df_geno_c <- aggregate(Genotype ~ agg_col, length, data = df_input)
    colnames(df_geno_c) <- c("agg_col","Genotype_n")
    df_geno_ls <- aggregate(Genotype ~ agg_col, function(x)paste0(x, collapse = '; '), data = df_input)
    
    df_study <- noWarnings(Reduce(left_join, list(df_control,df_case,df_geno_c, df_geno_ls, 
                                                  rmDup(df_input[, c('dataset', 'study')]))))
    df_study$Disease <- disease
    
    #######
    # get a disease overview for the disease,by phase
    #######
    
    df_1 <- aggregate(control_n ~ phase,sum, data=df_study)
    df_2 <- aggregate(case_n ~ phase,sum, data=df_study)
    # get genotypes and count
    df_3 <- aggregate(Genotype ~ phase, function(x)unique(paste0(x)), data = df_input)
    df_3$Genotype_n <- ''
    df_3$Genotype_list <- ''
    for(i in 1: length(df_3$Genotype)){
        x <- sort(unlist(df_3$Genotype[i]))
        df_3$Genotype_n[i] <- length(x)
        df_3$Genotype_list[i] <- paste0(x, collapse = ', ')
    }
    df_3 <- df_3[, -2] ## remove tmp genotype list
    
    # get studys and count
    df_4 <- aggregate(study ~ phase, function(x)unique(paste0(x)), data = df_input)
    df_4$study_n <- ''
    df_4$study_list <- ''
    for(i in 1: length(df_4$study)){
        x <- sort(unlist(df_4$study[i]))
        df_4$study_n[i] <- length(x)
        df_4$study_list[i] <- paste0(x, collapse = ', ')
    }
    df_4 <- df_4[, -2] ## remove tmp study list
    
    # get platforms and count
    df_5 <- aggregate(platform ~ phase, function(x)unique(paste0(x)), data = df_input)
    df_5$platform_n <- ''
    df_5$platform_list <- ''
    for(i in 1: length(df_5$platform)){
        x <- sort(unlist(df_5$platform[i]))
        df_5$platform_n[i] <- length(x)
        df_5$platform_list[i] <- paste0(x, collapse = ', ')
    }
    df_5 <- df_5[, -2] ## remove tmp platform list
    
    df_disease <- noWarnings(Reduce(left_join, list(df_1,df_2,df_3, df_4, df_5)))
    df_disease$disease <- disease
    
    
#     ## get the outliers
#     (rm_f <- grep('removed_samples', list.files(f_dir, full.names = T), value = T))
#     if(length(rm_f) !=1){
#        print(paste0('no removed_samples.tsv file found in ',f_dir))
#     }else{
#         df_rm <- read.delim(rm_f, comment.char = '#')
#         df_outlier <- filterContain(df_rm, 'category','outlier')
#         df_subset <- filterContain(df_rm, 'category','removed by subset')
#     }
    

    
    
    ## log and return values
    msg <- paste0('# ', input_f)
    returnlist <- list(df_pf, df_study, df_disease, msg)
    return(returnlist)
    
}


#########################################################
## loop to get summary of platform, studies, and disease input
#########################################################
dir.create(f_out, recursive = T, showWarnings = F)
df_pf_all <- NULL
df_study_all <- NULL
df_disease_all <- NULL
msg_all <- paste0('# ', Sys.Date())

for(disease in disease_ls){
    print(disease)
    x <- summarizeOneDisease(disease, df_pfinfo)
    df_pf <- x[[1]]
    df_study <- x[[2]]
    df_disease <- x[[3]]
    msg <- x[[4]]
    
    df_pf_all <- rbind(df_pf_all, df_pf)
    df_study_all <- rbind(df_study_all, df_study)
    df_disease_all <- rbind(df_disease_all,df_disease)
    msg_all <- paste0(msg_all, '\n', msg)
}

## remove duplicated platform
index <- which(duplicated(df_pf_all$platform))
if(length(index) >0){
    df_pf_all <- df_pf_all[-index, ]%>% droplevels()
}
## sort results:
df_study_all <- df_study_all[with(df_study_all, order(Disease, phase)), ]

## update all column names to first cap
colnames(df_pf_all) <- Hmisc::capitalize(colnames(df_pf_all))
colnames(df_study_all) <- Hmisc::capitalize(colnames(df_study_all))
colnames(df_disease_all) <- Hmisc::capitalize(colnames(df_disease_all))


## write to table
f_pf <- paste0(f_out, 'methods_platform_info', Sys.Date(), '.tsv')
writeTable(df_pf_all, f_out =f_pf, msg=msg_all)
f_study <- paste0(f_out, 'methods_study_list_info', Sys.Date(), '.tsv')
writeTable(df_study_all, f_out =f_study, msg=msg_all)
f_d <- paste0(f_out, 'methods_disease_input_info', Sys.Date(), '.tsv')
writeTable(df_disease_all, f_out =f_d, msg=msg_all)


#########################################################
## get outliers summary -- see manual input
#########################################################