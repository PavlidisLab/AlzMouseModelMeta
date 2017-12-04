#' 2016-07-27
### run for gene with NA
# done: intercept for all disease, all phasese
# runing: slope NA for AD, PD, late

# need: gene with no NA
#' slope PD, HD, early (ready to go)
#' slope HD, late (runing in otto)
#' slope AD, early (queue)

home_dir <- '/home/bzhuang/ND_project_combined/'

cat("
    #***************************************************************************#
    # PART 6B.2 get which gene with NA value
    #***************************************************************************#\n")
## find which gene are with NAs and copy df_all to expression.Rdata for MM
## this will avoid gene MM values that have been calculated before, so when run the MM, only the genes that
## has not done before in part6

 
rm(list=setdiff(ls(),'home_dir'))

model_ls = c('random_slope')
disease_ls <- c('AD', 'PD')
phase_ls <- c('late')

for(disease in disease_ls){ ## loop1 for disease
    source('config/config_wrappers.R')
    for(model in model_ls){ ## loop2 for models
        data_dir <- paste0(disease_dir, 'mixed_model/', model, '_include_NA/')  ## the one with NA
        out_dir <- data_dir
        input_dir <- paste0(disease_dir,'/mixed_model/', model,'/') ## where the original mixed model 
        
        for (phase in phase_ls){#loop 3 for phase
            ## load previous result and expression with NA, save the df_all to the expression Rdata for model
            (df_all_data <- paste0(input_dir, phase, '/mixed_model_results.Rdata'))
            load(df_all_data)
            (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))  ## get the one with NA
            load(exprdata)
            
            ## what are the genes in with NA that are not in previouse one?
            genes_done <- rownames(df_all)
            msg <- paste0('# ', Sys.Date(), '#\n Input previous results with NA omitted: ', df_all_data,
                          '#\n Input of current expression list : ', exprdata,
                          '\n# genes that are done')
            writeTable(df= data.frame(genes_done=genes_done), f_out = paste0(data_dir, phase, '/genes_done.tsv'),
                       msg = )
            save(array_dat_not_qn, array_dat, array_design, phase, msg, plot_dir,all_array_dat_bf_aggregation, all_df_gene_match,df_all,
                 file = exprdata)
            print(exprdata)
        }
    }
}




# disease_ls <- c('AD', 'PD', 'HD')
# phase_ls <- c('early', 'late')
full_report =T
NA_filter='0.3'
width = 1200
height = 1200


for(disease in disease_ls){## loop1 for disease
    for(model in model_ls){## loop2 for models
        cat("
        #***************************************************************************#
            # PART 6.1 get MM results (random slope and random intercept)
            #***************************************************************************#\n")
        source('config/config_wrappers.R')
        source('mixed_models/mixed_model.R')
        md_info <- paste0(disease_dir, 'config_files/dataset_info_mixed_model.tsv')  #dataset label, phase label, order of the datasets etc.
        datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
        (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
        
        
        #*********************#
        ### run mixed models: random intercept and random slope
        #*********************#
        data_dir <- paste0(disease_dir,'/mixed_model/', model,'_include_NA/') ## where the mixed model with NA
        (out_dir <- paste0(disease_dir, 'mixed_model/', model, '_include_NA/')) ## where the mixed model with NA
        for (phase in phase_ls){ ## loop3 for phases
            ## get genes to rm (which is already done)
            gene_df <- read.delim(paste0(data_dir, phase, '/genes_done.tsv'), comment.char = '#')
            rm_genes <- as.character(gene_df$genes_done)
            ##
            (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
            ### get the expression data and MM results
            x <- mixedModelAll(phase =phase, out_dir =out_dir, 
                               exprdata = exprdata, model=model,
                               full_report = full_report,
                               to_plot= F,rm_genes = rm_genes)   ## must add rm_genes
            assign(x = paste0(model, "_",phase,"_",disease), value = x)
        }
    }

}


