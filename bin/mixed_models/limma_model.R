#' 2016-10-31
#' input expression data for all samples in the phase (study corrected or not)
#' use limma to correcrt for study, cell types (optional)
#' 
#' 
#' 

# if(disease == 'AD'){
#     cell_types = c("Astrocyte" ,"DentateGranule", 'GabaSSTReln',"Microglia", 'Oligo', 'Pyramidal_Thy1',
#                    "GABAergic",  "Pyramidal")
# }else if(disease =='HD'){
#     cell_types = c("Astrocyte","Cholinergic" ,"Microglia" ,"Spiny", 'Oligo','ForebrainCholin')
# }
# 
# (cell_markers_f=max(list.files(paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/', disease, '/', phase), recursive = T, full.names = T,
#                                pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' )))

library(HelperFunctions)
library(limma)
limmaModel <- function(exprdata,
                       phase,
                       out_dir,
                       cell_markers_f,
                       cell_types = NULL,
                       include_study =F){
    #' cell_markers_f', e.g. cell_markers_f= paste0(home_dir, '/ND_results/cell_population/2016-10-26/cell_population_estimate_2016-10-26.tsv')
    #' cell_types: input cell types for correction
    print(Sys.time())
    
    load(exprdata)
    print(exprdata)
    ## dir to save rdata
    (result_dir <- paste0(out_dir, '/', phase,'/'))
    dir.create(result_dir, recursive = T, showWarnings = F)
    
    ###
    #process input cell propotions
    ###
    if(!is.null(cell_types)){
        cell_markers <-read.delim(cell_markers_f, comment.char = '#')
        
        ## double check sample names of cell type and expression matrix
        if(length(unique(as.character(cell_markers$Sample[which(cell_markers$Sample %in% array_design$Sample)]))) != length(unique(as.character(array_design$Sample)))){
            stop(print('samples names do not match'))
        }else{
            print('samples names match, continue')
        }
        
        ## get only the samples
        cell_markers <- cell_markers[which(cell_markers$Sample %in% array_design$Sample), ] %>%droplevels()
        
        ### get the cells
        cell_types =intersect(cell_types,levels(cell_markers$cell_type))
        print(paste0('Input cell types are ', paste0(cell_types, collapse = ', ')))
        
        
        cell_markers <- cell_markers[which(cell_markers$cell_type %in% cell_types), ]%>%droplevels()
        
        
        #get the PC1 scaled for each cell type (from long to wide table)
        df_markers <- cell_markers[, c('Sample', 'cell_type', 'PC1_scaled')]%>%droplevels()
        df_markers <- rmDup(df_markers)
        df_markers <- reshape2::dcast(df_markers,Sample~cell_type)
        rownames(df_markers) <- df_markers$Sample
        
        ## add cell pop in the design file
        array_design <- noWarnings(left_join(array_design, df_markers))
    }else{
        df_markers =NULL
    }
    
    
    ## make sure WT is the baseline (first level)
    array_design$Disease_stage <- factor(array_design$Disease_stage, levels = c('WT', "Disease"))
    ###########
    ## Limma
    ###########
    factor_ls <- c(cell_types)
    if(include_study){
        factor_ls <- c(factor_ls, 'Study')
    }
    
    (model_formular <- paste0('model.matrix( ~ Disease_stage +', paste0(factor_ls, collapse=' + '), ', data = array_design)' ))
    
    print(paste0('Model: ', model_formular))
    assign('design_matrix', eval(parse(text = model_formular)))
    array_fit <- lmFit(array_dat, design_matrix)
    array_eb_fit <- eBayes(array_fit)
    
    ########################################
    # loop for each contrast coef:
    #' 1. plot p value density, save top tables, save toptables after filter
    #######################################
    
    #*****************************#
    ### get disease coef
    #*****************************#
    target_coef <- 'Disease_stageDisease'
    print(paste0('coef: ', target_coef))
    array_toptable <- topTable(array_eb_fit, 
                               coef = target_coef,
                               number = Inf) # get all the probes
    ## process top table (pvalue for disease stage only)
    toptable_df_all <- cbind(rownames(array_toptable),
                             array_toptable[, c('logFC', 't','P.Value', 'adj.P.Val')])
    
    table_title <- c('geneSymbol','Estimate', 't.value','pvalue', 'adj.P.Val')
    colnames(toptable_df_all) <- table_title
    
    # get the other coefs (not the p values)
    if(length(factor_ls) >0){
        array_toptable <- topTable(array_eb_fit, number = Inf) 
        rm_coef <- c("Disease_stageDisease", "AveExpr", "F","P.Value","adj.P.Val") ## this p value is for the entire model
        array_toptable <- array_toptable[, setdiff(colnames(array_toptable), rm_coef)]
        array_toptable <- cbind(data.frame(geneSymbol =rownames(array_toptable)),
                                array_toptable)
    }
    df_all <- noWarnings(left_join(toptable_df_all, array_toptable))
    (msg <- paste0(msg, '\n ####### ', Sys.Date(), "#######\n Input expression: ", exprdata, "\nmodel: ",model_formular))
    
    save(array_dat,array_design, df_all, msg, df_markers, file=paste0(result_dir, "mixed_model_results.Rdata"))
    writeTable(df=NULL, msg = msg, f_out = paste0(result_dir, 'mixed_model_results_log_', Sys.Date(), '.txt'))
    
    returnlist=list(array_dat,array_design, df_all, msg)
    print(paste0(result_dir, "mixed_model_results.Rdata"))
    return(returnlist)
    
}
