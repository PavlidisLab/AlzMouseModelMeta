#2016-11-25
# updated 2017-02-01, 2017-02-20- add to select which ermineJ folder
#' summarize all the Go enrichment analysis data (add GO aspect to the results)
#' to plot GO terms
#' r object saved for all process results for corrected and non corrected enrichment results
#' gather all the GO enrichment analysis results- from all process
#' paste0('../../results/emineJ_', Sys.Date(),'.Rdata')


 


########################
## part 1 get the GO aspects
########################
## read the go terms downloaded from MGI
library(HelperFunctions)
library(dplyr)
(f <- '../configs/GO_annotation/gene_association.mgi.gz')
    
df <-read.delim(f, comment.char = '!',header = F) 

df <- df[, c(3,5,9,12,14,15)] %>%droplevels()

colnames(df) <- c('geneSymbol', 'ID', 'Aspect', 'gene_type', 'updated', 'source')

## biological process terms:
#go_bio <- df[which(df$Aspect == 'P'), ]%>%droplevels()

go_aspect <- rmDup(df[, c('ID', 'Aspect')])




##################
############## helper functions
##################

plotGo <- function(go_fdr, outdir, start_i = 1){
    ## make the tree plot of GO term relationships
    library(RamiGO)
    for(i in start_i:length(go_fdr)){
        print(i)
        goIDs <- unlist(go_fdr[[i]]$go_terms)
        print
        color <- c("lightblue")
        ## get results
        file_name <- paste0(outdir, Sys.Date(), '_', go_fdr[[i]]$disease, '_', go_fdr[[i]]$phase,'_',go_fdr[[i]]$regulation,'.svg')
        
        #pp <- getAmigoTree(goIDs=goIDs,color=color,filename=file_name,picType="svg")
        file_name <- paste0(outdir, Sys.Date(), '_', go_fdr[[i]]$disease, '_', go_fdr[[i]]$phase,'_',go_fdr[[i]]$regulation,'.png')
        pp <- getAmigoTree(goIDs=goIDs,color=color,filename=file_name)
        print(file_name)
    }
}


goTermsAspect <- function(go_aspect,disease_ls,phase_ls,regulation_ls,adj_cell_pop,keyword,fdr,filter = F,
                          folder_pre){
    #' keyword ='_all_processes' or ''     ## folder name contains
    #' adj_cell_pop ='_adj_cell_pop' or ''  ## folder name contains
    #' folder_pre = 'ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm', the ermineJ folder and the specific folder
    #' the script will select for the most recent dated folder
    #' filter=T: select sig pathways < fdr
    df_all <- NULL
    go_fdr <- vector('list')
    counter <- 1
    
    f_ls_all <- NULL
    for(disease in disease_ls){
        print(disease)
        folder <- paste0(home_dir,'/',disease, '_mouse_model_project/', folder_pre, adj_cell_pop, keyword,'/')
        (folder <- max(list.dirs(folder)))
        print(paste0('folder ', folder))
        
        (f_ls <- list.files(folder,pattern = '_regulation_mixed_model.tsv'))

        for(phase in phase_ls){# loop phase
            for(regulation in regulation_ls){# loop reg
                (f=paste0(folder,'/', grep(paste0(phase, '_', regulation), f_ls, value = T)))
                f_ls_all <- c(f_ls_all, f)
                print(f)
                df <- read.delim(f, comment.char = '#')
                dft <- data.frame(disease = disease, regulation =regulation, phase =phase, adj_cell_pop =adj_cell_pop)
                df <- cbind(dft, df)
                
                ## which go terms are under fdr
                go_terms <- as.character(df$ID[which(df$CorrectedPvalue <= fdr)])
                if(length(go_terms) >0){
                    x <- list(disease, phase, regulation, go_terms)
                    names(x) <- c('disease', 'phase', 'regulation', 'go_terms')
                    go_fdr[counter] <- list(x)
                    counter <- counter+1
                }
                
                if(is.null(df_all)){
                    df_all <- df
                }else{
                    df_all <- rbind(df_all,df)
                }
            }# loop reg
        }# loop phase
        print(disease)
    }# loop disease
    
    
#     # get AD_HD top 50
#     ## if combined
#     if(combined){
#         phase_ls = c('early', 'late')
#         regulation_ls = c('up', 'down')
#         folder <- paste0(home_dir,'ND_results/ermineJ/results_jackknife', adj_cell_pop, keyword, '/AD_HD/')
#         (folder <- max(list.dirs(folder)))
#         print(folder)
#         (f_ls <- list.files(folder,pattern = '_regulation_mixed_model.tsv'))
#         
#         for(phase in phase_ls){
#             for(regulation in regulation_ls){
#                 (f=paste0(folder,'/', grep(paste0(phase, '_', regulation), f_ls, value = T)))
#                 print(f)
#                 f_ls_all <- c(f_ls_all, f)
#                 df <- read.delim(f, comment.char = '#')
#                 dft <- data.frame(disease = 'combined', regulation =regulation, phase =phase, adj_cell_pop =adj_cell_pop)
#                 df <- cbind(dft, df)
#                 
#                 ## which go terms are under fdr
#                 go_terms <- as.character(df$ID[which(df$CorrectedPvalue <= fdr)])
#                 if(length(go_terms) >0){
#                     x <- list('combined', phase, regulation, go_terms)
#                     names(x) <- c('disease', 'phase', 'regulation', 'go_terms')
#                     go_fdr[counter] <- list(x)
#                     counter <- counter+1
#                 }
#                 
#                 
#                 if(is.null(df_all)){
#                     df_all <- df
#                 }else{
#                     df_all <- rbind(df_all,df)
#                 }
#             }
#         } 
#     }
   
    
    
    ## get the aspect type
    ej <- noWarnings(left_join(df_all, go_aspect)) 
    if(filter){ ## filter sig pathways < fdr
        ej <- ej[which(ej$CorrectedPvalue <=fdr), ]%>%droplevels()
    }
    
    
    
    
    ## process the columns
    tmpFun <- function(df){
        df <- df[which(df$disease != 'combined'), ]%>%droplevels() # rm combined
        df$P_adj <- round(df$CorrectedPvalue, digits = 3)
        df$MP_adj <- round(df$CorrectedMFPvalue, digits = 3)
        colnames(df)
        df <- df[, colnames(df)[c(1:5, 19:20, 8:9, 18, 15,10, 12:14, 7,16:17,6)]]
        return(df)
    }
    ej <- tmpFun(ej)
    
    return(list(ej, go_fdr, f_ls_all))
}




########################
######## run functions [adj- all processes, filter FDR<0.05]
########################
adj_cell_pop ='_adj_cell_pop'  ## folder name contains
keyword ='_all_processes'      ## folder name contains


fdr = 0.05

## get AD, HD  top 50
phase_ls = c('early', 'late')
regulation_ls = c('up', 'down')

tmp <- goTermsAspect(go_aspect,disease_ls,phase_ls,regulation_ls,adj_cell_pop,keyword,fdr, filter = T, folder_pre = folder_pre)
GO_adj <- tmp[[1]]    
go_fdr_adj <- tmp[[2]]
(outdir <- paste0(go_plot_out, 'GO_figures',adj_cell_pop,keyword, '/', Sys.Date(), '/'))
dir.create(outdir, showWarnings = F, recursive = T)
plotGo(go_fdr_adj, outdir, start_i = 1)


########################
## for AD HD [not adj -all processes, filter FDR<0.05]
########################
adj_cell_pop =''  ## folder name contains
keyword ='_all_processes'      ## folder name contains

fdr = 0.05

## get top 50
phase_ls = c('early', 'late')
regulation_ls = c('up', 'down')

tmp <- goTermsAspect(go_aspect,disease_ls,phase_ls,regulation_ls,adj_cell_pop,keyword,fdr, filter = T, folder_pre = folder_pre)
GO_non_adj <- tmp[[1]]    
go_fdr_non_adj <- tmp[[2]]
(outdir <- paste0(go_plot_out, 'GO_figures_',adj_cell_pop,keyword, '/', Sys.Date(), '/'))
dir.create(outdir, showWarnings = F, recursive = T)
plotGo(go_fdr_non_adj, outdir, start_i = 1)

# save(GO_non_adj, GO_adj, go_fdr_non_adj, go_fdr_adj, file = f_result_out)


########################
## for AD HD combined results [not adj -all processes, filter FDR<0.05]
########################
# adj_cell_pop_ls =''  ## folder name contains
# keyword_ls =''      ## folder name contains
# fdr = 0.05
# for(adj_cell_pop in adj_cell_pop_ls){
#     for(keyword in keyword_ls){
#         tmp <- goTermsAspect(go_aspect,disease_ls,phase_ls,regulation_ls,adj_cell_pop,keyword,fdr, combined = F)
#         #ej_all_process <- tmp[[1]]    
#         go_fdr <- tmp[[2]]
#         f_ls_all <- tmp[[3]]
#         msg <- paste0(sys.date(), '\nFDR = ', fdr, '\nINPUT files:\n\n', paste0(f_ls_all, collapse = '\n'))
#         (outdir <- paste0(go_plot_out, 'GO_figures_disease_combined/',adj_cell_pop,keyword, '/', Sys.Date(), '/'))
#         dir.create(outdir, showWarnings = F, recursive = T)
#         writeTable(df=NULL,msg = msg, f_out = paste0(outdir, 'inputs.txt'))
#          plotGo(go_fdr_all_process, outdir, start_i = 1) 
#     }
# }

########################
## for AD HD combined results [adj -all processes, filter FDR<0.05]
########################
# adj_cell_pop_ls ='_adj_cell_pop'  ## folder name contains
# keyword_ls =c('','_all_processes' )
# fdr = 0.05
# for(adj_cell_pop in adj_cell_pop_ls){
#     for(keyword in keyword_ls){
#         tmp <- goTermsAspect(go_aspect,disease_ls,phase_ls,regulation_ls,adj_cell_pop,keyword,fdr, combined = T)
#         #ej_all_process <- tmp[[1]]    
#         go_fdr <- tmp[[2]]
#         f_ls_all <- tmp[[3]]
#         msg <- paste0('FDR = ', fdr, '\nINPUT files:\n\n', paste0(f_ls_all, collapse = '\n'))
#         (outdir <- paste0(go_plot_out, 'GO_figures_disease_combined/',adj_cell_pop,keyword, '/', Sys.Date(), '/'))
#         dir.create(outdir, showWarnings = F, recursive = T)
#         writeTable(df=NULL,msg = msg, f_out = paste0(outdir, 'inputs.txt'))
#         plotGo(go_fdr_all_process, outdir, start_i = 1) 
#     }
# }






