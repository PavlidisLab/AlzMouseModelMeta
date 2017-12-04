#2016-11-1
## wrapper for limma model for combined samples (must run mixed models first)
#' results including limma and limma jackknife (stored in mixed_model, and mixed_model_jackknife folders)
source('mixed_models/limma_model.R')
source('helper_functions.R')
limmaModelJack <- function(model_keyword,
                           corrected_exp_input,
                           expr_data_keyword,
                           include_study,
                           disease_ls =c('AD','HD'),
                           model = c('random_intercept'),
                           phase_ls =c('early','late'),
                           model_ls =c('random_intercept'),
                           regulation_ls= c('up', 'down'),
                           xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz',
                           AD_cell_types = c("Astrocyte" ,"DentateGranule", 'GabaSSTReln',"Microglia", 'Oligo', 'Pyramidal_Thy1',
                                             "GABAergic",  "Pyramidal"),
                           HD_cell_types = c("Astrocyte","Cholinergic" ,"Microglia" ,"Spiny", 'Oligo','ForebrainCholin')){
    
    cat("
        #---------------------------------------------------------------------------#
        # PART 8.6: with expression corrected by study, use limma for model
        #---------------------------------------------------------------------------#\n")
    
    #*********************#
    #'with expression corrected by study, use limma for model
    ### 8.6.1 prepare combined expression data
    #' copy the corrected to the new expression
    #' and make a copy in mixed model jackknife
    #*********************#
    # 
    
    # disease_ls <- c('AD','HD')
    # model <- c('random_intercept')
    # phase_ls <- c('early', 'late')
    
    #corrected_exp_input <-'_include_NA_low_exp_rm' 
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    #expr_data_keyword = 'mixed_model_results_exp_corrected.Rdata' ## with type of expression data to start from

    for (disease in disease_ls){## loop 1 for disease
        print(disease)
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        for (phase in phase_ls){ ## loop 2 for phases
            (in_dir <- paste0(disease_dir, 'mixed_model/', model, corrected_exp_input,'/', phase,'/'))
            
            ## input the study corrected rdata and make a new expression.Rdata
            (rfile_in <- grep(expr_data_keyword, list.files(in_dir, recursive=T, full.names=T), value=T))
            load(rfile_in)
            print(paste0(disease, ':', phase))
            print(paste0('input expression : ', rfile_in))
            ######
            ## save in mixed model
            ######
            (out_dir <- paste0(disease_dir, 'mixed_model/', model, model_keyword,'/', phase, '/'))  ## name outdir by model
            dir.create(out_dir, showWarnings = F, recursive = T)
            save(array_dat, array_design, phase, msg, file = paste0(out_dir, 'expression.Rdata'))
            plotFunctions(na.omit(array_dat), out_dir, dataset="", array_design=array_design, prefix="", simple_heatmap = F)
            
            ######
            ## save in mixed model jackknife
            ######
            (out_dir <- paste0(disease_dir, 'mixed_model_jackknife/', model, model_keyword,'/', phase, '/'))  ## name outdir by model
            dir.create(out_dir, showWarnings = F, recursive = T)
            save(array_dat, array_design, phase, msg, file = paste0(out_dir, 'expression.Rdata'))
        }
    }
    
    
    
    #*********************#
    #'with expression corrected by study, use limma for model
    ### 8.6.2 JACKKNIFE: creat expression.Rdata when 1 of the study is removed, in subfolder, run1...runX
    #*********************#
    # see more examples in ('mixed_models/mixed_model.R')
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    
    
    # disease_ls <- c('AD','HD')
    # model <- c('random_intercept')
    # phase_ls <- c('early','late')
    
    
    
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    library(dplyr)
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        print(disease)
        
        for (phase in phase_ls){ ## loop3 by phase
            data_dir <- paste0(disease_dir,'/mixed_model_jackknife/', model,model_keyword, '/', phase, '/') ## where the mixed model input expression
            (exprdata <- paste0(data_dir, '/expression.Rdata'))
            
            load(exprdata)
            study <- levels(array_design$Study)
            for(i in 1: length(study)){
                
                load(exprdata)  ## must reload the array design and array data for every run
                
                (out_dir <- paste0(data_dir, '/run',i,'/'))
                plot_dir <- out_dir
                dir.create(out_dir, showWarnings = F)
                (rm_study <- study[i])
                print(rm_study)
                ## remove 1 study samples
                array_design <- array_design[which(array_design$Study != rm_study), ] %>% droplevels()
                
                print(levels(array_design$Study))
                array_dat <- array_dat[, array_design$Sample]%>% droplevels()
                
                ## check for gender
                df_tmp <- filterContain(array_design, column = 'Gender', value = 'F')
                if(nlevels(df_tmp$Study) == 1){
                    array_design$Gender=NA
                    print('rm gender only due to 1 dataset has female')
                }
                df_tmp <- filterContain(array_design, column = 'Gender', value = 'M')
                if(nlevels(df_tmp$Study) == 1){
                    array_design$Gender=NA
                    print('rm gender due to only 1 dataset has female')
                }
                ##
                rm_msg <- paste0(Sys.Date(), '\n Study removed: ', rm_study,
                                 '\n keep studies (',nlevels(array_design$Study), '): ' ,paste0(levels(array_design$Study), collapse = ', '))
                writeTable(df=NULL, f_out = paste0(out_dir, 'rm_sample.txt'), 
                           msg = rm_msg)
                msg <- paste0(msg, '\n', rm_msg)
                ## save r obj
                save(array_dat, array_design, phase, msg,rm_study,plot_dir, 
                     file = paste0(out_dir, 'expression.Rdata'))
                
                print(out_dir)
            }
        }
        
    }##loop1
    
    #*********************#
    #'with expression corrected by study, use limma for model [not jackknife]
    ### 8.6.3 run mixed models: random intercept, NA removed (all genes): correct for cell types
    #' updated 2016-10-30: use the new population estimation from study corrected expressions of all samples in the phase
    #' the input cell pop file is in /ND_results/cell_population/all_sample_estimation/'
    #*********************#
    # see more examples in ('mixed_models/mixed_model.R')
    # 
    #rm(list=setdiff(ls(),'home_dir'))

    
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    
    
    # disease_ls <- c('AD','HD')
    # model_ls <- c('random_intercept')
    # phase_ls <- c('early', 'late')
    
    
    # include_study = F
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        if(disease == 'AD'){
            cell_types = AD_cell_types
        }else if(disease =='HD'){
            cell_types = HD_cell_types
        }
        
        for(model in model_ls){ ## loop2 by model
            out_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword , '/') ## where the mixed model result
            data_dir <- paste0(disease_dir,'/mixed_model/', model,model_keyword, '/') ## where the mixed model input expression
            
            for (phase in phase_ls){ ## loop3 by phase
                (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
                print(exprdata)
                ## get the cell marker files
                cell_markers_f=max(list.files(paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/', disease, '/', phase), recursive = T, full.names = T,
                                              pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' ))
                print(cell_markers_f)
                ### get the expression data and limma results
                x <- limmaModel(exprdata,
                                phase,
                                out_dir,
                                cell_markers_f,
                                cell_types = cell_types,
                                include_study = include_study)
                assign(x = paste0(model, "_",phase,"_",disease), value = x)
            }
        }## loop2 end
    }##loop1
    
    cat("
        #***************************************************************************#
        #'with expression corrected by study, use limma for model
        # PART 8.6.4 compare the MM results to Fisher results and also prepare the ermineJ files
        # forrandome intercept
        #***************************************************************************#\n")
    
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    source('mixed_models/compare_mm_meta.R')
    
    
    # disease_ls <- c('AD','HD')
    # model_ls <- c('random_intercept')
    # phase_ls <- c('early','late')
    # model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected'
    
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        for (model in model_ls){## loop2 for models
            (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                      paste0(disease_dir, 'meta_analysis/meta_jack/')))
            mm_dir = paste0(disease_dir, 'mixed_model/',model,model_keyword,'/') ## mixed model dir(parent dir)
            mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls=phase_ls, compare_fisher =F)
        }## loop2 end
    }##loop1
    
    
    
    
    #*********************#
    #'with expression corrected by study, use limma for model
    ### 8.6.5 JACKKNIFE:  run mixed models: random intercept (fast mode): correct for cell types
    #*********************#
    # see more examples in ('mixed_models/mixed_model.R')
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    # source('mixed_models/limma_model.R')
    
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    # include_study = F
    # model <- 'random_intercept'
    # disease_ls <- c('HD','AD')
    # phase_ls <- c('early','late')
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        if(disease == 'AD'){
            cell_types = AD_cell_types
        }else if(disease =='HD'){
            cell_types = HD_cell_types
        }
        for (phase in phase_ls){ ## loop2 by phase
            (data_dir_all <- paste0(disease_dir,'/mixed_model_jackknife/', model,model_keyword, '/', phase, '/')) ## where the mixed model input expression
            
            ## get the cell marker files
            cell_markers_f=max(list.files(paste0(home_dir, '/ND_results/cell_population/all_sample_estimation/', disease, '/', phase), recursive = T, full.names = T,
                                          pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' ))
            print(cell_markers_f)
            ## get all the run files:
            run_ls <- grep('run', list.dirs(data_dir_all, recursive = F), value = T)
            
            for(i in 1:length(run_ls)){#loo for each run
                (data_dir <- run_ls[i])
                (out_dir <- run_ls[i])
                print(paste0('START: ', data_dir))
                (exprdata <- paste0(data_dir, '/expression.Rdata'))
                ### get the expression data and MM results
                x <- limmaModel(exprdata,
                                phase,
                                out_dir,
                                cell_markers_f,
                                cell_types = cell_types,
                                include_study = include_study)
                assign(x = paste0(model, "_",phase,"_",disease), value = x)
            } ## loop3
        }## loop2
    }##loop1
    
    
    
    cat("
        #***************************************************************************#
        # with expression corrected by study, use limma for model
        # PART 8.6.6 JACKKNIFE: get the up and down list of genes for all the runs: correct for cell types
        #***************************************************************************#\n")
    ## get the up, down regulation mixed model results
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    source('mixed_models/compare_mm_meta.R')
    
    
    # disease_ls <- c('AD','HD')
    # model_ls <- c('random_intercept')
    # phase_ls <- c('early','late')
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        for (model in model_ls){## loop2 for models
            (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                      paste0(disease_dir, 'meta_analysis/meta_jack/')))
            (mm_dir = paste0(disease_dir, 'mixed_model_jackknife/',model,model_keyword,'/')) ## mixed model dir(parent dir)
            for(phase in phase_ls){
                mm_dir_p <- paste0(mm_dir, '/', phase,'/')
                (mm_dir_p <- paste0(grep('done|run', list.dirs(mm_dir_p, recursive = F), value = T), '/'))
                for(f in mm_dir_p){  # for each jackfile
                    print(f)
                    mainCompareMM(jack_meta_folder_ls, f, phase_ls=phase, compare_fisher =F)
                }
            }
        }## loop2 end
    }##loop1
    
    
    
    cat("
        #***************************************************************************#
        # with expression corrected by study, use limma for model
        # PART 8.6.7 JACKKNIFE: summary of jackknife ranks and prep for ermineJ results: correct for cell types
        #***************************************************************************#\n")
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    source('mixed_models/mixed_model_jackknife_results.R')
    
    # disease_ls <- c('AD','HD')
    # model_ls <- c('random_intercept')
    # phase_ls <- c('early','late')
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    # regulation_ls= c('up', 'down')
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        for (model in model_ls){## loop2 for models
            (mm_jack_dir = paste0(disease_dir, 'mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
            (mm_dir = paste0(disease_dir, 'mixed_model/random_intercept',model_keyword,'/')) ## MM parent dir
            for(phase in phase_ls){
                compareJackMM(mm_jack_dir,mm_dir, phase, regulation_ls)
            }
        }## loop2 end
    }##loop1
    
    
    cat("
        #***************************************************************************#
        # with expression corrected by study, use limma for model
        # PART 8.6.8a JACKKNIFE: use the genes in jaccknife as background: correct for cell types
        #***************************************************************************#\n")
    #' background is the same as input jaccknife genes (low expr rm, and inclusing NAs)
    # 
    #rm(list=setdiff(ls(),'home_dir'))
    source('mixed_models/mixed_model_jackknife_results.R')
    
    # disease_ls <- c('AD','HD')
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    
    ## gene, geneID, GO terms of all genes
    df_gene_anno<- read.delim(paste0(home_dir, 'ND_results/gene_annotation/all_gene_GO.tsv'), comment.char = '#')
    
    for (disease in disease_ls){## loop 1 for disease
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        (mm_jack_dir = paste0(disease_dir, 'mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
        bg_folder_out = paste0(mm_jack_dir, '/ermineJ_background/')
        dir.create(bg_folder_out, showWarnings = F, recursive = T)
        rownames(df_gene_anno) <- df_gene_anno$geneSymbol
        ## get the mm jack results (which contains all the genes for background)
        (f_ls <- list.files(mm_jack_dir, pattern = 'mixed_model_results.tsv'))
        for(i in 1: length(f_ls)){
            (f=paste0(mm_jack_dir, f_ls[i]))
            df <- read.delim(f, comment.char = '#')
            genes <- as.character(df$geneSymbol)
            bg <- df_gene_anno[genes, ] %>% droplevels()
            ## save bg file
            (f_out <- paste0(bg_folder_out, 'bg_', f_ls[i]))
            writeTable(bg, f_out)
            print(f_out)
        }
    }##loop1
    
    
    cat("
        #***************************************************************************#
        # with expression corrected by study, use limma for model
        # PART 8.6.8b JACKKNIFE: make ermineJ sh: correct for cell types
        #***************************************************************************#\n")
    ## make the sh script to run from meta and jack files with lowly expressed genes removed as background
    #rm(list=setdiff(ls(),'home_dir'))
    # 
    source('ermineJ_preprocess/make_ermineJ_sh.R')
    
    # disease_ls <- c('AD','HD')
    # model_ls <- c('random_intercept')
    #model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected' ## the expression.Rdata is already study corrected
    # xml = '/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz'
    for (disease in disease_ls){## loop 1 for disease
        print(disease)
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        for (model in model_ls){## loop2 for models
            (input_folder <- paste0(disease_dir, 'mixed_model_jackknife/', model, model_keyword, '/'))
            bg_folder <- paste0(input_folder, '/ermineJ_background/')
            erminej_dir <- paste0(disease_dir,'/ermineJ/mixed_model_jackknife/',model,model_keyword, '/', Sys.Date(),'/')
            
            mkErminejSH(disease,input_folder, bg_folder, erminej_dir,xml=xml)
        }## loop2 end
    }##loop1
    
}



checkErmineJResults <- function(disease_ls = c('AD','HD'),
                                model_ls = c('random_intercept'),
                                model_keyword = '_include_NA_low_exp_rm_adj_cell_pop_study_corrected',
                                keyword_ls = c('mixed_model'),
                                regulation_ls = c('up', 'down'),
                                phase_ls = c('early','late'),
                                top_threshold = 200,  ## the top genes to annotate
                                mixed_model_only=F,
                                threshold = 50){
    #threshold = 50,   ## number of top pathways to look at
    cat("
        #***************************************************************************#
        # with expression corrected by study, use limma for model
        # PART 8.6.8c JACKKNIFE: after MM ermineJ results, look at the top genes in the top pathways: correct for cell types
        #***************************************************************************#\n")
    for (disease in disease_ls){## loop 1 for disease
        print(disease)
        disease_dir= paste0(home_dir, disease, "_mouse_model_project/")
        
        for (model in model_ls){## loop2 for models
            
            ## where the mixed model result files
            (jack_meta_folder <- paste0(disease_dir, 'mixed_model_jackknife/', model,model_keyword,'/'))
            ## make sure only greps the folder with date (otherwise will grep the jackknife folder)
            (ermineJ_folder <- max(grep('2016', 
                                        list.dirs(paste0(disease_dir, 'ermineJ/mixed_model_jackknife/', model,model_keyword,'/'), recursive = F),
                                        value =T)))
            source('result_explore/top_genes_for_ermineJ.R')
            
        }## loop2 end
    }##loop1
}




