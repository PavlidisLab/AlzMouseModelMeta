## setwd to bin/

setwd(here::here('bin'))

#---------------------------------------------------------------------------#
# PART 0.2.1: Biological GO annotation for all genes from all platforms used in 3 diseases----
# for biological pathway only GO terms
# saved in '/ND_results/gene_annotation/all_gene_GO.tsv' used for ermineJ
#---------------------------------------------------------------------------#

print('PART 0.2.1: Biological GO annotation for all genes from all platforms used in 3 disease')

rm(list=setdiff(ls(),'home_dir'))
source("./ermineJ_preprocess/get_all_gene_annotation.R")

#### get all platforms from all diseases
platform_ls <- NULL

source('config_wrappers.R')
md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')
df_info <- read.delim(md_info, comment.char = '#')
platform <- unique(as.character(df_info$Platform))
platform_ls <- sort(unique(c(platform_ls, platform)))


## define platform_folder in the config file
platform_folder <- platform_folder_biological_pathway_only

(f_out<-paste0(home_dir, '/ND_results/gene_annotation/'))
dir.create(f_out, showWarnings = F, recursive = T)
f_out <- paste0(f_out, 'all_gene_GO.tsv')

df <- geneAnno(platform_ls, platform_folder, f_out, ermineJ_format = T, df_return = T)




#---------------------------------------------------------------------------#
# PART 0.2.2: all GO annotation for all genes from all platforms ----
# including cellular component, molecular function, biological processes
# saved in '/ND_results/gene_annotation/all_gene_GO_all_processes.tsv'' used for ermineJ
#---------------------------------------------------------------------------#
print('PART 0.2.2: all GO annotation for all genes from all platforms')

rm(list=setdiff(ls(),'home_dir'))
source("./ermineJ_preprocess/get_all_gene_annotation.R")

#### get all platforms 
platform_ls <- NULL

source('config_wrappers.R')
md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')
df_info <- read.delim(md_info, comment.char = '#')
platform <- unique(as.character(df_info$Platform))
platform_ls <- sort(unique(c(platform_ls, platform)))

platform_folder <- platform_folder_biological_all_GO ## all disease linked to the same folder

f_out<-paste0(home_dir, '/ND_results/gene_annotation/')
dir.create(f_out, showWarnings = F, recursive = T)
f_out <- paste0(f_out, 'all_gene_GO_all_processes.tsv')

df <- geneAnno(platform_ls, platform_folder, f_out, ermineJ_format = T, df_return = T)


#---------------------------------------------------------------------------#
# PART 1B: PROCESS gemma and CEL(normalized files) for limma DEA ----
# all files in /all_data/
#---------------------------------------------------------------------------#
print('PART 1B: PROCESS gemma and CEL(normalized files) for limma DEA ')

#***************************************************************************#
#**1B.a  for preparing datasets for limma -- from gemma & CEL -- results in all data ----
#***************************************************************************#

rm(list=setdiff(ls(),'home_dir'))
 
source('config_wrappers.R')

rm_outlier <- T
to_subset <- T
model_test <- T
explore_method <- 'analysis_mode'
output_folder <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/')
platform_folder <- platform_folder_biological_all_GO

dir.create(output_folder, recursive=T, showWarnings=F)

start_row <- 1
# dataset_todo <-'GSE64398.1'
#to_do_ls <-1

## for gemma unfiltered profiles GSE1556

f_dataset_ls <- paste0('../configs/', disease,'_mouse_dataset_doc/prioritized_gemma_explore_data.tsv')
unfiltered_exp_folder <- gemma_data  ## gemma processed
df_tmp <-read.delim(f_dataset_ls, comment.char = '#')
platform_folder <- platform_folder_biological_all_GO
if(nrow(df_tmp) > 0){
    source('explore_data_general.R')
    exploreExpData()
    }


### for CEL input and illumina
(f_dataset_ls <- paste0('../configs/', disease,'_mouse_dataset_doc/prioritized_CEL_mouse_datasets_explore.tsv'))
unfiltered_exp_folder <- paste0(disease_dir, 'data_and_QC/GEO_data/normalized_matrix/')
source('explore_data_general.R')
exploreExpData()



###

#---------------------------------------------------------------------------#
# PART 2.1: BATCH CORRECTION if needed ----
#---------------------------------------------------------------------------#
### batch correction
rm(list=setdiff(ls(),'home_dir'))
setwd(here::here('bin'))


source('config_wrappers.R')
source('batchCorrection.R')

(dataset_ls <- batch_correction_config)  # datasets for batch correction

data_folder <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
for(dataset in dataset_ls){
    object_dir = paste0(data_folder, dataset, '/results/')
    batchCorrection(dataset, object_dir = object_dir)
}


##

#---------------------------------------------------------------------------#
# PART 2.2: GENDER GENE expression (to look at lowly expressed thresholds) ----
#---------------------------------------------------------------------------#

rm(list=setdiff(ls(),'home_dir'))
 

source('config_wrappers.R')
source("./low_exp_rm_analysis/rm_low_exp_by_gender_genes.R")

## define folders
datadir <-paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')  ## analysis mode
df_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  #dataset label, phase label, order of the datasets etc.
plot_dir <- paste0(disease_dir, 'results/gender_expression_threshold/')
method_t_ls <-c("median") 
all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do


## grep the batch_correction.Rdata, or .Rdata files (if both exist, use batch corrected matrix)
r_obj_ls=NULL
for (folder in grep('results', list.dirs(datadir), value = T)){
    (r_obj <- grep('batch_correction.Rdata', list.files(folder, recursive=T, full.names=T), value=T))
    if(length(r_obj) <1){
        (r_obj <- grep('.Rdata', list.files(folder, recursive=T, full.names=T), value=T))
    }
    r_obj_ls <- c(r_obj_ls, r_obj)
}

# if(!all_dataset){
#     dataset_todo <-c('GSE63617.1', 'GSE63617.2')
#     dataset_todo <-c('GSE50521')
#     r_obj_ls <- grep(paste0(dataset_todo, collapse = '|'), r_obj_ls, value = T)
# }


df <- loopExpThreshold(method_t_ls, r_obj_ls, plot_dir, width=1000, height=800, df_return =T)


##

    #---------------------------------------------------------------------------#
    # PART 3.1 LIMMA DEA
    #---------------------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    # PART 3.1A LIMMA DEA for gemma and CEL and illumina
    #---------------------------------------------------------------------------#
## output put all limma (gemma and CEL) in 1 folder

### limma for Fisher only -CEL and gemma for 
rm(list=setdiff(ls(),'home_dir'))


source('config_wrappers.R')
q_threshold <- 0.1
datadir <-paste0(limma_dir, 'prioritized_limma/')
dir.create(datadir, showWarnings=F)
object_dir <- paste0(disease_dir, 'data_and_QC/all_data/explore_data/analysis_mode/')
f_model <- paste0('../configs/', disease,'_mouse_dataset_doc/limma_prioritized_for_meta_analysis.tsv')
combined_genotype =F   ## use original genotypes

plot_heatmap=T
write_table=T


start_row <- 1

#dataset_todo <-c('GSE48622')
#to_do_ls <-1
source('limma_DEA.R')


## 

    #---------------------------------------------------------------------------#
    # PART 6: Mixed Models include NAs
    #---------------------------------------------------------------------------#


    #***************************************************************************#
    # PART 6.1 get MM results (random intercept)
    #***************************************************************************#
## see 'mixed_models/mixed_model.R' for more examples
## no rm genes, run all genes
#*********************#
### 6.1.1a prepare combined expression data; , lowly expressed rm
#*********************#


rm(list=setdiff(ls(),'home_dir'))
setwd(here::here('bin'))
source('mixed_models/mixed_model.R')
NA_filter='0.3'  ## include genes with <0.3 are NAs
width = 1200
height = 1200

low_exp_rm=T
affy_threshold = 6
filter_method ='median'

source('config_wrappers.R')
model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')

model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder

for (disease in disease_ls){## loop 1 for disease
    # source('config_wrappers.R')
    
    ## #dataset label, phase label, order of the datasets etc.
    md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
    (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
    
    for (model in model_ls){## loop2 for models
        out_dir <- paste0(disease_dir, 'mixed_model/', model, model_keyword,'/')  ## name outdir by model
        for (phase in phase_ls){ ## loop 3 for phases
            x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
                                   width = width, height = height,
                                   affy_threshold = affy_threshold, 
                                   filter_method =filter_method,
                                   low_exp_rm=low_exp_rm)
        }
    }
}


###


#*********************#
### 6.1.2 run mixed models: random intercept, NA removed (all genes), lowly expressed rm
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print('6.1.2 run mixed models: random intercept, NA removed (all genes), lowly expressed rm')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/mixed_model.R')
full_report =T # to run fast estimateCI=F full_report = F for jackknife
estimateCI=T  # to run fast estimateCI=F full_report = F for jackknife


model_keyword = '_include_NA_low_exp_rm'  ## to specify which model folder
disease_stage_only =T  ## only use disease stage as fixed effect
tmp_result_output=F
rm_genes = ''
NA_filter='0.3'  ## NA_filter='0.3' to include NA
REML=F


model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for(model in model_ls){ ## loop2 by model
        out_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword , '/') ## where the mixed model result
        data_dir <- paste0(disease_dir,'/mixed_model/', model,model_keyword, '/') ## where the mixed model input expression
        
        for (phase in phase_ls){ ## loop3 by phase
            (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
            ### get the expression data and MM results
            print(Sys.time())
            x <- mixedModelAll(phase =phase, out_dir =out_dir, 
                               exprdata = exprdata, model=model,
                               full_report = full_report,
                               to_plot= F, rm_genes = rm_genes, NA_filter = NA_filter,
                               REML =REML,
                               estimateCI=estimateCI,
                               disease_stage_only =disease_stage_only,
                               tmp_result_output=tmp_result_output)
            assign(x = paste0(model, "_",phase,"_",disease), value = x)
        }
    }## loop2 end
}##loop1

#####

# see below to get an array with study (i.e. intercept) corrected

    #*********************#
    # PART 6.1.3 save a corrected array expression for gender and intercept (not plotting yet), see 11.3.2 for heatmap plotting
    #*********************#
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('helper_functions.R')
source('mixed_models/corrected_gender_intercept.R')


model_ls <- c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm')  ## this correct all the genes if intercept is provided
phase_ls <- c('early', 'late')
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')

correct_gender=F # dont need to correct gender when use disease only as FE

f_original_genotype <-'../configs/all_sample_design.tsv'
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            for(phase in phase_ls){
                (mm_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/', phase,'/'))
                correctInterceptGender(mm_dir, f_original_genotype,correct_gender=correct_gender)
            }
        }## loop 3 end
    }## loop2 end
}##loop1


####

    #***************************************************************************#
    # PART 6.1.4 prepare the ermineJ files for up and down gene lists
    # for random intercept
    #***************************************************************************#


rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/compare_mm_meta.R')


model_ls <- c('random_intercept')
phase_ls <- c('early','late')
model_keyword <- '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        mm_dir = paste0(disease_dir, '/mixed_model/',model,model_keyword,'/') ## mixed model dir(parent dir)
        mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls=phase_ls, compare_fisher =F)
    }## loop2 end
}##loop1


####

    #---------------------------------------------------------------------------#
    # PART 6C: JACKKNIFE: Mixed Models include NAs, low exp rm
    #---------------------------------------------------------------------------#


    #***************************************************************************#
    # PART 6C.1 JACKKNIFE: create expression Rdata for all studies
    #***************************************************************************#
## see 'mixed_models/mixed_model.R' for more examples
## no rm genes, run all genes


#*********************#
### 6C.1.1 JACKKNIFE: prepare combined expression data
#*********************#
print('6C.1.1 JACKKNIFE: prepare combined expression data')
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/mixed_model.R')
source('config_wrappers.R')

NA_filter='0.3'  ## include genes with <0.3 are NAs
width = 1200
height = 1200


model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')


model_keyword = '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder


low_exp_rm=T
affy_threshold = 6
filter_method ='median'
for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    ## #dataset label, phase label, order of the datasets etc.
    md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
    (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
    
    for (model in model_ls){## loop2 for models
        out_dir <- paste0(disease_dir, 'mixed_model_jackknife/', model, model_keyword,'/')  ## name outdir by model
        for (phase in phase_ls){ ## loop 3 for phases
            x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
                                   width = width, height = height,
                                   affy_threshold = affy_threshold, 
                                   filter_method =filter_method,
                                   low_exp_rm=low_exp_rm)
        }
    }
}



#####

#*********************#
### 6C.1.2 JACKKNIFE: creat expression.Rdata when 1 of the study is removed, in subfolder, run1...runX
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print('6C.1.2 JACKKNIFE: creat expression.Rdata when 1 of the study is removed, in subfolder, run1...runX')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')


model_keyword = '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder
library(dplyr)
for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for(model in model_ls){ ## loop2 by model
        
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
    }## loop2 end
}##loop1



#*********************#
### 6C.1.3 JACKKNIFE:  run mixed models: random intercept (fast mode)
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print('6C.1.3 JACKKNIFE:  run mixed models: random intercept (fast mode)')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/mixed_model.R')

model_keyword = '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder

full_report =T
rm_genes = ''
NA_filter='0.3'  ## NA_filter='0.3' to include NA
REML=F
estimateCI=F ## speed up, dont need to estimate 95%CI
disease_stage_only =T
tmp_result_output=F

model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for(model in model_ls){ ## loop2 by model
        for (phase in phase_ls){ ## loop3 by phase
            
            (data_dir_all <- paste0(disease_dir,'/mixed_model_jackknife/', model,model_keyword, '/', phase, '/')) ## where the mixed model input expression
            ## get all the run files:
            run_ls <- grep('run', list.dirs(data_dir_all, recursive = F), value = T)
            
            for(i in 1:length(run_ls)){#loop4 for each run
                (data_dir <- run_ls[i])
                print(paste0('START: ', data_dir))
                print(Sys.time())
                (exprdata <- paste0(data_dir, '/expression.Rdata'))
                ### get the expression data and MM results
                x <- mixedModelAll(phase =phase, out_dir =data_dir, 
                                   exprdata = exprdata, model=model,
                                   full_report = full_report,
                                   to_plot= F, rm_genes = rm_genes, NA_filter = NA_filter,
                                   REML =REML,
                                   estimateCI=estimateCI,
                                   disease_stage_only =disease_stage_only,
                                   tmp_result_output=tmp_result_output)
                assign(x = paste0(model, "_",phase,"_",disease,'_', i), value = x)
            } ## loop4
        }## loop3
    }## loop2 end
}##loop1



####
    #***************************************************************************#
    # PART 6C.2 JACKKNIFE: get the up and down list of genes per jackknife run
    #***************************************************************************#
## get the up, down regulation mixed model results
print('6C.2 JACKKNIFE: get the up and down list of genes per jackknife run')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/compare_mm_meta.R')


model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')
model_keyword <- '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        
        
        (mm_dir = paste0(disease_dir, '/mixed_model_jackknife/',model,model_keyword,'/')) ## mixed model dir(parent dir)
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




    #***************************************************************************#
    # PART 6C.3 JACKKNIFE: summary of jackknife ranks and prep for ermineJ results
    #***************************************************************************#
print('PART 6C.3 JACKKNIFE: summary of jackknife ranks and prep for ermineJ results')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/mixed_model_jackknife_results.R')
load('../configs/genenames.Rdata')

model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')
model_keyword <- '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder
regulation_ls= c('up', 'down')


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (mm_jack_dir = paste0(disease_dir, '/mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
        (mm_dir = paste0(disease_dir, '/mixed_model/random_intercept',model_keyword,'/')) ## MM parent dir
        for(phase in phase_ls){
            df_out_dir <- paste0(home_dir, '/results/ND_results/jackknife_rank_tables/')
            df_all <- compareJackMM(mm_jack_dir,mm_dir, phase, regulation_ls, return_df = T, df_out_dir=df_out_dir, 
                                    notes= 'Linear mixed model, with disease as fixed effect, and study as random effect',
                                    file_pre = '_before_MGP')
        }
    }## loop2 end
}##loop1






    #***************************************************************************#
    # PART 6C.4.1 JACKKNIFE: use the genes in jaccknife as background
    #***************************************************************************#
#' background is the same as input jaccknife genes (low expr rm, and inclusing NAs)
print('PART 6C.4.1 JACKKNIFE: use the genes in jaccknife as background')
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/mixed_model_jackknife_results.R')
source('config_wrappers.R')


model_keyword <- '_include_NA_low_exp_rm'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder
GO_annotation_dir <- '../configs/GO_annotation/all_gene_GO' # file prefix for GO annotation files
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    (mm_jack_dir = paste0(disease_dir, '/mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
    
    
    for(process_keyword in process_ls){
        bg_folder_out = paste0(mm_jack_dir, '/ermineJ_background',process_keyword, '/')
        df_gene_anno<- read.delim(paste0(GO_annotation_dir,process_keyword,'.tsv'), comment.char = '#')
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
    }## end loop for bio process only or all process
    
}##loop1




    #**************************#
    # PART 6D save the results of top genes (pvalue is raw p, not up or down p)
    # mixed model : no MGP correction
    #**************************#
#' rdata in mm_results_cell_adj
#' summerized top genes, top jack genes after cell type corrections
print('PART 6D save the results of top genes (pvalue is raw p, not up or down p)')
rm(list=setdiff(ls(),'home_dir'))

source('summary_tables/summary_MM_top_genes.R')
source('config_wrappers.R')
mgi_info <- paste0('../configs/ND_files/mgi_annotation_precessed_omim.tsv')

model_ls <- c('mixed_model_jackknife')  ## specify which mixed model results

f_out <- paste0(home_dir, '/ND_results/tables/top_genes/top_genes_non_adj/')
f_rdata_out <- paste0(home_dir, '/ND_results/mm_results/', Sys.Date(), '/')  ## save the r data
threshold =50
phase_ls <- c('early', 'late')

mm_folder ='random_intercept_include_NA_low_exp_rm'

## must run this for the mixed model first (because the estimate and std are from full mixed model)
for(phase in phase_ls){
    summarizeTopGenes(disease_ls, phase =phase, f_out,f_rdata_out = f_rdata_out, threshold =threshold, top_genes = T, mm_folder= mm_folder,
                      mgi_info =mgi_info )
}

## get the top genes for jackknife, and save rdata
(f_out <- paste0(home_dir, '/ND_results/tables/top_genes/jack_top_genes_non_adj/'))
mm_jack_folder ='random_intercept_include_NA_low_exp_rm'
threshold =50
for(phase in phase_ls){
    summarizeTopJackGenes(disease_ls, phase =phase, f_out,f_rdata_out = f_rdata_out, threshold =threshold, top_genes = T, 
                          mm_jack_folder=mm_jack_folder)
}




    #---------------------------------------------------------------------------#
    ### PART 7 (previously 3.3.2c)  MGP estimation, Cell population estimation from mixed model----
    #' array must be study (intercept)corrected data
    #' cell type profiles are used for the mixed model
    #---------------------------------------------------------------------------#
print('PART 7 (previously 3.3.2c)  MGP estimation, Cell population estimation from mixed model')
## output: '/MGP_estimation/
## recorded done the cell type markers for future heatmaps
library(markerGeneProfile)
library(ogbox)
library(HelperFunctions)

rm(list=setdiff(ls(),'home_dir'))

source('config_wrappers.R')

# load('../configs/mouseMarkerGenes.Rdata') ## load the markergenes
# load('../configs/mouseMarkerGenes.rda') ## load the markergenes
load('../configs/mouseMarkerGenesCombined.rda') ## load the markergenes
# ogbox::loadGithub('oganm/neuroExpressoAnalysis/data/mouseMarkerGenesCombined.rda')
mouseMarkerGenes = mouseMarkerGenesCombined

phase_ls = c('early', 'late')


all_dataset <- T  # if to do all datasets, if not define dataset_todo to specify which dataset to do
design_group = 'Genotype'  # design_group = 'Original_genotype'
sigTest = wilcox.test 
wt_only = T  ## only compare disease to WT

threshold = -10 # filter probes exp less than threshold (input is study corrected value, and low expression filtered)
mixed_model =T
removeNegatives = F  ## ogan suggested use F, the new version don't use this argument (sets the param (markerGeneProfile::fullEstimate(removeMinority = )

for (disease in disease_ls){
    for(phase in phase_ls){
        source('config_wrappers.R')
        ## define marker genes
        if(disease =='AD'){
            genes = mouseMarkerGenes$Hippocampus
            print('Input Hippocampus markers')
        }else{
            genes = mouseMarkerGenes$Striatum
            print('Input striatum markers')
            
        }
        
        (file_ls <- paste0(disease_dir, '/MGP_estimation/', phase,'/'))
        
        
        (r_ob <- list.files(paste0(disease_dir, '/mixed_model/random_intercept_include_NA_low_exp_rm/', phase,'/'), 
                            full.names=T, recursive=T, pattern= "mixed_model_results_exp_corrected.Rdata"))
        (outDir <-paste0(file_ls,Sys.Date(), '/'))
        dir.create(outDir,showWarnings=F, recursive=T)
        
        source("MGP_estimation/estimate_cell_population.R")
    }
}



    #---------------------------------------------------------------------------#
    # PART 8: mixed model, and JACCKNIFE mixed model: correct for cell types
    #---------------------------------------------------------------------------#
#' for late HD and late AD


    #***************************************************************************#
    # PART 8.1 get MM results (random intercept) for all samples : correct for cell types
    # saved in 'mixed_model/random_intercept_include_NA_adj_cell_pop'
    #***************************************************************************#
## see 'mixed_models/mixed_model.R' for more examples
## no rm genes, run all genes
#*********************#
### 8.1.1 prepare combined expression data: correct for cell types, filter lowly expressed
#*********************#

print('8.1.1 prepare combined expression data: correct for cell types, filter lowly expressed')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

source('mixed_models/mixed_model.R')
NA_filter='0.3'  ## include genes with <0.3 are NAs
width = 1200
height = 1200

model_ls <- c('random_intercept')
phase_ls <- c('early','late')

model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'

low_exp_rm=T
affy_threshold = 6
filter_method ='median'

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    
    ## #dataset label, phase label, order of the datasets etc.
    md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
    (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
    
    for (model in model_ls){## loop2 for models
        out_dir <- paste0(disease_dir, 'mixed_model/', model, model_keyword,'/')  ## name outdir by model
        for (phase in phase_ls){ ## loop 3 for phases
            x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
                                   width = width, height = height,
                                   affy_threshold = affy_threshold, 
                                   filter_method =filter_method,
                                   low_exp_rm=low_exp_rm)
        }
    }
}


#*********************#
# correct for cell types
### 8.1.2 run mixed models: random intercept, NA removed (all genes): correct for cell types
#' updated 2016-10-30: use the new population estimation from study corrected expressions of all samples in the phase
#' input MGP estimations 
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print('8.1.2 run mixed models: random intercept, NA removed (all genes): correct for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

source('mixed_models/mixed_model.R')
source('mixed_models/mixed_model_adj_celltype.R')
full_report =T # to run fast estimateCI=F full_report = F
estimateCI=F  # to run fast estimateCI=F full_report = F


model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'

disease_stage_only =T  ## only use disease stage as fixed effect
tmp_result_output=F
rm_genes = ''
NA_filter='0.3'  ## NA_filter='0.3' to include NA
REML=F



model_ls <- c('random_intercept')
phase_ls <- c('early', 'late')

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    
    for(model in model_ls){ ## loop2 by model
        out_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword , '/') ## where the mixed model result
        data_dir <- paste0(disease_dir,'/mixed_model/', model,model_keyword, '/') ## where the mixed model input expression
        
        for (phase in phase_ls){ ## loop3 by phase
            (exprdata <- paste0(data_dir, phase, '/expression.Rdata'))
            ## get the cell marker files
            cell_markers_f=max(list.files(paste0(disease_dir, '/MGP_estimation/', phase,'/'), 
                                          recursive = T, full.names = T,
                                          pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' ))
            
            print(cell_markers_f)
            ### get the expression data and MM results
            cat(paste0('\n', Sys.time(), '\n PROCESSING ', disease, ': ', phase, '\n'))
            x <- mixedModelAllCellType(phase =phase, out_dir =data_dir, 
                                       exprdata = exprdata, model=model,
                                       full_report = full_report,
                                       to_plot= F, rm_genes = rm_genes, NA_filter = NA_filter,
                                       REML =REML,
                                       estimateCI=estimateCI,
                                       disease_stage_only =disease_stage_only,
                                       tmp_result_output=tmp_result_output,
                                       cell_markers_f=cell_markers_f)
            assign(x = paste0(model, "_",phase,"_",disease), value = x)
        }
    }## loop2 end
}##loop1






    #***************************************************************************#
    # -- lowly expressed genes removed : correct for cell types
    # PART 8.1.4 prepare the ermineJ files for up and down gene lists
    # for random intercept
    #***************************************************************************#

print('PART 8.1.4 prepare the ermineJ files for up and down gene lists')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/compare_mm_meta.R')


model_ls <- c('random_intercept')
phase_ls <- c('early','late')
model_keyword <- '_include_NA_low_exp_rm_adj_cell_pop'  ## if with NA: model_keyword = '_include_NA'  ## to specify which model folder


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        mm_dir = paste0(disease_dir, 'mixed_model/',model,model_keyword,'/') ## mixed model dir(parent dir)
        mainCompareMM(jack_meta_folder_ls, mm_dir, phase_ls=phase_ls, compare_fisher =F)
    }## loop2 end
}##loop1




    #***************************************************************************#
    # PART 8.1.5 correct for cell types and intercept (not plotting yet)
    #***************************************************************************#

print('PART 8.1.5 correct for cell types and intercept (not plotting yet)')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

source('helper_functions.R')
source('mixed_models/corrected_cell_type_intercept.R')


model_ls <- c('random_intercept')
model_keyword_ls <- c('_include_NA_low_exp_rm_adj_cell_pop')  ## this correct all the genes if intercept is provided
phase_ls <- c('early', 'late')
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')

correct_gender=F # dont need to correct gender when use disease only as FE
#cell_markers_f ="/home/bzhuang/ND_project_combined//ND_results/cell_population/2016-10-26/cell_population_estimate_2016-10-26.tsv"

f_original_genotype <- '../configs/all_sample_design.tsv'
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    if(disease == 'AD'){
        cell_types = c("Astrocyte" ,"DentateGranule", 'GabaSSTReln',"Microglia", 'Oligo', 'Pyramidal')
    }else if(disease =='HD'){
        cell_types = c("Astrocyte","ForebrainCholin" ,"Microglia" ,"Spiny",'Oligo')
    }
    
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            for(phase in phase_ls){
                print(phase)
                cell_markers_f=max(list.files(paste0(disease_dir, '/MGP_estimation/', phase,'/'), 
                                              recursive = T, full.names = T,
                                              pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' ))
                print(cell_markers_f)
                
                (mm_dir <- paste0(disease_dir, 'mixed_model/', model,model_keyword,'/', phase,'/'))
                correctInterceptCellType(mm_dir, f_original_genotype,correct_gender=correct_gender,cell_markers_f=cell_markers_f, cell_types=cell_types)
            }
        }## loop 3 end
    }## loop2 end
}##loop1






    #***************************************************************************#
    # PART 8.2 JACKKNIFE: create expression Rdata for all studies : correct for cell types
    #***************************************************************************#
## see 'mixed_models/mixed_model.R' for more examples
## no rm genes, run all genes
#*********************#
### 8.2.1 JACKKNIFE: prepare combined expression data : correct for cell types
#*********************#
print(' 8.2.1 JACKKNIFE: prepare combined expression data : correct for cell types')

rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/mixed_model.R')
NA_filter='0.3'  ## include genes with <0.3 are NAs
width = 1200
height = 1200


model_ls <- c('random_intercept')
phase_ls <- c('early','late')


model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'


low_exp_rm=T
affy_threshold = 6
filter_method ='median'
for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    ## #dataset label, phase label, order of the datasets etc.
    md_info <- paste0('../configs/', disease,'_mouse_dataset_doc/dataset_info_mixed_model.tsv')  
    datadir <-paste0(limma_dir, 'prioritized_limma/')  ## dir for limma robjects
    (rfile_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
    
    for (model in model_ls){## loop2 for models
        out_dir <- paste0(disease_dir, 'mixed_model_jackknife/', model, model_keyword,'/')  ## name outdir by model
        for (phase in phase_ls){ ## loop 3 for phases
            x <- prepareAllSamples(rfile_ls, md_info, phase, out_dir, NA_filter=NA_filter,
                                   width = width, height = height,
                                   affy_threshold = affy_threshold, 
                                   filter_method =filter_method,
                                   low_exp_rm=low_exp_rm)
        }
    }
}





#*********************#
### 8.2.2 JACKKNIFE: creat expression.Rdata when 1 of the study is removed, in subfolder, run1...runX
# correct for cell types
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print(' 8.2.2 JACKKNIFE: creat expression.Rdata when 1 of the study is removed, in subfolder, run1...runX')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')

model_ls <- c('random_intercept')
phase_ls <- c('early','late')


model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'
library(dplyr)
for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for(model in model_ls){ ## loop2 by model
        
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
    }## loop2 end
}##loop1




#*********************#
### 8.2.3 JACKKNIFE:  run mixed models: random intercept (fast mode): correct for cell types
#*********************#
# see more examples in ('mixed_models/mixed_model.R')
print(' 8.2.3 JACKKNIFE:  run mixed models: random intercept (fast mode): correct for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('mixed_models/mixed_model.R')
source('mixed_models/mixed_model_adj_celltype.R')

model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'

full_report =T
rm_genes = ''
NA_filter='0.3'  ## NA_filter='0.3' to include NA
REML=F
estimateCI=F ## speed up, dont need to estimate 95%CI
disease_stage_only =T
tmp_result_output=F

model <- 'random_intercept'
phase_ls <- c('early','late')

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (phase in phase_ls){ ## loop2 by phase
        (data_dir_all <- paste0(disease_dir,'/mixed_model_jackknife/', model,model_keyword, '/', phase, '/')) ## where the mixed model input expression
        
        ## get the cell marker files
        cell_markers_f=max(list.files(paste0(disease_dir, '/MGP_estimation/', phase,'/'), 
                                      recursive = T, full.names = T,
                                      pattern ='mixed_model_cell_proportion_estimation_scaled.tsv' ))
        print(cell_markers_f)
        ## get all the run files:
        run_ls <- grep('run', list.dirs(data_dir_all, recursive = F), value = T)
        
        for(i in 1:length(run_ls)){#loo for each run
            (data_dir <- run_ls[i])
            cat(paste0('\n##############################\n', Sys.time(), '\n PROCESSING ', disease, ': ', phase,': run ', i, '\n##############################\n'))
            print(paste0('START: ', data_dir))
            (exprdata <- paste0(data_dir, '/expression.Rdata'))
            ### get the expression data and MM results
            x <- mixedModelAllCellType(phase =phase, out_dir =data_dir, 
                                       exprdata = exprdata, model=model,
                                       full_report = full_report,
                                       to_plot= F, rm_genes = rm_genes, NA_filter = NA_filter,
                                       REML =REML,
                                       estimateCI=estimateCI,
                                       disease_stage_only =disease_stage_only,
                                       tmp_result_output=tmp_result_output,
                                       cell_markers_f=cell_markers_f)
            assign(x = paste0(model, "_",phase,"_",disease,'_', i), value = x)
        } ## loop3
    }## loop2
}##loop1




    #***************************************************************************#
    # PART 8.3.1 JACKKNIFE: get the up and down list of genes for all the runs: correct for cell types
    #***************************************************************************#
## get the up, down regulation mixed model results
print(' PART 8.3.1 JACKKNIFE: get the up and down list of genes for all the runs: correct for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/compare_mm_meta.R')
source('config_wrappers.R')


model_ls <- c('random_intercept')
phase_ls <- c('early','late')
model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'

for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (jack_meta_folder_ls <- c(paste0(disease_dir, 'meta_analysis/low_exp_rm/'),
                                  paste0(disease_dir, 'meta_analysis/meta_jack/')))
        (mm_dir = paste0(disease_dir, '/mixed_model_jackknife/',model,model_keyword,'/')) ## mixed model dir(parent dir)
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




    #***************************************************************************#
    # PART 8.3.2 JACKKNIFE: summary of jackknife ranks and prep for ermineJ results: correct for cell types
    #***************************************************************************#
print(' PART 8.3.2 JACKKNIFE: summary of jackknife ranks and prep for ermineJ results: correct for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/mixed_model_jackknife_results.R')
source('config_wrappers.R')
load('../configs/genenames.Rdata')

model_ls <- c('random_intercept')
phase_ls <- c('early','late')
model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'
regulation_ls= c('up', 'down')


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    for (model in model_ls){## loop2 for models
        (mm_jack_dir = paste0(disease_dir, '/mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
        (mm_dir = paste0(disease_dir, '/mixed_model/random_intercept',model_keyword,'/')) ## MM parent dir
        for(phase in phase_ls){
            df_out_dir <- paste0(home_dir, '/ND_results/jackknife_rank_tables_adj_cell_pop/')
            df_all <- compareJackMM(mm_jack_dir,mm_dir, phase, regulation_ls, return_df = T, df_out_dir=df_out_dir,
                                    notes ='Linear mixed model with marker gene profiles correction, with disease as fixed effect, and study as random effect, and each cell type as fixed effect',
                                    file_pre = '_after_MGP')
        }
    }## loop2 end
}##loop1





    #***************************************************************************#
    # PART 8.4.1 JACKKNIFE: use the genes in jaccknife as background: corrected for cell types
    #***************************************************************************#
#' background is the same as input jaccknife genes (low expr rm, and inclusing NAs)
print(' PART 8.4.1 JACKKNIFE: use the genes in jaccknife as background: corrected for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('mixed_models/mixed_model_jackknife_results.R')

source('config_wrappers.R')



model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'
GO_annotation_dir <- '../configs/GO_annotation/all_gene_GO' # file prefix for GO annotation files
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories


## gene, geneID, GO terms of all genes
## for biological process only


for (disease in disease_ls){## loop 1 for disease
    source('config_wrappers.R')
    (mm_jack_dir = paste0(disease_dir, 'mixed_model_jackknife/random_intercept',model_keyword,'/')) ## jackknife MM parent dir
    
    for(process_keyword in process_ls){
        bg_folder_out = paste0(mm_jack_dir, '/ermineJ_background',process_keyword, '/')
        df_gene_anno<- read.delim(paste0(GO_annotation_dir,process_keyword,'.tsv'), comment.char = '#')
        
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
    }
}##loop1




    #***************************************************************************#
    # PART 8.4.2 JACKKNIFE: make ermineJ sh: correct for cell types ----
    #***************************************************************************#
#' This will generate all the enrichment results.
#' define variable `xml` in the `config_wrappers.R` 
#' `go_daily-termdb.rdf-xml.gz` is the gene ontology database, 
#' the latest `go_daily-termdb.rdf-xml.gz` can be downloaded from http://archive.geneontology.org/latest-termdb/


## make the sh script to run from meta and jack files with lowly expressed genes removed as background
rm(list=setdiff(ls(),'home_dir'))
print('PART 8.4.2 JACKKNIFE: make ermineJ sh: correct for cell types')
source('ermineJ_preprocess/make_ermineJ_sh.R')
source('config_wrappers.R')

model_ls <- c('random_intercept')
model_keyword = '_include_NA_low_exp_rm_adj_cell_pop'

maxsize = 500
iteration= '200000'
## get both biological and all process enrichment
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories
ermineJShs = c()
for (disease in disease_ls){## loop 1 for disease
    print(disease)
    for (model in model_ls){## loop2 for models
        (input_folder <- paste0(disease_dir, '/mixed_model_jackknife/', model, model_keyword, '/'))
        for(process in process_ls){
            bg_folder <- paste0(input_folder, '/ermineJ_background',process,'/')
            erminej_dir <- paste0(disease_dir,'/ermineJ/mixed_model_jackknife/',model,model_keyword,process, '/', Sys.Date(),'_geneset_', maxsize, '/')
            y = mkErminejSH(disease,input_folder, bg_folder, erminej_dir,xml=xml, maxsize=maxsize,iteration = iteration)
            ermineJShs = c(ermineJShs,y)
        }
        
    }## loop2 end
}##loop1



### run ermineJ sh files
# set java stuff
ermineR:::findJava()
# ermineJHome = system.file("ermineJ-3.1", package = "ermineR")
Sys.setenv(ERMINEJ_HOME = here::here('ermineJ-3.0.2'))

ermineJShs %>% lapply(function(x){
    system2(x)
})


    #***************************************************************************#
    # PART 8.4.3 JACKKNIFE: after MM ermineJ results ----
    # look at the top genes in the top pathways: correct for cell types
    # see part 11 for the cell pop adj
    #***************************************************************************#
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
print('PART 8.4.3 JACKKNIFE: after MM ermineJ results, look at the top genes in the top pathways: correct for cell types')


model_ls <- c('random_intercept')
model_keyword_ls = '_include_NA_low_exp_rm_adj_cell_pop'
keyword_ls <- c('mixed_model')
regulation_ls <- c('up', 'down')
phase_ls <- c('early','late')


top_threshold = 200  ## the top genes to annotate
mixed_model_only=F
threshold = 50   ## number of top pathways to look at
fdr_threshold = 0.05  # or sig_paths to look at if more than threshold
process_ls <- c('', '_all_processes')  # '' is the biolofical process only, and '_all_processes' are with all 3 pathway categories

for (disease in disease_ls){## loop 1 for disease
    print(disease)
    source('config_wrappers.R')
    
    for (model in model_ls){## loop2 for models
        for(model_keyword in model_keyword_ls){ ## loop3, with or no NA
            ## where the mixed model result files
            (jack_meta_folder <- paste0(disease_dir, 'mixed_model_jackknife/', model,model_keyword,'/'))
            
            for(process in process_ls){
                ## make sure only greps the folder with date (otherwise will grep the jackknife folder)
                (ermineJ_folder <- max(grep('201', 
                                            list.dirs(paste0(disease_dir, 'ermineJ/mixed_model_jackknife/', model,model_keyword,process, '/'), recursive = F),
                                            value =T)))
                source('result_explore/top_genes_for_ermineJ.R')
            }
            
        }## loop 3 end
    }## loop2 end
}##loop1




    #**************************#
    # PART 8.6.1 save the results of top genes (pvalue is raw p, not up or down p)
    # mixed model and jackknife results :corrected for cell types
    #**************************#
#' rdata in mm_results_cell_adj
#' summerized top genes, top jack genes after cell type corrections
print('PART 8.6.1 save the results of top genes (pvalue is raw p, not up or down p)')
rm(list=setdiff(ls(),'home_dir'))

source('summary_tables/summary_MM_top_genes.R')
source('config_wrappers.R')
mgi_info <- paste0('../configs/ND_files/mgi_annotation_precessed_omim.tsv')

model_ls <- c('mixed_model_jackknife')  ## specify which mixed model results

f_out <- paste0(home_dir, '/ND_results/tables/top_genes/top_genes_cell_adj/')
f_rdata_out <- paste0(home_dir, 'ND_results/mm_results_cell_adj/', Sys.Date(), '/')  ## save the r data
threshold =50
phase_ls <- c('early', 'late')

mm_folder ='random_intercept_include_NA_low_exp_rm_adj_cell_pop'

## must run this for the mixed model first (because the estimate and std are from full mixed model)
for(phase in phase_ls){
    summarizeTopGenes(disease_ls, phase =phase, f_out,f_rdata_out = f_rdata_out, threshold =threshold, top_genes = T, mm_folder= mm_folder,
                      mgi_info =mgi_info )
}

## get the top genes for jackknife, and save rdata
(f_out <- paste0(home_dir, '/ND_results/tables/top_genes/jack_top_genes_cell_adj/'))
mm_jack_folder ='random_intercept_include_NA_low_exp_rm_adj_cell_pop'

threshold =50
for(phase in phase_ls){
    summarizeTopJackGenes(disease_ls, phase =phase, f_out,f_rdata_out = f_rdata_out, threshold =threshold, top_genes = T, 
                          mm_jack_folder=mm_jack_folder)
}



    #**************************#
    # PART 8.6.2 get the top pathways:correct for cell types
    #**************************#
#' top pathways in doc_tables/pathways_cell_adj
print('PART 8.6.2 get the top pathways:correct for cell types')
rm(list=setdiff(ls(),'home_dir'))
source('config_wrappers.R')
source('summary_tables/summary_MM_top_genes.R')


model_ls <- c('mixed_model_jackknife')  ## specify which mixed model results

mm_jack_folder ='random_intercept_include_NA_low_exp_rm_adj_cell_pop'

f_out <- paste0(home_dir, '/ND_results/ND_summary/doc_tables/pathways_cell_adj/')
summarizeTopPathways(disease_ls, f_out,model_ls,
                     phase_ls=c('early', 'late'),
                     regulation_ls = c('up', 'down'),
                     threshold =50, mm_jack_folder=mm_jack_folder)
