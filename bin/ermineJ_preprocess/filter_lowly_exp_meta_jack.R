# 2016-05-26
# Updated 2016-05-26
#' @require ermineJ_preprocess/filter_out_low_exp_for_enrichment.R, filterLowExprGenes() function
#' use meta and jack files (has Go terms) and limma object (expression matrix) to remove genes that dont have any probes
#' that are expressed, and return a meta/ jack file with genes with at least 1 probes expressed
#' 
#' @example 
#  
# rm(list=setdiff(ls(),'home_dir'))
# disease = 'HD'
# source('config/config_wrappers.R')
# source("./ermineJ_preprocess/filter_lowly_exp_meta_jack.R")
# 
# ## meta, jack folders
# (jack_meta_folder <- paste0(disease_dir, '/meta_analysis/meta_jack/'))
# r_object_dir <-paste0(limma_dir, 'prioritized_limma/')# the limma object
# df_info <- paste0(disease_dir, 'config_files/dataset_info.tsv') # for early, intermediate and late files
# out_dir <- paste0(disease_dir, '/meta_analysis/low_exp_rm/')
# affy_threshold = 6
# agilent_threshold = -10
# filter_method = "median"
# 
# makeFilteredJackmeta(jack_meta_folder = jack_meta_folder, r_object_dir=r_object_dir, df_info=df_info, 
#                      out_dir=out_dir, 
#                     affy_threshold = affy_threshold, agilent_threshold = agilent_threshold, 
#                     filter_method = filter_method)

 
source('helper_functions.R')
source("./ermineJ_preprocess/filter_out_low_exp_for_enrichment.R")

makeFilteredJackmeta <- function(jack_meta_folder, r_object_dir, df_info, out_dir, 
         affy_threshold = 6, agilent_threshold = -10, filter_method = "median"){
    
    dir.create(out_dir, recursive = T, showWarnings = F)
    ## get the disease phases from dataset info
    dataset_df <- read.delim(df_info, comment.char = "#")
    (disease_phase_ls <- setdiff(union(levels(dataset_df$Phase), levels(dataset_df$Extra_phase)), ""))
    
    ## loop 1 for each disease stat
    for (i in 1: length(disease_phase_ls)){
        (disease_phase <- disease_phase_ls[i])
        print(paste0("DISEASE PHASE: ", disease_phase))
        ## grep the corresponding meta and jack file (including up and down, meta and jack)
        (f_m_j_ls <- grep(disease_phase, list.files(jack_meta_folder, pattern = 'jackknife.tsv|meta_genes.tsv', full.names = T), value = T) )
        dir.create(out_dir, recursive = T, showWarnings = F)
        
        (f_ls_all <- list.files(path = r_object_dir, recursive = T, pattern= '.Rdata', full.names=T))
        ## a df with 'Dataset', 'Timepoint', 'Phase', 'Order', 'Genotype', 'File'(limma objectdir), file_labels(for plt)
        info_df <- getRObjectLabels(f_ls_all, df_info)
        
        ## get the studies corresponding to the disease phase
        (selected_study <- filterContain(info_df, "Phase", disease_phase))
        if(nrow(selected_study) == 0){  # for more phases
            selected_study <- filterContain(info_df, "Extra_phase", disease_phase)
        }
        
        (f_ls <- selected_study$File) # contain a list of object dirs of the specified phase
        
        ## get the overlap of genes that are expressed in these matrix (no genotype is filtered)
        x <- filterLowExprGenes(file_dir_ls = f_ls, affy_threshold = affy_threshold, agilent_threshold = agilent_threshold, 
                                filter_method = filter_method, return_array = F)
        df_all <- x[[1]]
        msg_filter <- x[[3]]
        
        ## loop 2 for a disease stat, for each meta and jack files of that disease stat
        ## now filter the correponding jack/meta file (by phase), and save a new file
        for(f_m_j in f_m_j_ls){
            (prefix <- grep('.tsv', unlist(strsplit(f_m_j, split = "/")), value = T))
            (f_o <- paste0(out_dir, "/low_expr_rm_", prefix))
            df <- read.delim(f_m_j, comment.char = "#")
            
            (original_rows <-  nrow(df))
            ## get the gene name as row names
            row.names(df) <- df$geneSymbol
            
            
            ## remove the lowly exp genes
            exp_gene_list <- intersect(df_all, df$geneSymbol)
            df_filtered <- df[exp_gene_list, ]
            (low_rows <- original_rows - nrow(df_filtered))
            
            ## rm NA
            rm_col <- intersect(c('adj_combined_max_p','Fisher'), colnames(df_filtered))
            rm_row <- which(is.na(df_filtered[, rm_col]))
            if(length(rm_row)>0){df_f <- df_filtered[-rm_row, ]}
            (na_rows <- nrow(df_filtered) - nrow(df_f))
            
            ## order by col
            df_f <- df_f[order(df_f[, rm_col]), ]
            
            ## log
            msg <- paste0("# ", Sys.Date(),
                          "\n# Input file: ", f_m_j,
                          "\n# Genes before filter: ",  original_rows,
                          "\n# Lowly expressed gene count: ", low_rows,
                          "\n# NA filtered: ", na_rows,
                          "\n# Final gene count: ", nrow(df_f),"\n#\n",
                          msg_filter)
            
            ## write table
            writeTable(df_f,f_out = f_o, msg = msg) 
            print(paste0("OUTPUT FILE: ", f_o))
            
            ## write ermineJ background file with
            dir.create(paste0(out_dir, "/ermineJ_background/"), recursive = T, showWarnings = F)
            (f_e_o <- paste0(out_dir, "/ermineJ_background/bg_low_expr_rm_", prefix))
            
            ## ermineJ background must have geneID (same as geneSymbol), geneSymbol, geneNames, GOterms
            df_erminej_bg <- cbind(data.frame(geneID = df_f$geneSymbol),
                                   df_f[, c('geneSymbol', 'geneNames', 'GOTerms')])
            
            
            writeTable(df_erminej_bg,f_out = f_e_o, msg = "") 
            print(paste0("BACKGOURND OUTPUT FILE: ", f_e_o))
        }
    }
    print("DONE!")
}





