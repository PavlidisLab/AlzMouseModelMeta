#' Author: Beryl Zhuang
#' created: 2016-01-10
#' updated: 2016-05-10
#' 2017-04-12: publish_codes
#' 
#' given a list of GEO IDs, download the raw CEL files, make a metadata file similar to gemma style
#' USAGE:
# source('GEO_download_CEL.R')
# GEO <- c("GSE15058", "GSE66598", "GSE74995") # 2016-03-16 other tissue with WT early phase
# CEL_dir <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/'
# downloadCEL(GEO, CEL_dir)
# ##get the sample info
# for(gse in GEO){
#     getSampleMeta(gse, CEL_dir)
# }




getSampleMeta <- function(gse, CEL_dir){
    ## get the metadata and make it similar to gemma style
    library('GEOquery')
    gse_df <- getGEO(gse)
    for(i in 1:length(gse_df)){  # loop to save every platforms
        df <- gse_df[[i]]
        # get the sample meta data
        samples <- as.data.frame(df@phenoData@data)
        colnames(samples)[which(colnames(samples)=="geo_accession")] = 'ExternalID'
        samples$Bioassay <- paste0(gse, '_Biomat_0___BioAssayImplId=99999Name=', gsub(' ' , '_', samples$title))
        
        file_name <- paste0(CEL_dir, '/', gse, '/sample_metadata_', i,'.tsv')
        dir.create(paste0(CEL_dir, '/', gse), recursive = T,showWarnings = F)
        write.table(samples, file = file_name, row.names = F,
                    sep ='\t', quote = F)
        print(paste0(file_name, ' saved'))
    }

}

downloadCEL <- function(geo, CEL_dir, meta_data = T){
  library('GEOquery')
  dir.create(CEL_dir, showWarnings=F, recursive = T)

  for (i in geo){
    print(paste0('Downloading ', i))
    res <- try(getGEOSuppFiles(i, makeDirectory = TRUE, baseDir = CEL_dir))
    if(inherits(res, 'try-error')){
        print(paste0('error, ', i, '; continue to the next GEO'))
    }
    if(meta_data){
        getSampleMeta(gse = i, CEL_dir= CEL_dir)
    }

  }
}

