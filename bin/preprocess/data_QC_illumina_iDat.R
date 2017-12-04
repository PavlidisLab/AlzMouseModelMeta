#' 2016-06-07
#' Process illumnia idat and non normalized data by limma package
#' The neqc functions performs normexp background correction using negative controls, then quantile
#' normalizes and finally log2 transforms [33]. It also automatically removes the control probes, leaving
#' only the regular probes 
#' 
#' all tar, gz files must be unzipped

 
source("helper_functions.R") ## for saveSampleInfo()
library(RColorBrewer)
library(HelperFunctions) # helpers
library(limma)




iDatDataQC <- function(file_dir, out_dir, data_dir, 
                           design_f,bgx='',
                           width=1600, height=1000,
                           sep="\t"){
    # bgx use absolute path
    #'
    #' INPUT:
    #'  - file_dir: dir with all iDat files and "SampleInfo.txt"; or non-normalized data
    #'  - file_dir, if only 1 file for all expression data, the design_f is required
    #'  - design_f: file for all experimental designs
    #'  - out_dir: out_dir for the QC results and plots
    #'  - data_dir: output of matrix
    #'  - bgx is for idat data, non-normalized data don't need this
    #'  - sep: sept for the non normalized files
    #'  
    #'  OUTOUT:
    #'  - boxplot before and after, and save NEQC matrix
    current_wd <- getwd()
    #*****************************#
    # set parameters
    #*****************************#
    dataset <- getGSEID(file_dir)
    file_dir <- paste0(file_dir,'/')
    (plot_pre <- paste0(out_dir,  "/", dataset, "_"))
    if(data_dir ==""){
        data_pre <-paste0(out_dir,  "/", dataset, "_")
    }else{data_pre <-paste0(data_dir,  "/", dataset, "_")}
    
    dir.create(out_dir, showWarnings=F, recursive = T)
    dir.create(data_dir,showWarnings=F, recursive = T)
    
    
    ## if IDat files and has no sampleinfo , generate the sampleinf
    (file_name <- paste0(file_dir, 'SampleInfo.txt'))
    if(any(grepl('.idat', list.files(file_dir), ignore.case = T)) & !file.exists(file_name)){
        saveSampleInfo(file_dir, design_f)
    }
    
    #*****************************#
    # read raw data
    #*****************************#
    print(paste0("Load raw data of ", dataset),quote=F)
    setwd(file_dir)
    
    #*****************************#
    # read design and raw data for both iDat and other illumina
    #*****************************#
    
    ## read the design file and load the rawdata
    if(file.exists(file_name) & any(grepl('.idat', list.files(file_dir), ignore.case = T))){ # idat file has sampleinfo.txt, where other illumina dont
        ### IDAT files
        print('INPUT is idat files')
        array_design <- read.delim(file_name, comment.char = '#')
        (illumina_f <- as.character(array_design$FileName))
        bgx <- paste0(current_wd,'/', bgx)
        rawdata <- read.idat(illumina_f, bgxfile = bgx, dateinfo=T, tolerance=0)
        idat =T
    }else{
        ## grep the non normalized file
        array_design <- processDesign(dataset, design_f)
        (illumina_f <- grep("GSE.*N|normalized.*.gz", list.files(file_dir), value=T))
        rawdata <- read.ilmn(illumina_f,sep=sep)
        idat=F
        print('INPUT is non-normalized files')
    }
    
    
    if(idat){
        ## get the metadata(iDAT)
        df <- rawdata@.Data[3][[1]]
        writeTable(df, f_out = paste0(file_dir, 'scan_metadata.tsv'))
    }
    
    
    #*****************************#
    # Update the sample names to GSM based on SampleID (has ExternalID and SampleID)
    #*****************************#
    if(!idat){ ## map the sample names in the non idat
        (sampleid <- grep("SampleID.txt", list.files(file_dir), value=T))
        if(length(sampleid) >0 ){
            sampleid <- read.delim(sampleid, comment.char = '#')
            (sample_match<- as.character(sampleid$ExternalID))
            (names(sample_match) <- as.character(sampleid$SampleID))
        }else{
            df_tmp <- data.frame(SampleID = dimnames(rawdata$E)[2][[1]])
            writeTable(df_tmp, f_out = paste0(file_dir, 'sample_name_list.tsv'))
            stop(paste0('No sampleID found in ', file_dir))
        }
    }else{ ## for iDat sample names from file name
        (sample_match<- as.character(array_design$ExternalID))
        

        if(any(grep('\\.idat', as.character(array_design$FileName)))){## '.idat is in the filename, remove it
            (names(sample_match) <-gsub('\\.idat','',as.character(array_design$FileName)))
        } else{
            (names(sample_match) <-as.character(array_design$FileName))
            }
    }
    
    ## update the sample names in the raw data to GSM
    dimnames(rawdata$E)[2][[1]] <- mapValue(dimnames(rawdata$E)[2][[1]], sample_match)
    ## get the matrix sample names
    (sample_order <- dimnames(rawdata$E)[2][[1]])
    
    ## reorder by sample names (GSM id)
    rownames(array_design) <- array_design$ExternalID
    array_design <- array_design[sample_order, ]
    
    ## sample names for plot label
    if('label' %in% colnames(array_design)){
        (sample_names <- as.character(array_design$label))
    }else{
        (sample_names <- as.character(array_design$Sample))
    }
    
    
    ## set colour for samples
    cols <- noWarnings(brewer.pal(length(sample_names), "Set1"))
    
    #*****************************#
    # RAW: exp value plots for Raw data
    #*****************************#
    print("Plot rawdata",quote=F)
    ## plot the raw data
    #1. boxplot
    png(filename= paste0(plot_pre, 'boxplot_rawdata.png'), width = width, height = height)
    par(mar = c(25, 5, 4, 2)+ 0.1)
    boxplot(log2(rawdata$E), col= cols, names = sample_names,las=2,  ylab = "unnormalised intensity values", 
            main = paste0(dataset, ": raw data")) #plot boxplots for each sample
    dev.off()
    
    
    #*****************************#
    # Background correction and 
    # normalization and save matrix
    #*****************************#
    #The neqc functions performs normexp background correction using negative controls, then quantile
    # normalizes and finally log2 transforms
    
    
    ## normalize data by neqc
    print("NORMALIZE DATA")
    eset<-neqc(rawdata)
    array_dat <- eset$E
    
    #     ## sample to sample corr (normalized)
    #     corr_mat<-cor(array_dat)
    #     res<-apply(corr_mat,1,mean)
    
    ## output the normalized matrix
    (output <- paste0(data_pre, "neqc_normalized.data.tsv" ))
    
    # update the column names and probe names
    df <- as.data.frame(array_dat)
    (colnames(df) <- array_design$ExternalID)
    ## update the probe names to probe ID
    probes_id <- eset@.Data[2][[1]]
    rownames(df) <- probes_id$Probe_Id
    df <- cbind(data.frame(Probe = rownames(df)), df)
    
    
    tmp_msg <- paste0("#", Sys.Date(), "\n# ", dataset,
                      "\n# Background corrected and normalized by neqc, log2 transformed\n#")
    sink(output, type="output")
    writeLines(tmp_msg)
    sink()
    
    noWarnings(write.table(df, file = output, row.names = F,
                           sep ='\t', quote = F, append = T))
    print(paste0("Saved normalized expression file ", output))
    
    
    #*****************************#
    # normalized NEQC QC: plots and cluster
    #*****************************#
    print("plot normalized")
    #1. boxplot
    png(filename= paste0(plot_pre, 'boxplot_NEQC_normalized.png'), width = width, height = height)
    par(mar = c(25, 5, 4, 2)+ 0.1)
    boxplot(eset$E, col= cols, names = sample_names,las=2,  ylab = "NEQC normalized intensity values", 
            main = paste0(dataset, ": NEQC normalized")) #plot boxplots for each sample
    dev.off()
    
    #3. cluster
    dis_method <- "maximum"
    hclust_method <- "complete"
    distance <- dist(t(array_dat),method=dis_method)
    clusters <- hclust(distance, method = hclust_method)
    clusters$labels <- sample_names
    
    ## colour the genotypes
    if ("Genotype" %in% colnames(array_design)){
        (genotypes <- levels(as.factor(array_design$Genotype)))
        n <- nlevels(as.factor(array_design$Genotype))  # number of genotypes
        dend <- as.dendrogram(clusters)
        
        (genotype_col <- brewer.pal(n=7, "Set1")[1:n]) #specify colors of the genotype
        (names(genotype_col) <- genotypes)
        (label_col <- HelperFunctions::mapValue(array_design$Genotype, genotype_col))
        (names(label_col) <- sample_names) # rename the list by sample names
        
        # loading the package for branch color
        noWarnings(library(dendextend))
        # get sample names in the order of the tree
        (tree_order <- clusters$labels[clusters$order])
        # Assigning the labels of dendrogram object with new colors:
        (sample_col <- HelperFunctions::mapValue(tree_order, label_col, verbose=F))      
        
        labels_colors(dend) <- sample_col
        
        png(filename= paste0(plot_pre, 'cluster_neqc_normalized.png'), width = width, height = height)  
        par(mar = c(25, 5, 4, 2)+ 0.1)
        plot(dend, main = paste0(dataset, " cluster dendrogram: NEQC normalized\n distance method: ", dis_method,
                                 "; hierarchical cluster method: ", hclust_method))
        legend("topright", legend = genotypes, fill = genotype_col, title="Genotype", box.col="transparent",  cex=0.8)
        dev.off()
        
    } else {
        png(filename= paste0(plot_pre, 'cluster_neqc_normalized.png'), width = width, height = height)
        par(mar = c(25, 5, 4, 2)+ 0.1)
        plot(clusters, main = paste0(dataset, " cluster dendrogram: NEQC normalized\n distance method: ", dis_method,
                                     "; hierarchical cluster method: ", hclust_method))
        dev.off()
    }
    setwd(current_wd)
}