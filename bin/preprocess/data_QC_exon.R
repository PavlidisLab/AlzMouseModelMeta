#' Author: Beryl Zhuang
#' created: 2016-06-30
#' take template from data_QC, specialized for affy exon arrays, only do rle QC
#' for platform: GPL6096
#' #' publish_codes
#' ################################################################################################
#' Quality control procedures of raw affy CEL files (modified from Sanja Rogic's script for FASD project )
#' 
#' INPUT: a folder dir with CEL files, exon arrays
#' OUTPUT: heatmap, boxplot, and a quality report 
#'  
#' NOTES:
#'    - required packages: oligo
#' USAGE:
# source('preprocess/data_QC_exon.R')
# file_dir <-'/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE50521/' ## contains *.CEL files
# out_dir <- '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/'
# dir.create(out_dir, showWarnings=F)
# dataExonArrayQC(file_dir, out_dir=out_dir)


#        eset<-oligo::rma(rawdata,background=TRUE, normalize=TRUE, target = "core")

####################
# preprocessing:
####################
#'1. 

#' @examples 
# source('preprocess/data_QC.R')
# 
# file_dir <-'/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE8425/'
# 
# design_f='./data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv'
# out_dir <- paste0('./data_and_QC/GEO_data/CEL_QC/')
# data_dir <- paste0('./data_and_QC/GEO_data/normalized_matrix/')
# 
# 
# dataExonArrayQC(file_dir, design_f=design_f, out_dir = out_dir, data_dir =data_dir, result_dir =out_dir, chip_img = F, rle =T)



#' ----------------------------------------------------- # 
#' ################################################################################################



#########################
library(oligo)
library(RColorBrewer)
library(HelperFunctions) # helpers

 
source("helper_functions.R")

#' 
#' INPUT:
#'  - file_dir: e.g. '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE13691.1/'
#'  - design_f
#'  - outdir: outdir for results
#'  -  chip_img: plot chip img for each sample, default = F
#' OUTPUT
#'  - heatmap, boxplot, and a quality report



dataExonArrayQC <- function(file_dir, design_f, 
                            out_dir = file_dir, width = 1600, height=1000, data_dir ="", result_dir ="", 
                            chip_img = F, rle =T, mk_sample_info =T){
    
    dir.create(out_dir, recursive = T,showWarnings = F)
    
    dir.create(out_dir, recursive = T,showWarnings = F)
    current_wd <- getwd()
    
    #*****************************#
    ## make the sample info 
    # need SampleInfo.txt in the CEL file dir to read the files
    #*****************************#
    if(mk_sample_info){
        saveSampleInfo(file_dir, design_f)
    }
    
    
    ###############################
    #  main QC function from Sanja
    ###############################
    maindataExonArrayQC <- function(file_dir, out_dir, dataset, width, height, data_dir, result_dir, chip_img, rle, mk_sample_info){
        #'
        #' INPUT:
        #'  - file_dir: dir with all CEL files and "SampleInfo.txt"
        #'  - out_dir: out_dir for the QC results and plots
        #'  - dataset: dataset for the plot and file names
        #'  - chip_img: plot chip img for each sample
        #'  - rle: use the NUSE and RLE QA step (e.g. GPL6246: rle = F)
        #'  OUTOUT:
        #'  - heatmap, boxplot, and a quality report, and save RMA matrix
        
        #*****************************#
        # set parameters
        #*****************************#
        (plot_pre <- paste0(out_dir,  "/", dataset, "_"))
        if(result_dir ==""){
            result_pre <-paste0(out_dir,  "/", dataset, "_")
        }else{result_pre <-paste0(result_dir,  "/", dataset, "_")}
        if(data_dir ==""){
            data_pre <-paste0(out_dir,  "/", dataset, "_")
        }else{data_pre <-paste0(data_dir,  "/", dataset, "_")}
        
        dir.create(out_dir, showWarnings=F)
        dir.create(data_dir,showWarnings=F)
        
        print("Load raw data",quote=F)
        setwd(file_dir)
        exonCELs <- list.celfiles(file_dir, full.names = TRUE)
        rawdata <- read.celfiles(exonCELs)
        
        ## get sample design info
        df <- read.delim('SampleInfo.txt', comment.char = '#')
        df <- left_join(data.frame(FileName = sampleNames(rawdata)), df)
        rownames(df) <- df$FileName
        rawdata@phenoData@data <- df
        ## sample names for plot label
        #(sample_names <- gsub(",", "\n", rawdata$label)) #break the names into multiple lines
        (sample_names <- as.character(df$label))
        
        ## set colour for samples
        cols <- noWarnings(brewer.pal(length(sample_names), "Set1"))
        
        ## save the batch info
        batch_df <- data.frame(ExternalID=rawdata@phenoData@data$ExternalID, 
                               Batch_info =rawdata@protocolData@data[, 'dates'])
        write.table(batch_df, file = "batchinfo.tsv", sep='\t',quote=FALSE, row.names =F )
        
        
        
        #*****************************#
        # RAW: exp value plots for Raw data (not background corrected, not normalized)
        #*****************************#
        eset_not_normalized <-oligo::rma(rawdata,background=F, normalize=F, target = "core")
        print("Plot rawdata",quote=F)
        ## plot the raw data
        #1. boxplot
        png(filename= paste0(plot_pre, 'boxplot_rawdata.png'), width = width, height = height)
        par(mar = c(25, 5, 4, 2)+ 0.1)
        boxplot(eset_not_normalized, col= cols, names = sample_names,las=2,  ylab = "unnormalised intensity values", 
                main = paste0(dataset, ": raw data")) #plot boxplots for each sample
        dev.off()
        
        #2. intensity plot
        png(filename= paste0(plot_pre, 'intensity_rawdata.png'), width = width, height = height)
        hist(eset_not_normalized, col=cols, main = paste0(dataset, ": raw data"))
        dev.off()
        
        
        #*****************************#
        # RMA background correction and save matrix
        #*****************************#
        
        ## background correction and then normalize data by RMA
        print("save normalized data")
        eset<-oligo::rma(rawdata,background=TRUE, normalize=TRUE, target = "core")
        array_dat <- as.data.frame(exprs(eset))
        
        ## get the design file
        array_design <- rawdata@phenoData@data
        
        ## sample to sample corr (normalized)
        corr_mat<-cor(array_dat)
        res<-apply(corr_mat,1,mean)
        
        ## output the normalized matrix
        output <- paste0(data_pre, "RMA_normalized.data.tsv" )
        
        df <- as.data.frame(array_dat)
        (colnames(df) <- rawdata$ExternalID)
        df <- cbind(data.frame(Probe = rownames(df)), df)
        
        tmp_msg <- paste0("#", Sys.Date(), "\n# ", dataset,
                          "\n# RMA normalized, log2 transformed\n#")
        sink(output, type="output")
        writeLines(tmp_msg)
        sink()
        
        noWarnings(write.table(df, file = output, row.names = F,
                               sep ='\t', quote = F, append = T))
        print(paste0("Saved normalized expression file ", output))
        
        
        #*****************************#
        # normalized RMA QC: plots and cluster
        #*****************************#
        print("plot RMA normalized")
        #1. boxplot
        png(filename= paste0(plot_pre, 'boxplot_rma_normalized.png'), width = width, height = height)
        par(mar = c(25, 5, 4, 2)+ 0.1)
        boxplot(eset, col= cols, names = sample_names,las=2,  ylab = "RMA normalized intensity values", 
                main = paste0(dataset, ": RMA normalized")) #plot boxplots for each sample
        dev.off()
        
        #2. intensity plot
        png(filename= paste0(plot_pre, 'intensity_rma_normalized.png'), width = width, height = height)
        hist(eset, col=cols, main = paste0(dataset, ": RMA normalized"))
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
            n <- nlevels(as.factor(rawdata@phenoData@data$Genotype))  # number of genotypes
            dend <- as.dendrogram(clusters)
            
            (genotype_col <- brewer.pal(n=7, "Set1")[1:n]) #specify colors of the genotype
            (names(genotype_col) <- genotypes)
            (label_col <- HelperFunctions::mapValue(array_design$Genotype, genotype_col))
            (names(label_col) <- array_design$label) # rename the list by sample names
            
            # loading the package for branch color
            noWarnings(library(dendextend))
            # get sample names in the order of the tree
            (tree_order <- clusters$labels[clusters$order])
            # Assigning the labels of dendrogram object with new colors:
            (sample_col <- HelperFunctions::mapValue(tree_order, label_col, verbose=F))      
            
            labels_colors(dend) <- sample_col
            
            png(filename= paste0(plot_pre, 'cluster_rma_normalized.png'), width = width, height = height)  
            par(mar = c(25, 5, 4, 2)+ 0.1)
            plot(dend, main = paste0(dataset, " cluster dendrogram: RMA normalized\n distance method: ", dis_method,
                                     "; hierarchical cluster method: ", hclust_method))
            legend("topright", legend = genotypes, fill = genotype_col, title="Genotype", box.col="transparent",  cex=0.8)
            dev.off()
            
        } else {
            png(filename= paste0(plot_pre, 'cluster_rma_normalized.png'), width = width, height = height)
            par(mar = c(25, 5, 4, 2)+ 0.1)
            plot(clusters, main = paste0(dataset, " cluster dendrogram: RMA normalized\n distance method: ", dis_method,
                                         "; hierarchical cluster method: ", hclust_method))
            dev.off()
        }
        
        
        if(rle){
            #*****************************#
            # Optional: RAW QC: NUSE and RLE  
            #*****************************#
            print("ALL RLM ANALYSIS",quote=F)
            ## Perform probe-level metric calculations on the CEL files
            print("RAW data: perform probe level metric calculations on the CEL files")
            rawdata.plm<-fitProbeLevelModel(rawdata,background=TRUE, normalize=TRUE, target="core")
            
            #*****************************#
            # RAW QC: loop to create image for each sample
            #*****************************#
            ## creat images for each sample # only work at lab computer...
            if(chip_img){
                print("Plot chip images for each experiment",quote=F)
                for(i in 1: length(rawdata$label)){
                    png(filename= paste0(plot_pre, rawdata$label[i], '_rawdata.png'), width = width, height = height)
                    image(rawdata.plm, which=i, add.legend=TRUE)
                    dev.off()
                }
            }
            
            
            ## RLE (Relative Log Expression) plots should have values close to zero. to spot outlier
            rle<-RLE(rawdata.plm,type="stats")
            RLE(rawdata.plm) # for plotting boxplot
            png(filename= paste0(plot_pre, 'RLE_rawdata.png'), width = width, height = height)
            par(mar = c(25, 5, 4, 2)+ 0.1)
            RLE(rawdata.plm, names = sample_names,las=2, main = paste0(dataset, ": RLE \n(Relative Log Expression), raw data")) # for plotting boxplot
            dev.off()
            
        }else{print("skipped RLE and chip img")}
    }
    
    ############################################################
    
    #3':5' intensity ratio (recommended by Affy)
    #ratio<-as.data.frame(ratios(qc)) #3:5 ratios for all samples
    #ratio[ratio$"actin3/actin5">3 & ratio$"gapdh3/gapdh5">1.25,]
    
    #for which samples the scale factor is not within 3-fold from the mean (affy recommendation)
    #scale<-qc@scale.factors
    #which(scale>3*mean(scale) | scale<mean(scale)/3)
    
    #background (should agree between samples - no specific recommendation)
    #background<-qc@average.background
    #which(background>1.5*mean(background) | background<mean(background)/1.5)
    
    #percentage of probes called present; should be in agreement between the samples; no specific recommendation from Affy
    
    #present<-qc@percent.present
    #which(present>1.2*mean(present) | present<mean(present)/1.2)
    
    
    ###############################
    #  main Function
    ###############################

    ## make dirs
    dir.create(out_dir, showWarnings=F)
    dir.create(data_dir, showWarnings=F)
    dir.create(result_dir, showWarnings=F)
    
    
    dataset <- getGSEID(file_dir)
    print(paste0("Processing ", dataset))
    maindataExonArrayQC(file_dir, out_dir, dataset, width, height, data_dir, result_dir, chip_img, rle)
    
    setwd(current_wd)

}








