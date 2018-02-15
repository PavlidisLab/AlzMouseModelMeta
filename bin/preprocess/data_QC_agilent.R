# need more doc
# 2016-02-05
# updated 2016-06-29
# this script is to process two coloured agilent platforms.
# background are corrected by "normexp", "rma" method. offset = 50. matrix is log transformed and quantile normalized.
# reads only green channel


#' OFFSET ( default offset = 50)
#' The offset can be used to add a constant to the intensities before log-transforming, 
#' so that the log-ratios are shrunk towards zero at the lower intensities. This may eliminate 
#' or reverse the usual 'fanning' of log-ratios at low intensities associated with local background subtraction. 
#' The higher the offset value, the value are pushed closer to 0 ( the density curve is less spread)
#' 

#' 2016-06-29
#' choose input from different channels: green(cy3) only, red only(cy5) output are expression, if both, output is ratio
#' channel: 

## need to detach oligo package (some function names overlapped)
if('package:oligo' %in% search()){
    if('package:pd.moex.1.0.st.v1' %in% search()){
        detach("package:pd.moex.1.0.st.v1", unload=TRUE)
    }
    detach('package:oligo', unload = TRUE, character.only = T) 
}

library(limma)
 
source("helper_functions.R") ## for saveSampleInfo()



#' @param file_dir file dir with all the raw files e.g. '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE13691.1/'
#' @param f_target (str) name of the file contain the list of sample dirs and exp design
#'          e.g. "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/other_raw/GSE63617.1/SampleInfo.txt.1.tsv"
#' @param out_dir (str) dir for plots
#' @param data_dir (str) dir for output matrix
#' @param green_only: if it's one colored, green_only =T, if it's two coloured, green_only =F, default = F (two color)
#' @param return_value: to return dataframes
#' @note
#' the dir has a file (f_target) with a list of samples dirs, and other experimental design files
#' @examples
# source("preprocess/data_QC_agilent.R")
# file_dir <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/other_raw/GSE63617.1/"
# f_target <-"/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/other_raw/GSE63617.1/SampleInfo.txt.1.tsv"
# out_dir <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/other_QC/"  # plots
# data_dir <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/other_raw/other_normalized_in_progress/"    # matix
# aligentDataQC(file_dir, f_target = f_target,agilent_source="agilent.median", channel='both',
#               out_dir = out_dir, width = 1600, height=1000, data_dir =data_dir)


aligentDataQC <- function(file_dir, f_target, agilent_source="agilent.median", 
                          channel=c('both','green', 'red'),
                   out_dir = file_dir, width = 1600, height=1000, data_dir ="", return_value =F){
    current_wd <- getwd()
    #*****************************#
    # set parameters
    #*****************************#
    setwd(file_dir)
    channel=channel[1]
    ## get only the root file name (input can be a full file dir)
    tmp <- unlist(strsplit(f_target, split = "/"))
    (f_target <- tmp[length(tmp)])
    
    ## check if the file list is in the same dir
    if(all(grepl(f_target, list.files(file_dir), fixed = T))){
        stop(paste0(f_target, " is not found in ", file_dir))
    }
    
    (dataset <- getGSEID(file_dir))
    print(paste0("Input file list: ", f_target))
    
    
    ## save dirs
    plot_pre <- paste0(out_dir,  "/", dataset, "_")
    if(data_dir ==""){
        data_pre <-paste0(out_dir,  "/", dataset, "_")
    }else{data_pre <-paste0(data_dir,  "/", dataset, "_")}
    
    dir.create(out_dir, showWarnings=F,recursive = T)
    dir.create(data_dir,showWarnings=F,recursive = T)
    
    #*****************************#
    # read raw data
    #*****************************#
    print(paste0("READING RAW DATA"))
    targets<- readTargets(f_target)
    
    ## the green (Cy3) channel only be read for one channel, for two channels green.only =F
    ## source="agilent.median" use median foreground;  Background estimates are always medians.
    ## if want to read cy5 (red channel), then both channels must be read
    if(channel =='green'){
        green_only =T
    }else{green_only =F}
    
    rawdata<- read.maimages(targets, source=agilent_source, green.only=green_only)
    
    (sample_names <- as.character(rawdata$targets$label))
    
    ## set colour for samples
    cols <- noWarnings(brewer.pal(length(sample_names), "Set1"))
    
    #*****************************#
    # plot raw data
    #*****************************#
    ## plot the raw data
    #1. boxplot (raw)
    png(filename= paste0(plot_pre, 'boxplot_rawdata.png'), width = width, height = height)
    par(mar = c(25, 5, 4, 2)+ 0.1)
    boxplot(as.matrix(rawdata), col= cols, names = sample_names,las=2,  
            ylab = "unnormalised intensity values", 
            main = paste0(dataset, ": background corrected")) #plot boxplots for each sample
    dev.off()
    
    #2. intensity plot
    png(filename= paste0(plot_pre, 'intensity_rawdata.png'), width = width, height = height)
    densPlotMultiLines(na.omit(as.matrix(rawdata)), xlab = "intensity", 
                       main = paste0(dataset, ": rawdata"))
    dev.off()
    
    
    #*****************************#
    # background correction, normalization, 
    # and calculate the average intensity (aggregate multiple readings of the same probe )
    #*****************************#
    ## Subtract the background signal.use RMA for normalization (not log2 transformed yet)
    ## offset and method is recommended by limma user guide
    raw_BGcorrected_norm <- limma::backgroundCorrect(rawdata, method="normexp", offset=50, 
                                              normexp.method = "rma")
    
    ## Then normalize and log-transformed the data.
    corrected_norm_pre <- limma::normalizeBetweenArrays(raw_BGcorrected_norm,method="quantile")
    
    ## Finally calculate the average intensity values from the probes of each gene.   
    corrected_norm <- avereps(corrected_norm_pre,ID=corrected_norm_pre$genes$ProbeName)

    
    ## if only want to extract the red (cy5) channel expression only
    if(channel == 'red'){
        corrected_norm_r <- normalizeBetweenArrays(raw_BGcorrected_norm$R,method="quantile")
        corrected_norm_r <- avereps(corrected_norm_r,ID=corrected_norm_pre$genes$ProbeName)
        ## log2 transformation
        if(max(range(corrected_norm_r)) >100){
            corrected_norm_r <- log2(corrected_norm_r)
        }
        ## get the expression matrix from red channel
        array_dat <- as.matrix(corrected_norm_r)
    }else{
        ## get the expression matrix
        array_dat <- as.matrix(corrected_norm)
    }

    

    #*****************************#
    # array_dat and design
    #*****************************#

    
    ## get the design file
    array_design <- rawdata$targets
    
    #*****************************#
    # RMA plots and cluster
    #*****************************#
    #1. boxplot
    png(filename= paste0(plot_pre, 'boxplot_rma_normalized.png'), width = width, height = height)
    par(mar = c(25, 5, 4, 2)+ 0.1)
    boxplot(array_dat, col= cols, names = sample_names,las=2,  ylab = "RMA normalized intensity values", 
            main = paste0(dataset, ": RMA normalized")) #plot boxplots for each sample
    dev.off()
    
    #2. intensity plot
    png(filename= paste0(plot_pre, 'intensity_rma_normalized.png'), width = width, height = height)
    densPlotMultiLines(array_dat, xlab = "intensity", 
                       main = paste0(dataset, ": background corrected, RMA normalized"))
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
        (n <- length(genotypes))  # number of different genotypes
        dend <- as.dendrogram(clusters)
        
        (genotype_col <- brewer.pal(n=7, "Set1"))[1:n] #specify colors of the genotype
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
    
    
    
    #*****************************#
    # output the normalized matrix
    #*****************************#    
    output <- paste0(data_pre, "RMA_normalized.data.tsv" )
    
    df <- as.data.frame(array_dat)
    (colnames(df) <- rawdata$targets$ExternalID)
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
    # Plot MA plots for raw and processed
    #*****************************# 
    
#     # Now some MA plots of only one replicate per condition:
#     noWarnings(library(affy))
#     png(filename= paste0(plot_pre, 'MA_raw.png'), width = width, height = height)
#     mva.pairs(na.omit(as.matrix(rawdata))) # Before BG correction
#     dev.off()
#     png(filename= paste0(plot_pre, 'MA_BG_corrected.png'), width = width, height = height)
#     mva.pairs(na.omit(as.matrix(raw_BGcorrected_norm))) # bg corrected by rma
#     dev.off()
#     png(filename= paste0(plot_pre, 'MA_rma_normalized.png'), width = width, height = height)
#     mva.pairs(na.omit(as.matrix(corrected_norm))) # rma normalized, quantile normalized
#     dev.off()
    
    
    #*****************************#
    # return values (optional)
    #*****************************#
    if (return_value){
        returnlist <- list(rawdata, raw_BGcorrected_norm, corrected_norm)
        print("A list of values are returned: 1. rawdata, 2. BG_corrected, 3. quantile_normalized")
        return(returnlist)
    }
    setwd(current_wd)
}







