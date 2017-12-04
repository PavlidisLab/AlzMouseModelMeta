#' 2016-06-07
#' Process illumnia  non normalized data by limma package
#' data is assumed background corrected from beadarray studio.
#' this script will log2 transform and quantile normalization

#GPL6333



 
source("helper_functions.R")

#'@examples
# file_dir <- "/home/bzhuang/HD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE19676/"
# design_f='/home/bzhuang/HD_mouse_model_project/data_and_QC/design_and_annotation/shorten_experimental_design_sample_names.tsv'
# 
# out_dir <- "/home/bzhuang/HD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/"  # plots
# data_dir <- "/home/bzhuang/HD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/" 
# dataIlluminaQC(file_dir, design_f, 
#                            out_dir = out_dir, width = 1600, height=1000, 
#                            data_dir =data_dir, result_dir =out)


##########

if('package:oligo' %in% search()){
detach('package:oligo', unload = TRUE, character.only = T) 
}

library(lumi)
 
source("helper_functions.R")
source('preprocess/data_QC.R')


dataIlluminaQC <- function(file_dir, design_f, 
                            out_dir = file_dir, width = 1600, height=1000, 
                           data_dir ="", result_dir =""){
    
    dir.create(out_dir, recursive = T,showWarnings = F)
    dataset <- getGSEID(file_dir)
    ###############################
    #  main QC function from Sanja
    ###############################
    maindataIlluminaQC <- function(file_dir, out_dir, dataset, width, height, data_dir, result_dir){
        #'
        #' INPUT:
        #'  - file_dir: dir with all CEL files and "SampleInfo.txt"
        #'  - out_dir: out_dir for the QC results and plots
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
        
        dir.create(out_dir, showWarnings=F,recursive = T)
        dir.create(data_dir,showWarnings=F,recursive = T)
        
        print("Load illumina non-normalized data",quote=F)
        setwd(file_dir)
        (illumina_f <- grep("GSE.*N|normalized.*.gz", list.files(file_dir), value=T))
        rawdata <- lumiR.batch(illumina_f)
        
        
        ## get sample design info, and match with the sample names in the array
        df <- read.delim('SampleID.txt', comment.char = '#')
        df <- noWarnings(left_join(data.frame(SampleID = sampleNames(rawdata)), df))
        array_design <- processDesign(dataset, design_f)
        array_design$label <- paste0(array_design$ExternalID,',',array_design$Sample)
        array_design <- noWarnings(left_join(df, array_design))
        
        rawdata@phenoData@data <- array_design
        sampleNames(rawdata) <- array_design$label
        ## sample names for plot label
        (sample_names <- as.character(array_design$label))
        
        ## set colour for samples
        cols <- noWarnings(brewer.pal(length(sample_names), "Set1"))
        
        #*****************************#
        # RAW: exp value plots for Raw data (not background corrected, not normalized)
        #*****************************#
        print("Plot rawdata",quote=F)
        ## plot the raw data
        #1. boxplot
        png(filename= paste0(plot_pre, 'boxplot_rawdata.png'), width = width, height = height)
        par(mar = c(25, 5, 4, 2)+ 0.1)
        boxplot(log2(exprs(rawdata)), col= cols, names = sample_names,las=2,  ylab = "unnormalised intensity values", 
                main = paste0(dataset, ": raw data")) #plot boxplots for each sample
        dev.off()
        
        #2. intensity plot
        png(filename= paste0(plot_pre, 'intensity_rawdata.png'), width = width, height = height)
        density(rawdata, col=cols, main = paste0(dataset, ": raw data"))
        dev.off()
        
        
        #*****************************#
        # Illumina assume background is corrected from beadstudio output
        #*****************************#
        
        ## We suppose the BeadStudio output data has been background corrected.
        
#         ##Variance stabilizing transform
#         Variance stabilization is critical for subsequent statistical inference to identify
#         differential genes from microarray data. We devised a variance-stabilizing trans-
#             formation (VST) by taking advantages of larger number of technical replicates
#         available on the Illumina microarray. 
        
        
        print("save normalized data")
        eset<-lumiT(rawdata,method = 'log2')
        ## Do quantile between microarray normaliazation
        eset <- lumiN(eset)
        array_dat <- as.data.frame(exprs(eset))
        


        
        ## sample to sample corr (normalized)
        corr_mat<-cor(array_dat)
        res<-apply(corr_mat,1,mean)
        
        ## output the normalized matrix
        output <- paste0(data_pre, "illumina_normalized.data.tsv" )
        
        df <- as.data.frame(array_dat)
        (colnames(df) <- rawdata$ExternalID)
        df <- cbind(data.frame(Probe = rownames(df)), df)
        
        tmp_msg <- paste0("#", Sys.Date(), "\n# ", dataset,
                          "\n# quantile normalized, log2 transformed\n#")
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
        boxplot(array_dat, col= cols, names = sample_names,las=2,  ylab = "RMA normalized intensity values", 
                main = paste0(dataset, ": quantile normalized")) #plot boxplots for each sample
        dev.off()
        
        #2. intensity plot
        png(filename= paste0(plot_pre, 'intensity_rma_normalized.png'), width = width, height = height)
        density(eset, col=cols, main = paste0(dataset, ": quantile normalized"))
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
    }
    
    ###############################
    #  main Function
    ###############################
    mainFunction <- function(file_dir, out_dir, design_f, width, height, data_dir, result_dir){
        ## make dirs
        dir.create(out_dir, showWarnings=F)
        dir.create(data_dir, showWarnings=F)
        dir.create(result_dir, showWarnings=F)
        
        
        dataset <- getGSEID(file_dir)
        print(paste0("Processing ", dataset))
        maindataIlluminaQC(file_dir, out_dir, dataset, width, height, data_dir, result_dir)
    }
    
    mainFunction(file_dir, out_dir, design_f, width, height, data_dir, result_dir)
}








