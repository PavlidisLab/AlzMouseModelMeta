#' Author: Beryl Zhuang
#' created: 2016-01-27
#' update: 2016-05-10
#' publish_codes
#' ################################################################################################
#' Quality control procedures of raw affy CEL files (modified from Sanja Rogic's script for FASD project )
#' 
#' INPUT: a folder dir with CEL files
#' OUTPUT: heatmap, boxplot, and a quality report 
#'  
#' NOTES:
#'    - required packages: simpleaffy, affyPLM and save RMA matrix
#' USAGE:
# source('preprocess/data_QC.R')
# file_dir <-'../AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE8425/' ## contains *.CEL files
# out_dir <- '../AD_mouse_model_project/data_and_QC/GEO_data/other_QC/'
# dir.create(out_dir, showWarnings=F)
# dataQC(file_dir, out_dir=out_dir)

####################
# preprocessing:
####################
#'1. 

#' @examples 
# source('preprocess/data_QC_affy.R')
# 
# file_dir <-'/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE14499/'
# design_f=paste0('../configs/AD_mouse_dataset_doc/shorten_experimental_design_sample_names.tsv')
# out_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/CEL_QC/')
# data_dir <- paste0(home_dir, '/AD_mouse_model_project/data_and_QC/GEO_data/normalized_matrix/')
# 
# 
# dataQC(i, design_f=design_f, out_dir = out_dir, data_dir =data_dir, result_dir =out_dir,chip_img = F, nuse_rle =T)



#' ----------------------------------------------------- # 
#' ################################################################################################



#########################
library(simpleaffy)
library(affyPLM)
library(RColorBrewer)
library(HelperFunctions) # helpers

 
source("helper_functions.R") ## for saveSampleInfo()

#' 
#' INPUT:
#'  - file_dir: e.g. '/home/bzhuang/AD_mouse_model_project/data_and_QC/GEO_data/CEL_raw/GSE13691.1/'
#'  - design_f
#'  - outdir: outdir for results
#'  -  chip_img: plot chip img for each sample, default = F
#' OUTPUT
#'  - heatmap, boxplot, and a quality report

###############################
#  main function
###############################
dataQC <- function(file_dir, design_f, 
                   out_dir = file_dir, width = 1600, height=1000, data_dir ="", result_dir ="", 
                   chip_img = F, nuse_rle =T, mk_sample_info =T){
    
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
    mainDataQC <- function(file_dir, out_dir, dataset, width, height, data_dir, result_dir, chip_img, nuse_rle, mk_sample_info){
        #'
        #' INPUT:
        #'  - file_dir: dir with all CEL files and "SampleInfo.txt"
        #'  - out_dir: out_dir for the QC results and plots
        #'  - dataset: dataset for the plot and file names
        #'  - chip_img: plot chip img for each sample
        #'  - nuse_rle: use the NUSE and RLE QA step (e.g. GPL6246: nuse_rle = F)
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

        
        #*****************************#
        # load files, and set working dir to CEL file dir
        #*****************************#
        print("Load raw data",quote=F)
        setwd(file_dir)
        rawdata<-read.affy(covdesc="SampleInfo.txt")
        
        ## sample names for plot label
        #(sample_names <- gsub(",", "\n", rawdata$label)) #break the names into multiple lines
        (sample_names <- as.character(rawdata$label))
        
        ## set colour for samples
        cols <- noWarnings(brewer.pal(length(sample_names), "Set1"))
        
        ## save the batch info
        batch_df <- data.frame(ExternalID=rawdata@phenoData@data$ExternalID, Batch_info =rawdata@protocolData$ScanDate)
        write.table(batch_df, file = "batchinfo.tsv", sep='\t',quote=FALSE, row.names =F )
        
        #*****************************#
        # RAW: exp value plots for Raw data
        #*****************************#
        print("Plot rawdata",quote=F)
        ## plot the raw data
        #1. boxplot
        png(filename= paste0(plot_pre, 'boxplot_rawdata.png'), width = width, height = height)
        par(mar = c(25, 5, 4, 2)+ 0.1)
        boxplot(rawdata, col= cols, names = sample_names,las=2,  ylab = "unnormalised intensity values", 
                main = paste0(dataset, ": raw data")) #plot boxplots for each sample
        dev.off()
        
        #2. intensity plot
        png(filename= paste0(plot_pre, 'intensity_rawdata.png'), width = width, height = height)
        hist(rawdata, col=cols, main = paste0(dataset, ": raw data"))
        dev.off()
        
        #*****************************#
        # RMA background correction and save matrix
        #*****************************#
        
        ## normalize data by RMA
        print("save normalized data")
        eset<-rma(rawdata, normalize=TRUE, background=TRUE)
        array_dat <- exprs(eset)
        
        ## get the design file
        array_design <- rawdata@phenoData@data
        
        ## sample to sample corr (normalized)
        corr_mat<-cor(exprs(eset))
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
        
        


        if(nuse_rle){ #---- this step does not work with GPL6246
            #*****************************#
            # Optional: RAW QC: NUSE and RLE  
            #*****************************#
            print("ALL NUSE+RLM ANALYSIS",quote=F)
            ## Perform probe-level metric calculations on the CEL files
            print("RAW data: perform probe level metric calculations on the CEL files")
            rawdata.plm<-fitPLM(rawdata)
            # affyPLM also provides more informative boxplots
            
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
            
            
            
            ## NUSE
            # We can also use NUSE (Normalised Unscaled Standard Errors). The median standard error should be 1 for most genes.
            nuse<-NUSE(rawdata.plm,type="stats")
            png(filename= paste0(plot_pre, 'NUSE_rawdata.png'), width = width, height = height)
            par(mar = c(25, 5, 4, 2)+ 0.1)
            NUSE(rawdata.plm, names = sample_names,las=2,  main = paste0(dataset, ": NUSE \n(Normalised Unscaled Standard Errors), raw data")) # for plotting boxplot
            dev.off()
            
            ## RLE (Relative Log Expression) plots should have values close to zero. to spot outlier
            rle<-RLE(rawdata.plm,type="stats")
            RLE(rawdata.plm) # for plotting boxplot
            png(filename= paste0(plot_pre, 'RLE_rawdata.png'), width = width, height = height)
            par(mar = c(25, 5, 4, 2)+ 0.1)
            RLE(rawdata.plm, names = sample_names,las=2, main = paste0(dataset, ": RLE \n(Relative Log Expression), raw data")) # for plotting boxplot
            dev.off()
            
            #*****************************#
            # RAW QC: more QC
            #*****************************#
            ## the 3'/5' ratio for RNA quality
            ## the scale factor
            ## average background
            
            print("ALL QC ANALYSIS",quote=F)
            
            #quality analysis III - all arrays - looking at the QC statistics directly
            qc<-qc(rawdata)
            
            ## the 3'/5' ratio for RNA quality
            ratio<-as.data.frame(ratios(qc))
            ratio.score<-apply(ratio,1,function(x){if(x["actin3/actin5"]>3 && x["gapdh3/gapdh5"]>1.25) 1 else 0})
            
            ## the scale factor
            scale<-qc@scale.factors
            avg<-mean(scale)
            scale.score<-as.numeric(scale>3*avg | scale<avg/3)
            
            ## average background
            background<-qc@average.background
            avg<-mean(background)
            background.score<-as.numeric(background>1.5*avg | background<avg/1.5)
            
            ## percent present 
            present<-qc@percent.present
            avg<-mean(present)
            present.score<-as.numeric(present>1.2*avg | present<avg/1.2)
            
            ## check the RNA degradation and plot
            rnadeg<-AffyRNAdeg(rawdata)
            png(filename= paste0(plot_pre, 'RNA_degradation_rawdata.png'), width = width, height = height)
            plotAffyRNAdeg(rnadeg, col=cols)
            dev.off()
            sumdeg<-data.frame(t(summaryAffyRNAdeg(rnadeg)))
            colnames(sumdeg) <- c("RNA_degradation_slope", "RNA_degradation_p_value")
            sumdeg$RNA_deg_outlier <- "No"
            sumdeg$RNA_deg_outlier[which(sumdeg[,1]>mean(sumdeg[,1])+2*sd(sumdeg[,1]))] <- "Yes"
            
            #*****************************#
            # Report for all scores
            #*****************************#
            
            ## report on 
            #' Raw data: NUSE IQR.NUSE    RLE IQR.RLE, ratio, scale, background, present, sumdeg
            #' RMA normalized: CORR
            #' see FASD sup for threshold detail for the scores
            #' 
            QA<-data.frame(QA_NUSE=round(nuse[1,],digits=3),QA_IQR_NUSE=round(nuse[2,],digits=3),
                           QA_RLE=round(rle[1,],digits=3),QA_IQR_RLE=round(rle[2,],digits=3),
                           QA_sample_corr_normalized=round(res,digits=3),row.names=row.names(pData(eset)))
            
            QA_score<-apply(QA,1,function(x){count=0;if(x["QA_NUSE"]>1.05) count<-count+1;
                                             if(x["QA_IQR_NUSE"]>2*mean(nuse[2,])) count<-count+1; 
                                             if( abs(x["QA_RLE"])>0.05) count<-count+1;
                                             if(x["QA_IQR_RLE"]>2*mean(rle[2,])) count<-count+1;
                                             if(x["QA_sample_corr_normalized"]<0.9) count<-count+1;count})
            
            QA<-cbind(ExternalID=rawdata$ExternalID, Sample=rawdata$Sample, QA,QA_score)
            
            ## combine the QC info
            QA<-cbind(QA, ratio, scale, background, present, sumdeg)
            
            all<-ratio.score+scale.score+present.score+background.score
            
            qc.report<-data.frame(RATIO=ratio.score,SCALE=scale.score,BACKGROUND=background.score,PRESENT=present.score,TOTAL_SCORE=all,row.names=row.names(ratio))
            QA<-cbind(data.frame(Dataset = rep(dataset, nrow(QA))), data.frame(File = rownames(QA)), QA,qc.report)
            
            write.table(QA,file=paste0(result_pre, "QA_QC_report.tsv"),sep='\t',quote=FALSE, row.names =F)
        }else{print("Skipped NUSE, RLE and other QC")}
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
    mainDataQC(file_dir, out_dir, dataset, width, height, data_dir, result_dir, chip_img, nuse_rle)

    

    setwd(current_wd)


}









