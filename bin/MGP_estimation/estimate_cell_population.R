#' 2016-02-16
#' update 2017-01-16, 2017-02-16 added library
#' new pacakge installed from Ogan; works locally now
#' #' 
#' devtools::install_github('oganm/homologene')
#' devtools::install_github('oganm/ogbox')
#' devtools::install_github('oganm/markerGeneProfile')
#' use wrapper
#' ################################################################################################
library(dplyr)
library(ogbox)
library(markerGeneProfile)
# this function is a generic function that looks for the most variable probeset
# of a gene. unlike the previous one, it takes in objects and outputs objects 
#' from Ogan https://github.com/oganm/neuroExpressoAnalysis/blob/master/R/mostVariable.R
mostVariable = function(allDataPre,genes = 'Gene.Symbol', threshold = 6, threshFun = max){
    list[,exprData]= sepExpr(allDataPre)
    rowmax = apply(exprData, 1, threshFun)
    discludeGenes = (rowmax<threshold)
    allDataPre = allDataPre[!discludeGenes,]
    exprData = exprData[!discludeGenes,]
    
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    if (class(allDataPre)[1]=='data.table'){
        allDataPre = allDataPre[!duplicated(allDataPre[,genes, with=F]),]
    } else {
        allDataPre = allDataPre[!duplicated(allDataPre[,genes]),]
        
    }
    allDataPre = allDataPre[!allDataPre[,genes]=='',]
    return(allDataPre)
}




#' to process data exp, rm probes don't map to a gene, select the most variable probeset for the gene
#' @param r_ob (dir) R object dir contains array_dat, array_design and annotation for loading
#' @return a list: array_dat (df), array_design (df), dataset (str)
#' 
preProcessMat <- function(r_ob, probename = T, threshold = 6){
    #' probename =T input expression is by probes (F: input is by genes, i.e mixed model input)
    #' threshold" filter probes/genes with expression lower than 6
    load(r_ob)
    dataset <- levels(array_design$Dataset)
    
    ## must make sure the matrix sample order is the same as the design sample order!!
    array_dat <- array_dat[, as.character(array_design$Sample)]
    
    
    ### remove rows with no variation
    # input an array_dat, and remove rows with no variation (ie all readings are the same)
    # if not removed, there will be problem with PCA
    
    
    #     
    #     # set threshold by the min(array_matrix), for affy platforms threshold = 6, for agilent = -10(setting is)
    #     if(min(array_dat) <0){
    #         threshold = -10
    #     }else{
    #         threshold = 6
    #     }
    
    
    # add dataset id to the sample names in array and design if sample name dont have GSE id
    if(!all(grepl('GSE', array_design$Sample))){
        array_design$Sample <- paste0(dataset, "_", array_design$Sample)
        colnames(array_dat) <- paste0(dataset, "_", colnames(array_dat))
    }else{
        print('sample name contains GSE id, continue...')
    }
    
    
    # add probename and gene symbol to array_dat if expression by value
    if(probename){
        df <- data.frame(ProbeName = rownames(array_dat))
        df <- noWarnings(left_join(df, annotation[, 1:2]))
        array_dat <- cbind(df, array_dat)
    }else{
        #GeneSymbols to array_dat
        df <- data.frame(GeneSymbols=row.names(array_dat))  ## the genesymbols must be the first column to make it work
        array_dat <- cbind(df, array_dat)
        array_dat <- na.omit(array_dat)
    }
    
    
    # rm probes not mapped to a gene
    array_dat <- excludeMatch(array_dat, "GeneSymbols", "")
    
    # rm probes mapped to multiple genes
    index <- grep('\\|', array_dat$GeneSymbols)
    if(length(index >0 )){array_dat <- array_dat[-index, ]}
    
    # mostVariable, filter the low exp probes (threshold = 6), function from ogan
    # if multiple probes are mapped to the same gene, the most variable probset will be selected 
    # as the exp of the gene, threshold filter out low exp probes
    array_dat <- mostVariable(array_dat, genes = "GeneSymbols", threshold = threshold)
    array_dat <- array_dat[, setdiff(colnames(array_dat), 'ProbeName')]
    
    returnlist <- list(array_dat, array_design, dataset, threshold)
    return(returnlist)
}



cellPopAmountTest <- function (estimates, groups, out_f, sigTest = wilcox.test, wt_only = T){
    #' with the result of cellTypeEstimate, get teh result of cell type estimates (scale to 0-1) and test results
    #' wt only, only get results that compared with wt
    #' input non-scaled estimate
    ## the groups
    
    estimate_non_scaled = estimates
    
    
    (groupNames = as.character(unique(groups[[1]])))
    ## make sure WT is the last
    if('WT' %in%groupNames){
        groupNames = c(setdiff(groupNames, 'WT'), 'WT')
    }
    (comparisons = combn(groupNames, 2))
    
    if (typeof(estimates) != "list") {
        estimates = list(estimates)
    }
    estimates %<>% lapply(scale01)  ## all scale to 0-1 (for plot)
    
    df_all <- data.frame(sample = NULL, cell_type= NULL, PC1= NULL, PC1_scaled = NULL, group = NULL, adj_p = NULL)  ## for plot cell amount estimate
    comp_p = data.frame(cell_type=NULL, group=NULL, group_comp=NULL, raw_p=NULL) ## for test adj pvalues
    
    ## get the wilcox test result for each comparison
    if (!is.null(comparisons)) {
        for (i in 1:len(estimates)) { ## loop for each cell type
            (cell_type= names(estimates)[i])
            
            ## estimate for each cell type
            (frame = data.frame(cell_type= names(estimates)[i], PC1 = estimate_non_scaled[[i]], PC1_scaled = estimates[[i]], group = groups[[i]]))
            frame = cbind(data.frame(Sample = rownames(frame), frame))
            df_all <- rbind(df_all, frame)
            
            for (j in 1:ncol(comparisons)) { ## get test result for each comparison
                ## raw p
                test_p = sigTest(estimates[[i]][groups[[i]] %in% comparisons[1, j]], estimates[[i]][groups[[i]] %in%comparisons[2, j]])$p.value
                
                comp = data.frame(cell_type=cell_type, group=comparisons[1, j], group_comp=comparisons[2, j], raw_p= test_p)
                comp_p <- rbind(comp_p, comp)
            } ## end for each comparison in the cell type
        } ## end loop for each cell type
    }
    
    ## adj p for all cell types (same method as in Ogan's) for all comparisons
    comp_p$adj_p_all <- p.adjust(comp_p$raw_p, method='BH')
    
    ## combine the test results
    ## adj p for all cell types (same method as in Ogan's) for WT comparisons
    (index <- union(which(comp_p$group == 'WT'),which(comp_p$group_comp == 'WT'))) ## grep the Wt comparison for p adj
    comp_p$adj_p_wt <- NA
    comp_p$adj_p_wt[index] <- p.adjust(comp_p$raw_p[index], method='BH')
    
    ## join the pvalue
    df_tpm <- filterContain(comp_p, column='group_comp', value='WT')
    df_all <- noWarnings(left_join(df_all, df_tpm[, c('cell_type', 'group','adj_p_all', 'adj_p_wt')])) 
    
    
    
    returnlist=list(df_all, comp_p)
    return(returnlist)
    
}

#---------------------------------------------------------------------------#
# PART 2. main function wrapper for cell population estimate
#---------------------------------------------------------------------------#

if(is.null(genes)){
    ## do some quick fix of the markers 
    newcholin = c("Ecel1","Slc18a3","Crabp2","Zic3","Trpc3","Slc5a7","Lhx8") ## 2016-11-01
    mouseMarkerGenes$Striatum$ForebrainCholin = newcholin
    mouseMarkerGenes$Striatum$Dopaminergic  = mouseMarkerGenes$Midbrain$Dopaminergic    ## some dopaminergic expressed
    
    # do not remove markers from HD 2017-02-01        
    #         ## "Slc1a2","Slc1a3" are down regulated in HD striatum
    #         mouseMarkerGenes$Striatum$Astrocyte= setdiff(mouseMarkerGenes$Striatum$Astrocyte, c("Slc1a2","Slc1a3"))
    #         print(paste0('removed "Slc1a2","Slc1a3" from markers'))
    
    ## define marker genes
    if(disease =='AD'){
        genes = mouseMarkerGenes$Hippocampus
        print('Input Hippocampus markers')
    }else{
        genes = mouseMarkerGenes$Striatum
        print('Input striatum markers')
        
    }
}


if(!mixed_model){
    # ## load and preprocess data
    tmp <- preProcessMat(r_ob, probename = T, threshold = threshold)
    array_dat <- tmp[[1]]
    array_design <- tmp[[2]]
    dataset <- tmp[[3]]
    keyword  <- dataset
    (threshold <- tmp[[4]])
    
    ## make outdir
    timepoint <- levels(array_design$Timepoint)
    (outDir <- paste0(outDir, "/", dataset,"_", timepoint, "/"))
    ## genotypes for each sample
    if(design_group %in% colnames(array_design)){
        (groups <- as.character(array_design[, design_group]))
        
        print(paste0('groups are ', paste0(groups, collapse=', ')))
    }else{
        print(paste0(design_group, ' is not in array design'))
        return(NULL)
    }
    
}else{
    #############
    ## for mixed model rdata input
    #############
    print('input is mixed model')
    print(r_ob)
    tmp <- preProcessMat(r_ob, probename = F, threshold = threshold)
    array_dat <- tmp[[1]]
    array_design <- tmp[[2]]
    (threshold <- tmp[[4]])
    
    keyword  <- 'mixed_model'
    
    ## genotypes for each sample
    (groups <- as.character(array_design[, 'Disease_stage']))
    array_dat <- na.omit(array_dat)
}



print(paste0("outdir is ", outDir))
dir.create(outDir, showWarnings = F, recursive = T)


#' @param estimateFile is the matrix of estimate of each cell pop in each sample (dimentionaless)
#' @param comparisons if multiple genotypes are in the dataset, the comparisons will be a matrix with two rows (of the 2 genetype to compare
#' and columns are the number of comparison to be made
#' e.g. (comp <- rbind(c("wt", "wt"), c("mut1", "mut2"))) # not tested yet
#' @param groups (str list), e.g. genotypes corresponding to the exp matrix
#' @param exprData (df) exp matrix with probe, gene info (as char), and expression (as double/numeric). If the column is not
#'          an expression of the sample, it cannot be numeric (double)
#' @param geneColNames (str) the column name in the exprData that are gene symbols (matched with genes)
#' @param genes the markers genes in a list, each list correponding to a cell type, e.g.
#'      genes <- './brainCellTypeSpecificGenes-master/analysis/01.Gene Selection/FinalGenes/GabaDeep/Hippocampus/'
#'      genes <- puristOut(genes)
#' @param seekConsensus (bool), whether to remove marker genes with PC1 the are negatively correlations
#'      if T, remove genes with PC1<0
#' @param outDir is dir for plots and rotation table
#' @param geneTransform if using mouse datasets, = NULL, if using human, then need to ask Ogan...
#' 


## make the plot....    (Ogan's violin plots) only compare against WT, 
## the p value adj is for all cell types (for all comparisons)

fullEstimate(exprData = array_dat,
             genes=genes,
             geneColName="GeneSymbols",
             groups = groups,
             outDir=outDir,
             groupRotations=T,
             geneTransform = NULL,
             comparisons = 'all',
             estimateFile = paste0(outDir,'estimatefile.tsv'), ## the estimatefile is needed for plots
             pAdjMethod = "BH",
             outlierSampleRemove=F,
             seekConsensus=F  # remove probes when PC1 is negative (neg corr)
             # removeNegatives = removeNegatives
)


## get the estimates etc
x <-mgpEstimate(exprData = array_dat,
                genes=genes,
                geneColName="GeneSymbols",
                groups = groups,
                geneTransform = NULL,
                outlierSampleRemove=F,
                seekConsensus=F,  # remove probes when PC1 is negative (neg corr)
                # removeNegatives = removeNegatives,
                tableOut = paste0(outDir, "/", names(genes),
                                  " rotTable.tsv"), indivGenePlot = paste0(outDir,
                                                                           "/", names(genes), " indivExp", ".png"))
estimates <- x[[1]]
groups=x[[2]]
save(estimates, x, groups,sigTest,wt_only, file = paste0(outDir,'/', keyword, 'estimates.Rdata' ))

y <- cellPopAmountTest(estimates, groups, sigTest,wt_only)

df_all <- y[[1]]
comp_p <- y[[2]]


df_all$dataset = keyword
comp_p$dataset = keyword

writeTable(df_all, msg=paste0('# ', Sys.Date()), f_out= paste0(outDir,'/', keyword, '_cell_proportion_estimation_scaled.tsv' ))
writeTable(comp_p, msg=paste0('# ', Sys.Date()), f_out= paste0(outDir,'/', keyword, '_cell_proportion_test_pvalue.tsv' ))
print(paste0(outDir,'/', keyword, '_cell_proportion_estimation_scaled.tsv' ))



