#2016-11-12


#  
# run the function in check_gene.R
cat(paste0('\n#############\n Must load all_ranks, all_ranks_adj\n'))





comparePublication <- function(df, threshold, phase, disease, gene_list_dir, prefix ='', out_dir,
                               write_out =T,
                               f_exp){

    
    ############
    ### get the top genes in my list
    ############
    (index <- which(df$disease == disease & df$phase == phase & df$up_jack <=threshold))
    df_up <- df[index, c('geneSymbol', 'pvalue_adj','up_jack')] %>%droplevels()
    df_up$regulation ='up'
    colnames(df_up)[3] <- 'Rank'
    
    (index <- which(df$disease == disease & df$phase == phase & df$down_jack <=threshold))
    df_down <- df[index, c('geneSymbol', 'pvalue_adj','down_jack')] %>%droplevels()
    df_down$regulation ='down'
    colnames(df_down)[3] <- 'Rank'
    
    df <- rbind(df_up, df_down)
    rownames(df) <- df$geneSymbol
    df <- orderCol(df, c('regulation', 'Rank'))
    
    ############
    #### get the paper gene list
    ############
    gene_dir <- paste0(gene_list_dir, '/', phase, '/')
    (f_ls <- list.files(gene_dir))
    (dataset_ls <- grep('GSE', unlist(strsplit(f_ls, '_')), value = T))
    f_ls <- paste0(gene_dir, f_ls)
    
    ## add NA genes if not included in the dataset    
    if (!is.null(f_exp)){
        df_exp <- read.delim(f_exp, comment.char = '#')
    }else{
        df_exp <- NULL
    }
    
    
    dd_sum_all <- NULL
    dd_all <- NULL
    for(dataset in dataset_ls){
        
        f <- grep(dataset, f_ls, value = T)
        print(f)
        dd <- read.delim(f, comment.char = '#')
        dd <- dd[, c('geneSymbol', 'regulation')]%>%droplevels()
        
        index <- which(dd$geneSymbol == '')
        if(length(index) >0){
            dd <- dd[-index, ]%>%droplevels()
        }
        dd <- rmDup(dd)
        
        ### add the NA genes of the top genes
        if(!is.null(df_exp)){
            
           
            
            if(dataset == "GSEKuhn"){
                (df_exp_NA <- which(apply(df_exp[, grep('GSE9857|GSE9375|GSE10202|GSE7958', colnames(df_exp))], 1, function(x) all(is.na(x))) == T))
                
            }else{
                (df_exp_NA <- which(is.na(df_exp[, grep(dataset, colnames(df_exp))[1]])))
            }
            
            
            if(length(df_exp_NA) >0){
                (na_genes <- setdiff(as.character(df_exp$geneSymbol)[df_exp_NA], dd$geneSymbol))
                df_exp_indi <- data.frame(geneSymbol = na_genes, regulation = 'no expression' )
                dd <- rbind(dd, df_exp_indi)
            }
        }
        
        ### get the summary
        (dd_sum <- data.frame(t(data.frame(summary(dd$regulation)))))
        tmp <- setdiff( c('up', 'down'),colnames(dd_sum))
        if(length(tmp) >0){
            tmp_f <- data.frame(x = 0)
            colnames(tmp_f) <- tmp
            dd_sum <- cbind(dd_sum, tmp_f)
        }
        
        rownames(dd_sum) <- dataset
        dd_sum$dataset <- dataset
        
        dd_sum <- dd_sum[, c('dataset', 'up','down')]
        if(is.null(dd_sum_all)){
            dd_sum_all <- dd_sum
        }else{
            dd_sum_all <- rbind(dd_sum_all,dd_sum)
        }
        
        ## change the col name
        colnames(dd)[2] <- dataset
        
        
        if(is.null(dd_all)){
            dd_all <- noWarnings(left_join(df, dd))
        }else{
            dd_all <- noWarnings(left_join(dd_all,dd))
        }
        
    }
    
    for(i in 5:ncol(dd_all)){
        dd_all[, i] <- as.character(dd_all[, i])
        dd_all[which(is.na(dd_all[, i])), i] <- ""
    }
    
    ## reorder regulations (up then down)
    dd_all <- dd_all[with(dd_all, order(regulation, decreasing = T)), ]
    
    ## add the summary of total studies show the gene
    tmpFun <- function(x){
        tmp <- table(as.character(x))
        index <- grep('up|down', names(tmp))
        if(length(index) == 1){
            return(as.numeric(tmp[index]))
            
        }else if (length(index) ==0){
            return('0')
        }else{
            return(paste0(paste(names(tmp)[index], tmp[index]), collapse=', '))
        }
    }
    
    dd_all$total_studies <- apply(dd_all[, grep('GSE', colnames(dd_all))], 1, function(x) tmpFun(x))
    
    
    ## add the total DE genes to the table
    dd_sum_all$total_genes <- dd_sum_all$up + dd_sum_all$down
    tmp <- data.frame(t(dd_sum_all[, c("total_genes")]))
    colnames(tmp) <- dd_sum_all$dataset
    colnames(dd_all)
    
    dd_all <- rbindAllColumns(dd_all,tmp)
    
    
    if(write_out){
        dir.create(out_dir, showWarnings = F, recursive = T)
        f_out <- paste0(out_dir,disease, '_', phase,'_', threshold,'_comparison_', prefix,Sys.Date(), '.tsv')
        writeTable(dd_all, f_out = f_out)
        
        f_out <- paste0(out_dir,disease, '_', phase,'_publication_summary_', Sys.Date(), '.tsv')
        writeTable(dd_sum_all, f_out = f_out)
        print(f_out)
    }
    
    return(list(dd_all, dd_sum_all))
}





###################
## RUN function
###################

