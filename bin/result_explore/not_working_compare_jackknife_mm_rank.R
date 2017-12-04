f1 <- '/home/bzhuang/ND_project_combined/PD_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/early_down_regulation_mixed_model.tsv'
f2 <- '/home/bzhuang/ND_project_combined/PD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm/jackknife_early_down_regulation_mixed_model_results.tsv'
df1 <- read.delim(f1, comment.char = '#')
df2  <- read.delim(f2, comment.char = '#')

(my_order_list <- as.character(df1$geneSymbol))
(paper_list <- as.character(df2$geneSymbol)[1:1000]) ## top

#my_order_list <- as.character(df1$geneSymbol)
#paper_list <- as.character(df2$geneSymbol) ## top




plotAUC <- function(my_order_list, paper_list, f_save, dataset="", prefix ="", df_return =F){
    (len<-length(paper_list))
    tp <- vector()    # true positive
    fp <- vector()
    (tn=length(my_order_list)-len)  # true negative (any probes/genes not in pub paper)
    for(i in 1:length(my_order_list)){
        (tp[i]=length(intersect(my_order_list[1:i],paper_list))/len)   # if the gene list in my result is in the pub list, and rate of overlap
        (fp[i]=(i-tp[i])/tn)
    }
    
    (sum <- sum(tp))
    # calculate AUC
    auc_sum <- sum/length(my_order_list)
    cat("AUC=",auc_sum,"\n")
    
    # when all the genes are recalled:
    if(max(tp) == 1){
        index <- min(which(tp ==1))
        fp_threshold <- fp[index]
    }else {
        fp_threshold = 1
    }
    
    
    graphics.off()
    png(filename=paste0(f_save, prefix,"AUC.png") , width = 600, height = 600)
    plot(fp,tp,pch=".", main = paste0(" AUC: ", round(auc_sum, 4)), 
         sub = paste0("fp = ", round(fp_threshold, 4), " when tp = 1"))
    abline(v= fp_threshold, col = "gray")
    dev.off()
    
    if (df_return){
        returnlist <- list(auc_sum, tp, fp, fp_threshold)
        print(paste0("returnlist: auc, tp, fp, fp value when recall =1"))
        return(returnlist)
    }
}
