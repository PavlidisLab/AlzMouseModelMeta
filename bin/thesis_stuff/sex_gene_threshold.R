#' 2017-01-16
#' plot the sex gene expression of all samples from AD and HD early and late (quantile normalized but
#' not corrected)-- set the base of the threshold 6 fto remove lowly expressed genes (the mesian is 5.1 now)
#' 
#' 
#' ################################################################################################

#' ----------------------------------------------------- #
#' USAGE:

# datadir <-'/home/bzhuang/AD_mouse_model_project/data_and_QC/all_data/explore_data//analysis_mode/'  ## dir for limma toptables
# df_info <- '/home/bzhuang/AD_mouse_model_project/config_files/dataset_info.tsv'  #dataset label, phase label, order of the datasets etc.
# plot_dir <- "/home/bzhuang/AD_mouse_model_project/data_and_QC/all_data/explore_data/gender_expression_threshold/"
# 
# (r_obj_ls <- grep('.Rdata', list.files(datadir, recursive=T, full.names=T), value=T))
# method_t_ls <-c("max", "median", "75quantile", "95quantile") 
# 
# df <- loopExpThreshold(method_t_ls, r_obj_ls, plot_dir, width=1000, height=800, df_return =T)
#' ----------------------------------------------------- # 
#' Notes: 


library(dplyr)
library(tidyr)


df_gene <- read.table(header=TRUE, text='
                          GenderGenes    GeneSymbols
                          F    Xist
                          M    Kdm5d
                          M    Rps4y1')






#############################
## functions
#############################

probeExpGender <- function(r_obj, df_gene){
    #' r_obj dir for object (array_dat, annotation etc)
    #' @return df_exp_t: exp of non expressed gender genes
    #'          df_exp: exp of all gender genes
    #'          dataset,platform, array_dat, 
    ###-------------------------
    ## probe expression
    ###-------------------------
    load(r_obj)

    ## make sure the probe is in the array dat too
    (probe_ls <- intersect(as.character(df_gene$GeneSymbols), row.names(array_dat)))
    
    ## get the probes and expression
    (df_exp <- array_dat[probe_ls, ] %>% droplevels)
    df_exp$ProbeName <- row.names(df_exp)
    df_exp <- gather(df_exp, ProbeName, Expression)

    colnames(df_exp) <- c("ProbeName", "Sample", "Expression")
    
    ## combine with the array_design, and gene names
    df_exp <- noWarnings(left_join(df_exp, array_design[, c("Sample", "Gender")]))
    df_exp$GeneSymbols <- df_exp$ProbeName
    df_exp <- noWarnings(left_join(df_exp, df_gene[, c("GenderGenes", "GeneSymbols")]))
    df_exp$Sample <- as.factor(df_exp$Sample)

    
    ## sort df_exp by gender
    df_exp <- df_exp[order(df_exp$Gender), ]
    
    return(df_exp)

}





df_all <- NULL

for (disease in c('AD', 'HD')){
    print(disease)
    for(phase in c('early', 'late')){
        print(phase)
        r_obj <- paste0(home_dir,'/', disease, "_mouse_model_project/mixed_model/random_intercept_include_NA/", phase, "/expression.Rdata")
        df_exp <- probeExpGender(r_obj, df_gene)
        if(is.null(df_all)){
            df_all <- df_exp
        }else{ df_all <- rbind(df_all, df_exp)}
        
    }
}

df_exp <- df_all

###-------------------------
## get the gender of samples and threshold cut off for non expressed geneder gene
###-------------------------
## get the expression other gendergenes (i.e. probes that represent no expression)
## if the sample is mixed gender (multiple samples combined, then ignore these samples)
index <-NULL
for (i in 1: nrow(df_exp)){
    if (as.character(df_exp$Gender)[i] != 'mixed'){
        if(as.character(df_exp$Gender)[i] != as.character(df_exp$GenderGenes)[i]){
            index <- c(index, i)
        }
    }
    
}
df_exp_t <- df_exp[index, ]%>%droplevels()
(threshold <- median(df_exp_t$Expression, na.rm = T))
print(paste0('threshold is ', threshold))

df_exp <- na.omit(df_exp)
df_exp <- orderCol(df_exp, "Gender")
df_exp$Gender <- as.factor(df_exp$Gender)


## label the x labels by gender
df_tmp <-df_exp[, c("Sample", "Gender")]
df_tmp <- df_tmp[!duplicated(df_tmp), ] %>%droplevels()
df_tmp <- noWarnings(left_join(data.frame(Sample = unique(df_exp$Sample)), df_tmp))
df_tmp$Gender <- gsub("F", "female", df_tmp$Gender)
df_tmp$Gender <- gsub("M", "male", df_tmp$Gender)
label_col <- mapValue(df_tmp$Gender, c(male="blue", female="red", mixed ="brown"))
## must reorder the factor levels same as the new sample order
# (otherwise the label of gender colour will not match)
df_exp$Sample <- factor(df_exp$Sample, levels = unique(df_exp$Sample))  
df_exp$ProbeName <- relevel(as.factor(df_exp$ProbeName), ref = 'Xist')

## save plots

p<- plotXY(df_exp, y="Expression", x="Sample", group ="ProbeName",plt_alpha = 1, 
           title = paste0("Expression of gender genes: median of non-expressed: ", threshold))
p + geom_hline(yintercept = threshold, show.legend = T)+
    theme(axis.text.x = element_text(colour=as.character(label_col))) +
    scale_x_discrete(labels = df_exp$Gender)

ggsave(filename= p_out,width = 6, height = 5)
print(p_out)
