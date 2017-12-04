#' 2016-03-07
#' summary of the meta and jack results, what are the overlaps

# use wrapper.R, and 'low_exp_rm_analysis/wrapper_filter_low.R'
source("helper_functions.R")

# 1. summary table for jackknife, how many FDR <0.05, FDR<0.1

##---------------------------------------------
## functions
##---------------------------------------------
## all the counts are by genes
summaryTopGenes <- function(i, keyword, files, files_label="", df_return =F){
    # keyword: meta or jackknife

    print(files[i])
    df <- read.delim(files[i], comment.char="#")
    (filename <- files_label[i])
    
    if(keyword == "jackknife"){key_col <- "adj_combined_max_p"
    }else if(keyword == "meta"){
        key_col <- "Fisher_BH_adjusted"
    }
    
    (index <- which(df[, key_col] <0.05))
    (FDR5 <- length(index))
    FDR5genes_ls <- as.character(unique(df$geneSymbol[index]))
    FDR5genes <- paste0(FDR5genes_ls, collapse= "; ")
    (index <- which(df[, key_col] <0.1))
    (FDR10 <- length(which(df[, key_col] <0.1)))
    FDR10genes_ls <- as.character(unique(df$geneSymbol[index]))
    FDR10genes <- paste0(FDR10genes_ls, collapse= "; ")
    
    (index <- which(df[, key_col] <0.2))
    (FDR20 <- length(which(df[, key_col] <0.2)))
    FDR20genes_ls <- as.character(na.omit(unique(df$geneSymbol[index])[1:10]))
    FDR20genes <- paste0(FDR20genes_ls, collapse= "; ")
    
    (index <- which(df[, key_col] <0.3))
    (FDR30 <- length(which(df[, key_col] <0.3)))
    FDR30genes_ls <- as.character(na.omit(unique(df$geneSymbol[index])[1:10]))
    FDR30genes <- paste0(FDR30genes_ls, collapse= "; ")
    
    (total_probes <- nrow(df))
    df_summary <- data.frame(Method = paste0(keyword," (" ,key_col, ")"),
                             file = filename, total_Genes = total_probes, 
                             fdr0.05 = FDR5,  
                             fdr0.1 = FDR10,
                             fdr0.2 = FDR20,
                             fdr0.3 = FDR30,
                             fdr0.05_genes = FDR5genes,
                             fdr0.1_genes = FDR10genes,
                             fdr0.2_top_genes = FDR20genes,
                             fdr0.3_top_genes = FDR30genes)
    
    if(df_return){
        returnlist <- list(df_summary, FDR5genes_ls, FDR10genes_ls)
        return(returnlist)
    }else{
        return(df_summary)}
}

## get the counts
loopSummaryTopGenes <- function(keyword, file_folder,fo_s, df_return=F, known_modifier= "", 
                                phase_ls=c("early", "late", "med")){
    # known_modifier: modifier genes given by Joerg
    (files_label <- grep(keyword, list.files(file_folder), value=T))
    ## get only the phases 
    (search_phase <- paste0(phase_ls, collapse="|"))
    (files_label<- grep(search_phase, files_label, value=T))
    (files <- paste0(file_folder, files_label))
    
    if(known_modifier != ""){
        df_known <- read.delim(known_modifier, comment.char="#")
        df_known$gene_modifier <- paste0(df_known$GeneSymbols, "(", df_known$Gene_name, ")")
    }
    
    df_all <- NULL
    for (i in 1:length(files)){
        x <- summaryTopGenes(i, keyword, files, files_label, df_return =T)
        df <-x[[1]]
        fdr_5<- x[[2]]
        fdr_10 <- x[[3]]
        if(known_modifier != ""){
            df$fdr0.05_gene_modifiers <- paste0(df_known$gene_modifier[which(df_known$GeneSymbols %in% fdr_5)], collapse="; ")
            df$fdr0.1_gene_modifiers <- paste0(df_known$gene_modifier[which(df_known$GeneSymbols %in% fdr_10)], collapse="; ")
        }
        df_all <- rbind(df_all, df)
        
    }

    ## save results 
    write.table(df_all, file = fo_s, row.names = F,
                sep ='\t', quote = F)
    
    print(fo_s)
    if(df_return){
        return(df_all)
    }
}

## get the venndiagram lists for FDR0.05 and 0.1
getSummaryList<- function(keyword, file_folder, phase_ls = c("early", "late", "med")){
    # phase_ls: the disease phase to be summerized
    # keyword is  "jackknife", "meta"
    # if a list is empty, do not return
    (files_label <- grep(keyword, list.files(file_folder), value=T))
    (files <- paste0(file_folder, files_label))

    cmd_5 <- NULL
    cmd_10 <- NULL
    phase_venn_ls <- NULL
    reg_venn_ls <- NULL
    regulation_ls <- c("up", "down")
    for(phase in phase_ls){
        for (reg in regulation_ls){
            (search_keyword <- paste0(phase, ".*", reg, "_regulation"))
            (file_l <- grep(search_keyword, files, value=T))
            if(length(file_l) == 1){
                df <- summaryTopGenes(i=1, keyword, file_l,df_return=T)
                (label <- paste0(phase, "_", reg))
                phase_venn_ls <- c(phase_venn_ls, phase)
                reg_venn_ls <- c(reg_venn_ls, reg)

                fdr_ls_5 <- df[[2]]
                fdr_ls_10 <- df[[3]]
                assign(paste0("fdr_ls_5_",label), eval(as.symbol("fdr_ls_5")))
                (cmd_5 <- c(cmd_5, paste0(label, "=", paste0("fdr_ls_5_",label))))
                assign(paste0("fdr_ls_10_",label), eval(as.symbol("fdr_ls_10")))
                (cmd_10 <- c(cmd_10, paste0(label, "=", paste0("fdr_ls_10_",label))))
            }else{
                print(paste0("Phase ", phase, " not found, skipped."))
            }
        }
    }
    
    ## get the venn list for plot
    (cmd <- paste0("list(", paste0(cmd_5, collapse=", "), ")") )
    assign("venn_list_5", eval(parse(text = cmd)))
    
    (cmd <- paste0("list(", paste0(cmd_10, collapse=", "), ")") )
    assign("venn_list_10", eval(parse(text = cmd)))
    returnlist <- list(venn_list_5, venn_list_10, phase_venn_ls, reg_venn_ls)
    return(returnlist)
}

## plot bar graph:
## plot bar graph for poster:
plotBarPoster <- function(df, key_col, plot_title, fo_s){
    ggplot(df, aes(x = phase, fill = regulation)) +
        geom_bar(stat="identity", aes_string(y=key_col), position="dodge") +
        geom_text(aes_string(x="phase", y=key_col,label=key_col), 
                  position = position_dodge(width=0.9), size=8) +  # asign count text to each bar
        ggtitle(plot_title) +
        theme_bw() +   ## white background
        ylab("Count") +
        xlab("") +
        theme(plot.title = element_text(size = 32), text =element_text(size = 20), # theme for poster size
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        ggsave(filename = fo_s, width = 10, height =6)
}
        

##---------------------------------------------
## function wrap(main function)
##---------------------------------------------


mainSummaryMetaJack <- function(file_folder, f_out_dir, prefix, phase_ls, known_modifier ="",
                                venn_plot =F){
    #' file folder: where the jackknife and meta files are 
    #' f_out_dir: result and plot output
    #' prefix: for plot names
    #' @return plots, tables and venn diagrams
       dir.create(f_out_dir, showWarnings=F, recursive=T)
       (search_phase <- paste0(phase_ls, collapse="|"))
       
       keyword <-"meta"
       (fo_s <-paste0(f_out_dir, "/summary_", prefix, keyword,"_", Sys.Date(), ".tsv"))
       df_all_meta <- loopSummaryTopGenes(keyword, file_folder,fo_s, df_return=T,
                                          known_modifier=known_modifier, phase_ls = phase_ls)
       
       keyword <- "jackknife" 
       (fo_s <-paste0(f_out_dir, "/summary_", prefix, keyword,"_", Sys.Date(), ".tsv"))
       df_all_jack <- loopSummaryTopGenes(keyword, file_folder,fo_s, df_return=T, 
                                          known_modifier=known_modifier, phase_ls = phase_ls)
       
       
       ## draw the bar graph
       df_meta <- df_all_meta[,c(2, 4, 5)]
       df_meta$method <- "Meta_signature_genes"
       (time_label <- grep(search_phase, unlist(strsplit(as.character(df_meta$file), split ="_")), value=T))
       (reg_label <- grep("up|down", unlist(strsplit(as.character(df_meta$file), split ="_")), value=T))
       
       df_meta$phase <- time_label
       df_meta$regulation <- reg_label
       df_meta$regulation <- factor(df_meta$regulation, levels = c("up", "down"))
       
       df_jack <- df_all_jack[,c(2, 4, 5)]
       df_jack$method <- "Jackknife_genes"
       (time_label <- grep(search_phase, unlist(strsplit(as.character(df_jack$file), split ="_")), value=T))
       (reg_label <- grep("up|down", unlist(strsplit(as.character(df_jack$file), split ="_")), value=T))
       
       df_jack$phase <- time_label
       df_jack$regulation <- reg_label
       df_jack$regulation <- factor(df_jack$regulation, levels = c("up", "down"))
       
       df_plot <- rbind(df_meta, df_jack)
       colnames(df_meta)[c(2,3)] <- c("FDR0.05", "FDR0.1")
       colnames(df_jack)[c(2,3)] <- c("FDR0.05", "FDR0.1")
       
       
       
       df<-df_meta
       key_col <- "FDR0.05"
       plot_title <-"Meta-analysis signature genes (FDR <0.05)"
       (fo_s <-paste0(f_out_dir, "/summary_", prefix, "meta_barplot_", Sys.Date(), ".png"))
       plotBarPoster(df, key_col, plot_title, fo_s)
       
       ##
       df<-df_jack
       key_col <- "FDR0.05"
       plot_title <-"Jackknife core genes (FDR <0.05)"
       (fo_s <-paste0(f_out_dir, "/summary_", prefix, "jack0.05_barplot_", Sys.Date(), ".png"))
       plotBarPoster(df, key_col, plot_title, fo_s)
       
       ##
       df<-df_jack
       key_col <- "FDR0.1"
       plot_title <-"Jackknife core genes (FDR <0.1)"
       (fo_s <-paste0(f_out_dir, "/summary_", prefix, "jack0.1_barplot_", Sys.Date(), ".png"))
       plotBarPoster(df, key_col, plot_title, fo_s)
       
       ##---------------------#
       ## find overlap genes and plot Venn (optional)
       ##----------------------#
       (files_dir <- paste0(f_out_dir, Sys.Date(), "_venn_plots/"))
       dir.create(files_dir, showWarnings=F)
       
       keyword_ls  <- c("jackknife", "meta")
       for(keyword in keyword_ls){
           vennlist <- getSummaryList(keyword, file_folder, phase_ls)
           vennlist_5 <- vennlist[[1]]
           vennlist_10 <- vennlist[[2]]
           (venn_combo <- vennlist[[3]]) ## phases are compared to other phases, but not the same phase
           (reg_combo <- vennlist[[4]]) ## for more than 2 phases
           
           combo_list_cmd <- NULL
           (x <- combn(1: (length(venn_combo)), m= 2))
           for (i in 1:dim(x)[2]){
               combo <- venn_combo[x[, i]]
               if(combo[1] != combo[2]){ # if the phases are different
                   new_cmd <- paste0("c(", paste0(x[, i], collapse=","), ")")
                   combo_list_cmd <- c(combo_list_cmd, new_cmd)
               }
           }
           ## if there are 3 phase groups: only compare the same direction ??
           if(length(reg_combo)/2 >2){
               for (i in c("up", "down")){
                   (new_cmd <- paste0("c(", paste0(grep(i, reg_combo), collapse =","), ")"))
                   combo_list_cmd <- c(combo_list_cmd, new_cmd)
               }
           }
           
           ## assign the venn combo
           combo_list_cmd <- paste0("list(", paste0(combo_list_cmd, collapse=", "), ")")
           assign("combo_ls", eval(parse(text = combo_list_cmd)))
           
           for(threshold in c(0.05, 0.1)){
               print(threshold)
               (msg <- paste0('# ', Sys.Date(),", ", keyword, " genes.", "\n# Threshold: ", threshold))
               if(threshold ==0.05){
                   vennlist_all <- vennlist_5
               }else{
                   vennlist_all <- vennlist_10
               }
               print("plotting venn list")
               for(i in 1:length(combo_ls)){
                   (combo <- combo_ls[[i]])
                   vennlist <- vennlist_all[combo]
                   
                   ## plot only when the list are not both empty
                   x = 0
                   for (j in 1:length(vennlist)){
                       x <- x + length(vennlist[[j]])
                   }
                   (title <- paste0(keyword,"_", paste0(names(vennlist), collapse="_"), "_FDR_", threshold))
                   if(all(x != 0, venn_plot)){ ## plot venn plot
                       graphics.off()
                       png(filename=paste0(files_dir, title, "_venn.png") , width = 600, height = 600)
                       vennDiagram2(vennlist, title = title, cat.cex =1.5, cex =1.5)
                       dev.off()
                   }
                   ##find overlap of all lists
                   (temp <- Reduce(intersect, vennlist))
                   msg <- paste0(msg, "\n#\n# ", title,
                                 "\n# TERMS OVERLAPPED for UP regulated terms: ", length(temp),"\n",
                                 paste0(temp, collapse="\n"))
               }
               f_out <-paste0(files_dir, keyword, "_overlap_threshold_", threshold, ".txt") 
               print(f_out)
               sink(f_out, type="output")
               writeLines(msg)
               sink()
           }
       }


}





