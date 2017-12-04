#' 2016-03-04
#' compare the results from ermineJ to see what pathways are overlapped
rm(list=setdiff(ls(),'home_dir'))
library("HelperFunctions")
library("dplyr")
library("VennDiagram")

# ermineJ ## input the ermineJ results

files_dir <- "/home/bzhuang/AD_mouse_model_project/ermineJ//2016-03-07//results/"
(files <- grep("meta", list.files(files_dir, full.names=T, pattern = "tsv"), value =T))
(files_label <- grep("meta", list.files(files_dir, full.names=F, pattern = "tsv"), value =T))

# select the top gene sets
top_label <- 30

## get all the top gene set labels
df_all <- NULL
for (i in 1:length(files)){
    print(files[i])
    df <- read.delim(files[i], comment.char="#")
    df <- df[1:top_label, -c(1)]
    df$label <- files_label[i]
    df_all <- rbind(df_all, df)
    print(paste0(files_label[i], " done"))
    tmp_ls <- list(as.character(df$Name))
    names(tmp_ls) <- files_label[i]
    assign(x = paste0("ls_", i), tmp_ls) 
}

write.table(df_all, file = paste0(files_dir, "top_", top_label, ".tsv"), row.names = F,
            sep ='\t', quote = F)



#------------#
# venndiagram - overlapped
#------------#
venn_list <- list(early_down = ls_1[[1]], early_up = ls_2[[1]], late_down = ls_3[[1]], late_up = ls_4[[1]])
(title= paste0("Up_regulated"))
graphics.off()
png(filename=paste0(files_dir, title, "_venn.png") , width = 600, height = 600)
vennDiagram2(venn_list[c(2,4)], title = title, cat.cex =1.5, cex =1.5)
dev.off()
temp <- intersect(ls_2[[1]], ls_4[[1]])

msg <- paste0('# ', Sys.Date(), "Top ", top_label, "biological pathways ",
#               "\n# Top 20 terms for early phases UP regulated: \n",
#               paste0(ls_2[[1]], collapse="; "),
#               "\n# Top 20 terms for late phases UP regulated: \n",
#               paste0(ls_4[[1]], collapse="; "),
              "\n# TERMS OVERLAPPED for UP regulated terms: ", length(temp),"\n",
              paste0(temp, collapse="\n"))


cat(msg)



(title= paste0("Down_regulated"))
graphics.off()
png(filename=paste0(files_dir, title, "_venn.png") , width = 600, height = 600)
vennDiagram2(venn_list[c(1,3)], title = title, cat.cex =1.5, cex =1.5)
dev.off()
temp <- intersect(ls_1[[1]], ls_3[[1]])
msg <- paste0(msg, 
#               "\n# Top 20 terms for early phases DOWN regulated: \n",
#               paste0(ls_1[[1]], collapse="\n"),
#               "\n# Top 20 terms for late phases DOWN regulated: \n",
#               paste0(ls_3[[1]], collapse="\n"),
              "\n# TERMS OVERLAPPED for DOWN regulated terms: ", length(temp),"\n",
              paste0(temp, collapse="\n"))



(title= paste0("early_up_late_down_regulated"))
graphics.off()
png(filename=paste0(files_dir, title, "_venn.png") , width = 600, height = 600)
vennDiagram2(venn_list[c(2, 3)], title = title, cat.cex =1.5, cex =1.5)
dev.off()
temp <- intersect(ls_2[[1]], ls_3[[1]])

msg <- paste0(msg,
              "\n# TERMS OVERLAPPED for early UP regulated and late DOWN regulated terms: ", length(temp),"\n",
              paste0(temp, collapse="\n"))


(title= paste0("early_down_late_up_regulated"))
graphics.off()
png(filename=paste0(files_dir, title, "_venn.png") , width = 600, height = 600)
vennDiagram2(venn_list[c(1,4)], title = title, cat.cex =1.5, cex =1.5)
dev.off()
(temp <- intersect(ls_1[[1]], ls_4[[1]]))
msg <- paste0(msg,
              "\n# TERMS OVERLAPPED for early DOWN regulated and late UP regulated terms: ", length(temp),"\n",
              paste0(temp, collapse="\n"))

f_out <-paste0(files_dir, "top_", top_label, "_terms_and_overlap.txt") 
sink(f_out, type="output")
writeLines(msg)
sink()

