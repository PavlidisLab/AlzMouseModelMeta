#2updated 2017-02-01 -- not in use anymore
## get the top genes and their associated enriched paths from the ermineJ analysis
# only for cell pop adj results
# for biological process and 
folder <- './11-25_results/'

dir.exists(folder)
(f_ls <- grep('adj_cell_pop.*analysis', list.dirs(folder), value = T))


(files <- list.files(f_ls, pattern = 'top_genes', full.names = T))

process_key_ls = c('','all_process')
disease_ls =c('AD', 'HD', 'combined')
phase_ls = c('early', 'late')
regulation_ls =c('up','down')

df_all <- NULL
for(process_key in process_key_ls){
    for(disease in disease_ls){
        for (phase in phase_ls){
            for(regulation in regulation_ls){
                
                (f = grep(paste0(disease, '.*', phase, '_', regulation), files, value = T))
                if(process_key == ''){
                    f=min(f)
                }else{(f=grep(process_key, f, value = T))}
                print(f)
                
                df <- read.delim(f, comment.char = '#')
                df <- cbind(data.frame(file =f), df)
                df$disease = disease
                df$phase = phase
                df$regulation = regulation
                df$process = process_key    
                if(is.null(df_all)){
                    df_all <- df
                }else{
                    df_all <- rbind(df_all, df)
                }
            }#reg
        }#phase
    }#disease
}#process


df <- na.omit(df_all)
df1 <- df[which(df$process == ''), ]
colnames(df1)[c(4,6)] <- c('count_bio','sig_paths_bio')
df2 <- df[which(df$process == 'all_process'), ]
colnames(df2)[c(4,6)] <- c('count_all','sig_paths_all')
library(dplyr)
sig_path_genes <- full_join(df1[, -c(1,5, 10)], df2[, -c(1,5,10)])
save(sig_path_genes, file = '../../results/ermineJ11-26_sig_genes.Rdata')
