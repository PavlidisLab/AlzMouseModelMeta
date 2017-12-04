load("C:/Users/User/Documents/lab_results/ND_project_combined/ND_results/mm_results/early_mm_results.Rdata")
load("C:/Users/User/Documents/lab_results/ND_project_combined/ND_results/mm_results/late_mm_results.Rdata")

library(HelperFunctions)
library(dplyr)
genes <- c('Trem2', 'Picalm','Cd33', 'Clu', 'Cr1', 'Abca7', 'Apoe',
           'Cd2ap','Epha1', 'Bin1') ## AD risk genes




df_all <- NULL

for (i in grep('mm_variable', ls(), value = T)){
  print(i)
  (cmd <- paste0("df <- filterContain(", i, "[1:10], column = 'geneSymbol', value = genes)"))
  eval(parse(text = cmd))
  df$data <- i
  if(is.null(df_all)){
    df_all <- df
  }else{
    df_all <- rbind(df_all, df)
  }
}
