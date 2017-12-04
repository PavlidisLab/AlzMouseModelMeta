# 2016-08-20
### this is to make the early PD with 'GSE60413|GSE60414', change the pink1 to NA 
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
rm(list=setdiff(ls(),'home_dir'))
f= paste0(home_dir, '/PD_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm/early/mixed_model_results.Rdata')
# f= '/home//bzhuang/ND_project_combined//PD_mouse_model_project/mixed_model/random_intercept_include_NA/early/mixed_model_results.Rdata'
# save(array_dat,array_design, df_all, file=f)

load(f)
source('mixed_models/mixed_model.R')
(index <- which(row.names(array_dat) == 'Pink1'))

## rm value for all the pink KO
grep('GSE60413|GSE60414', colnames(array_dat), value=T)
array_dat[index, grep('GSE60413|GSE60414', colnames(array_dat))] <- NA



(i <- which(row.names(array_dat) == 'Pink1'))

model <- c('random_intercept')
full_report = T
to_plot =F
REML = F

df_t <- mixedModelStudy(array_dat, array_design, i, model=model, 
                        full_report = full_report,
                        to_plot= F,
                        REML=REML)

## rm original pink result
i <- which(row.names(df_all) != 'Pink1')
df_all<- df_all[i, ]%>%droplevels()
df_all <- rbindAllColumns(df_all, df_t)

save(array_dat,array_design, df_all, msg,df_all_expression, file=f)



# 2016-08-20
### this is to make the early PD with 'GSE60413|GSE60414', change the pink1 to NA 
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
rm(list=setdiff(ls(),'home_dir'))
f= paste0(home_dir, '/PD_mouse_model_project/mixed_model/random_intercept_include_NA/early/mixed_model_results.Rdata')


load(f)
source('mixed_models/mixed_model.R')
(index <- which(row.names(array_dat) == 'Pink1'))

## rm value for all the pink KO
grep('GSE60413|GSE60414', colnames(array_dat), value=T)
array_dat[index, grep('GSE60413|GSE60414', colnames(array_dat))] <- NA



(i <- which(row.names(array_dat) == 'Pink1'))

model <- c('random_intercept')
full_report = T
to_plot =F
REML = F

df_t <- mixedModelStudy(array_dat, array_design, i, model=model, 
                        full_report = full_report,
                        to_plot= F,
                        REML=REML)

## rm original pink result
i <- which(row.names(df_all) != 'Pink1')
df_all<- df_all[i, ]%>%droplevels()
df_all <- rbindAllColumns(df_all, df_t)

save(array_dat,array_design, df_all, file=f)


# 2016-08-20
### rm move pink results in late PD(3/8 of datasets with pink1 ko)
setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
rm(list=setdiff(ls(),'home_dir'))
f= paste0(home_dir, '/PD_mouse_model_project/mixed_model/random_intercept_include_NA/late/mixed_model_results.Rdata')


load(f)

## rm original pink result
i <- which(row.names(df_all) != 'Pink1')
df_all<- df_all[i, ]%>%droplevels()

save(array_dat,array_design, df_all, file=f)


setwd(paste0(home_dir, "/git/bin/mouse_dataset_process/"))
rm(list=setdiff(ls(),'home_dir'))
f= paste0(home_dir, '/PD_mouse_model_project/mixed_model/random_intercept_include_NA_low_exp_rm//late/mixed_model_results.Rdata')


load(f)

## rm original pink result
i <- which(row.names(df_all) != 'Pink1')
df_all<- df_all[i, ]%>%droplevels()

save(array_dat,array_design, df_all, msg,df_all_expression, file=f)
