#' 2016-09-25
#' make all mm jack results into a r data, and save to 
#' 'ND_results/mm_results/mm_jack_each_disease.Rdata'
#' use wrapper_plots_and_tables.R to run this script


#' examples
# source('config_wrappers.R')
# 
# regulation_ls = c('up', 'down')
# 
# phase_ls =c('early', 'late')
# variable_prefix =''
# input_dir = 'mixed_model_jackknife/random_intercept_include_NA_low_exp_rm/'  # where to grab the mixed_model_results.tsv
# output_r ='mm_jack_each_disease.Rdata'  ## output r data name
# source('summary_tables/summary_rdata/mm_jack_summary_each_disease.R')
# 
# 
# variable_prefix = '_adj_cell'
# input_dir = 'mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/'  # where to grab the mixed_model_results.tsv
# output_r ='mm_jack_each_disease_adj_cell.Rdata'  ## output r data name
# source('summary_tables/summary_rdata/mm_jack_summary_each_disease.R')




for(disease in disease_ls){
  # source('config/config_wrappers.R')
  (jack <- paste0(disease_dir,'/', input_dir))
  (jack_ls <- list.files(jack, pattern = 'mixed_model_results.tsv', full.names = T))
  
  for(regulation in regulation_ls){
    for(phase in phase_ls){
      (f <- grep(paste0(phase, "_", regulation), jack_ls, value = T))
      df <- read.delim(f, comment.char = '#')
      assign(paste0(disease, '_', phase,'_', regulation,"_mm_jack", variable_prefix), df)
    }
  }
}

out_dir <- file.path(home_dir, 'ND_results/mm_results/')
dir.create(out_dir, showWarnings = F)
(f_out <- paste0(out_dir, '/',output_r))
(cmd <- paste0('save(', 
               paste0(grep('mm_jack',ls(), value = T), collapse = ','), 
               ', file = "', f_out, '")'))
print(cmd)
eval(parse(text =cmd))
