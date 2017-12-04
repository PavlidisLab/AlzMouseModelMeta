#' 2016-09-25
#' make all mm jack results into a r data, and save to 
#' 'ND_results/mm_results/mm_jack_each_disease.Rdata'




for(disease in disease_ls){
  # source('config/config_wrappers.R')
  (jack <- paste0(disease_dir, input_dir))
  (jack_ls <- list.files(jack, pattern = 'mixed_model_results.tsv', full.names = T))
  
  for(regulation in regulation_ls){
    for(phase in phase_ls){
      (f <- grep(paste0(phase, "_", regulation), jack_ls, value = T))
      df <- read.delim(f, comment.char = '#')
      assign(paste0(disease, '_', phase,'_', regulation,"_mm_jack", variable_prefix), df)
    }
  }
}

out_dir <- paste0(home_dir, 'ND_results/mm_results/')
dir.create(out_dir, showWarnings = F)
(f_out <- paste0(out_dir, output_r))
(cmd <- paste0('save(', 
               paste0(grep('mm_jack',ls(), value = T), collapse = ','), 
               ', file = "', f_out, '")'))
print(cmd)
eval(parse(text =cmd))
