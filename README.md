# AD_mouse_model_project
AD_mouse_model_project

## setwd('bin/') to run all the scripts

## notes: need to install ermineJ to run GO enrichment analysis

## 1. install all required packages
`wrappers/install_packages.R`

## 2. update the home_dir in `config_wrappers.R`
all the downloaded files and results will be in subfolders in `home_dir`
e.g. `home_dir <- 'C:/Users/bzhuang/Documents/AD/'`

## 3. download all raw data and data QC
`wrappers/wrapper_download_data_QC_for_AD.R`

## 4. run to get all results (ranked gene lists)
`wrappers/wrapper_for_all_results.R`

## 5. to plot and make tables
`wrappers/wrapper_plots_and_tables.R`
