# AD_mouse_model_project

Scripts, data and supplementary information for Zhuang et al. "Mega-analysis of gene expression in mouse models of Alzheimer’s disease".

## setwd('bin/') to run all the scripts

## notes: You need to install ermineJ to run GO enrichment analysis

## 1. install all required packages
`wrappers/install_packages.R`

## 2. download all raw data and data QC
`wrappers/wrapper_download_data_QC_for_AD.R`

## 3. run to get all results (ranked gene lists)
`wrappers/wrapper_for_all_results.R`

## 4. to plot and make tables
`wrappers/wrapper_plots_and_tables.R`
