#!/bin/bash
# ErmineJ script for AD, 2017-02-02
# Input folder: /home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop/
# Background folder: /home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//ermineJ_background_all_processes/
# Result out dir: /home/bzhuang/ND_project_combined/AD_mouse_model_project//ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop_all_processes/2017-02-02_geneset_500/
                  
#############
                  
# settings 
                  
#############
                  
iteration=200000
ORAthreshold=0.1
                  
gsrmethod=PRECISIONRECALL   # set method for GSR
                  
maxsize=500
minsize=10	# set the min class size
                  
xml=/home/bzhuang/ermineJ.data/go_daily-termdb.rdf-xml.gz
testmethod=GSR  #choose one of GSR (need to use -m for method), ORA(need to use -t for threshold, default =0.001), ROC, CORR


############
# jackknife_early_down_regulation_mixed_model_results.tsv
############


background=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//ermineJ_background_all_processes//bg_jackknife_early_down_regulation_mixed_model_results.tsv
infile=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//jackknife_early_down_regulation_mixed_model_results.tsv
outfile=/home/bzhuang/ND_project_combined/AD_mouse_model_project//ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop_all_processes/2017-02-02_geneset_500/ermineJ_jackknife_early_down_regulation_mixed_model_results.tsv
scorecol=2


sh $ERMINEJ_HOME/bin/ermineJ.sh -a $background -c $xml -n $testmethod -o $outfile -s $infile -e $scorecol -g BEST -x $maxsize -y $minsize -j -i $iteration -m $gsrmethod -t $ORAthreshold -l


############
# jackknife_early_up_regulation_mixed_model_results.tsv
############


background=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//ermineJ_background_all_processes//bg_jackknife_early_up_regulation_mixed_model_results.tsv
infile=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//jackknife_early_up_regulation_mixed_model_results.tsv
outfile=/home/bzhuang/ND_project_combined/AD_mouse_model_project//ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop_all_processes/2017-02-02_geneset_500/ermineJ_jackknife_early_up_regulation_mixed_model_results.tsv
scorecol=2


sh $ERMINEJ_HOME/bin/ermineJ.sh -a $background -c $xml -n $testmethod -o $outfile -s $infile -e $scorecol -g BEST -x $maxsize -y $minsize -j -i $iteration -m $gsrmethod -t $ORAthreshold -l


############
# jackknife_late_down_regulation_mixed_model_results.tsv
############


background=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//ermineJ_background_all_processes//bg_jackknife_late_down_regulation_mixed_model_results.tsv
infile=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//jackknife_late_down_regulation_mixed_model_results.tsv
outfile=/home/bzhuang/ND_project_combined/AD_mouse_model_project//ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop_all_processes/2017-02-02_geneset_500/ermineJ_jackknife_late_down_regulation_mixed_model_results.tsv
scorecol=2


sh $ERMINEJ_HOME/bin/ermineJ.sh -a $background -c $xml -n $testmethod -o $outfile -s $infile -e $scorecol -g BEST -x $maxsize -y $minsize -j -i $iteration -m $gsrmethod -t $ORAthreshold -l


############
# jackknife_late_up_regulation_mixed_model_results.tsv
############


background=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//ermineJ_background_all_processes//bg_jackknife_late_up_regulation_mixed_model_results.tsv
infile=/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop//jackknife_late_up_regulation_mixed_model_results.tsv
outfile=/home/bzhuang/ND_project_combined/AD_mouse_model_project//ermineJ/mixed_model_jackknife/random_intercept_include_NA_low_exp_rm_adj_cell_pop_all_processes/2017-02-02_geneset_500/ermineJ_jackknife_late_up_regulation_mixed_model_results.tsv
scorecol=2


sh $ERMINEJ_HOME/bin/ermineJ.sh -a $background -c $xml -n $testmethod -o $outfile -s $infile -e $scorecol -g BEST -x $maxsize -y $minsize -j -i $iteration -m $gsrmethod -t $ORAthreshold -l
