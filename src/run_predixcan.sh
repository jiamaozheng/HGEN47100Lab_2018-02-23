#!/bin/bash
#
# PrediXcan analysis (https://github.com/hakyimlab/PrediXcan) 
#
# sh run_predixcan.sh
#
##

# run predixcan
python ../softwares/PrediXcan/Software/PrediXcan.py \
--predict \
--dosages ../inputs/for_predixcan/ \
--dosages_prefix chr \
--samples samples.txt \
--weights ../datasets/gtex_v7_Cells_EBV-transformed_lymphocytes_imputed_europeans_tw_0.5_signif.db \
--output_prefix ../results/predixcan_results/geuvadis
