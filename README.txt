HGEN 47100 (Winter 2018) - Lab
-------------------------------

Introduction 
We will learn Matrix-eQTL and PrediXcan analysis during this lab. For more detailed information, please refer to the following manuscripts: Matrix-eQTL - https://academic.oup.com/bioinformatics/article/28/10/1353/213326; PrediXcan - https://www.nature.com/articles/ng.3367. 

Prerequisites 
1. Python2.7 (https://www.python.org/downloads/)
2. R 3.0+ (https://www.r-project.org)
3. data.table (https://github.com/Rdatatable/data.table): install.packages('data.table')
4. dplyr (https://github.com/tidyverse/dplyr): install.packages('dplyr')
5. ggpmisc (https://github.com/cran/ggpmisc): install.packages('ggpmisc)
6. ggplot2 (http://ggplot2.tidyverse.org): install.packages('ggplot2')
7. tidyr (http://tidyr.tidyverse.org): install.packages('tidyr)

Datasets/inputs/results/softwares/scripts 
Download a compressed file from S3("wget https://s3.amazonaws.com/imlab-open/HGEN47100/HGEN47100Lab_2018-02-23.zip"), unzip it("unzip HGEN47100Lab_2018-02-23.zip"), and navigate to the working "src" directory ("cd HGEN47100Lab_2018-02-23/src").
Please use tree command to display fold structure, which should be as below. 

HGEN47100Lab_2018-02-23/
├── README.txt
├── datasets
│   ├── GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
│   ├── README.txt
│   ├── gencode.v12.annotation.gtf.gz
│   ├── geuvadis.annot.txt.gz
│   ├── geuvadis.snps.txt.gz
│   └── gtex_v7_Cells_EBV-transformed_lymphocytes_imputed_europeans_tw_0.5_signif.db
├── inputs
│   ├── for_matrix-eqtl
│   ├── for_predixcan
│   └── inputs.zip
├── results
│   ├── correlations
│   ├── eqtls
│   ├── predixcan_results
│   └── results.zip
├── softwares
│   └── PrediXcan
└── src
    ├── prepare_inputs_for_matrix-eqtl.r
    ├── prepare_inputs_for_predixcan.r
    ├── run_correlation.r
    ├── run_matrix-eqtl.r
    └── run_predixcan.sh

Dataset Description (Adapted from https://github.com/hakyimlab/PredictDBPipeline/wiki/Tutorial) 
1) GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz: contains normalized expression quantification data in Lymphoblastoid Cell Lines from the GEUVADIS consortium. You can download it from here at https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/. Data normalization are detailed here at https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GeuvadisRNASeqAnalysisFiles_README.txt. 
2) gencode.v12.annotation.gtf.gz: a gene annotation file from gencode v12 which was used in the quantifying gene expression. You can download it from here at https://www.gencodegenes.org/releases/12.html. 
3) geuvadis.snps.txt.gz: genotype data of the samples. Only snps from HapMap are included. 
4) geuvadis.annot.txt.gz: the SNP annotation for the snps in genotype file. 
5) gtex_v7_Cells_EBV-transformed_lymphocytes_imputed_europeans_tw_0.5_signif.db: a prediction model used for PrediXcan analysis. You can download it from here at https://predictdb.org. 

Running
Step 1. Generate input files for Matrix-eQTL analysis. 
Given it will take around 10 hours to complete a job for all gene-snp pairs, we suggest that you look at only a few genes during this lab.

    For partial analysis (a few genes), please run a command such as "Rscript prepare_inputs_for_matrix-eqtl.r 8" in terminal. 

    For full analysis (all genes): please run "Rscript prepare_inputs_for_matrix-eqtl.r" command in terminal. 

The results will be automatically saved in the files ("geuvadis_expression_data.txt", "geuvadis_genotype_data.txt", "geuvadis_mRNA_location_data.txt", "geuvadis_snp_location_data.txt") under the directory of "HGEN47100Lab_2018-02-23/inputs/for_matrix-eqtl". 

Step 2. Run Matrix-eQTL analysis (Modified from http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html#cis).

	For default analysis (pvalue cutoff == 1e-6), please run "Rscript run_matrix-eqtl.r" command in terminal. 

	If you want to try a different pvalue cutoff, please run above command with an argument such as "Rscript run_matrix-eqtl.r 1e-10" in terminal. 

The results will be automatically saved in the files ("geuvadis_cis-eQTL.txt", "geuvadis_trans-eQTL.txt", "gevadis_eQTL_QQ-plot.png") under the directory of "HGEN47100Lab_2018-02-23/results/eqtls".

Step 3. Generate input files for PrediXcan analysis. 

	Please run "Rscript prepare_inputs_for_predixcan.r" command in terminal 

The results will be automatically saved in the files ("chr*.dosage.txt.gz", "samples.txt") under the directory of "HGEN47100Lab_2018-02-23/inputs/for_predixcan". 

Step 4. Run PrediXcan analysis (https://github.com/hakyimlab/PrediXcan)

	Please run "sh run_predixcan.sh" command in terminal 

The results will be automatically saved in the file "geuvadis_predicted_expression.txt" under the directory of "HGEN47100Lab_2018-02-23/results/predixcan_results".

Step 5. Run correlation analysis for predicted and observed gene expression. 

    For partial analysis (a few genes), please run a command such as "Rscript run_correlation.r 8" in terminal. 

    For full analysis (all genes): please run "Rscript run_correlation.r" command in terminal. 

The results will be automatically saved in the files ("*.png", "summary_correlation.csv") under the directory of "HGEN47100Lab_2018-02-23/results/correlations". 









