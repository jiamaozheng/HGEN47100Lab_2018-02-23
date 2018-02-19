Dataset Description (Adapted from https://github.com/hakyimlab/PredictDBPipeline/wiki/Tutorial): 

1) GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz: contains normalized expression quantification data in Lymphoblastoid Cell Lines from the GEUVADIS consortium. You can download it from here at https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/. Data normalization are detailed here at https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GeuvadisRNASeqAnalysisFiles_README.txt. 

2) gencode.v12.annotation.gtf.gz: a gene annotation file from gencode v12 which was used in the quantifying gene expression. You can download it from here at https://www.gencodegenes.org/releases/12.html. 

3) geuvadis.snps.txt.gz: genotype data of the samples. Only snps from HapMap are included. 

4) geuvadis.annot.txt.gz: the SNP annotation for the snps in genotype file. 

5) gtex_v7_Cells_EBV-transformed_lymphocytes_imputed_europeans_tw_0.5_signif.db: a prediction model used for PrediXcan analysis. You can download it from here at https://predictdb.org. 