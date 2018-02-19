# Matrix-eQTL Analysis (Adapted from http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis)
#
# Author: Jiamao Zheng (jiamaoz@yahoo.com) 02/18/2018
#
# Rscript run_matrix-eqtl.r
#
##

# install.packages('MatrixEQTL')
suppressWarnings(suppressMessages(library(MatrixEQTL)))

# setwd('~/Dropbox/Im_lab/2018/lab_session/HGEN47100Lab_2018-02-23/src/')

# arg variables
argv <- commandArgs(trailingOnly = TRUE)
pvalue <- as.numeric(argv[1])
if (length(argv) >1) {
  print('warning: the number of arguments should be less than 1')
  break
}

# choose model
print('choose modelLINEAR')
useModel = modelLINEAR

# set input file path 
print('set input file path')
SNP_file_name = '../inputs/for_matrix-eqtl/geuvadis_genotype_data.txt'
snps_location_file_name = '../inputs/for_matrix-eqtl/geuvadis_snp_location_data.txt'
expression_file_name = '../inputs/for_matrix-eqtl/geuvadis_expression_data.txt'
gene_location_file_name = '../inputs/for_matrix-eqtl/geuvadis_mRNA_location_data.txt'
# covariates_file_name = character() # no covariates 
# errorCovariance = numeric() # Error covariance matrix. Set to numeric() for identity.

# set output file path
print('set output file path')
output_file_name_cis = "../results/eqtls/geuvadis_cis-eQTL.txt"
output_file_name_tra = "../results/eqtls/geuvadis_trans_trans-eQTL.txt"

# only associations significant at this level will be saved
print('set pvalue cutoff')
if (length(argv) == 0) {
  pvOutputThreshold_cis = 1e-6;
  pvOutputThreshold_tra = 1e-6;
} else if (length(argv) == 1) {
  pvOutputThreshold_cis = pvalue;
  pvOutputThreshold_tra = pvalue;
}

# distance for local gene-SNP pairs
print('set distance for local gene-SNP pairs')
cisDist = 1e6;

# load genotype data 
print('load genotype data')
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile(SNP_file_name )

# load gene expression data 
print('load gene expression data')
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 20
gene$LoadFile(expression_file_name )

# # load covariates 
# cvrt = SlicedData$new()
# cvrt$fileDelimiter = "\t"
# cvrt$fileOmitCharacters = "NA"
# cvrt$fileSkipRows = 1
# cvrt$fileSkipColumns = 1
# cvrt$fileSliceSize = 10
# cvrt$LoadFile(covariates_file_name)

## Run the analysis
print('Run the analysis')
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  # cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel, 
  # errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE)

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')

## Plot the Q-Q plot of local and distant p-values
print('plot the Q-Q plot of local and distant p-values')
png('../results/eqtls/gevadis_eQTL_QQ-plot.png')
plot(me, pch = 16, cex = 0.7)
print('Done!')
dev.off()
