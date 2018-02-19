# To Prepare Inputs for PrediXcan Analysis (Software can be downloaded from https://github.com/hakyimlab/PrediXcan/tree/master/Software) 
#
# Author: Jiamao Zheng (jiamaoz@yahoo.com) 02/18/2018
#
# Rscript prepare_inputs_for_predixcan.r
#
##

# install.packages('data.table')
# install.packages('dplyr')
# install.packages('tidyr')
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))

# setwd('~/Dropbox/Im_lab/2018/lab_session/HGEN47100Lab_2018-02-23/src/')

# load snp 
print('load genotype file')
snps <- fread("zcat < ../datasets/geuvadis.snps.txt.gz") 

# prepare sample.txt 
print('prepare sample.txt file')
samples <- colnames(snps)[2:374]
sample_file <- data.frame(samples, samples)
write.table(sample_file, file = paste('../inputs/for_predixcan/', 'samples.txt', sep = ''), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# prepare genotype files (split 1 to 22)
print('split the whole genome genotype file into the 22 individual chromosome genotype files')
snp_anno <- fread('zcat < ../datasets/geuvadis.annot.txt.gz') %>% select(Chr, RSID_dbSNP137, Pos, Ref_b37, Alt, VariantID)
colnames(snp_anno) <- c('chr', 'rsid', 'pos', 'refAllele', 'effectAllele', 'Id')

final = snp_anno %>% inner_join(snps, by='Id')
for (i in 1:22) {
  print(paste('processing chr', i,  sep=''))
  
  final_temp <- final %>% filter(chr == i) %>% mutate(chr=paste('chr', i, sep=''))
  write.table(final_temp, file = gzfile(paste('../inputs/for_predixcan/chr', i, '.dosage.txt.gz', sep='')), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

print('done!')


