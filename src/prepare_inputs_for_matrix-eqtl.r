# To Prepare Inputs for Matrix-eQTL Analysis (http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis)
#
# Author: Jiamao Zheng (jiamaoz@yahoo.com) 02/18/2018
#
# Rscript prepare_inputs_for_matrix-eqtl.r

# install.packages('data.table')
# install.packages('dplyr')
# install.packages('tidyr')
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))

# setwd('~/Dropbox/Im_lab/2018/lab_session/HGEN47100Lab_2018-02-23/src/')

# arg variables
argv <- commandArgs(trailingOnly = TRUE)
gene_number <- as.numeric(argv[1])
if (length(argv) >1) {
  print('warning: the number of arguments should be less than 1')
  break
}

# load data 
print('load and prefilter data..')
gene_anno <- fread('zcat < ../datasets/gencode.v12.annotation.gtf.gz') %>% filter (V3=='gene') %>% filter(V1 != 'chrX' & V1 != 'chrY' & V1 != 'chrM') # only keep chrosomal 1-22 
snp_anno <- fread('zcat < ../datasets/geuvadis.annot.txt.gz') %>% filter(length(Ref_b37) > 1 & length(Alt) > 1) # remove snps with the size larger than 1 
gene <- fread('zcat < ../datasets/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz')
snp <- fread('zcat < ../datasets/geuvadis.snps.txt.gz')

# snp location
print('prepare snp location file..')
snps_location <- snp_anno %>% select(RSID_dbSNP137, Chr, Pos)
colnames(snps_location) <- c('rsid', 'chr', 'pos')

# snp 
print('prepare snp file..')
snps <- snp %>% 
  inner_join(snp_anno %>% select(VariantID, RSID_dbSNP137), by = c('Id'='VariantID')) %>% 
  select(one_of(c('RSID_dbSNP137', intersect(colnames(snp),colnames(gene))))) # common cases 

# gene location
print('prepare gene location file..')
gene_location <- gene_anno %>%
  mutate(V1=gsub('chr', '', V1)) %>% 
  select(V1, V2, V3, V4, V5, V9) %>% 
  rowwise() %>% 
  mutate(V2=strsplit(gsub("\"", "", strsplit(V9, ';')[[1]][5]), ' ')[[1]][3]) %>%  # gene symbol
  mutate(V3=strsplit(gsub("\"", "", strsplit(V9, ';')[[1]][1]), ' ')[[1]][2]) %>%  # ensemble id 
  select(V3, V2, V1, V4, V5) %>% 
  distinct(V2, .keep_all = TRUE)  # remove gene isoforms 
colnames(gene_location) = c('ensemble_id', 'gene_symbol', 'chr', 's1', 's2')

# gene expression 
print('prepare gene expression file..')
gene_expression <- gene %>% 
  select(-Gene_Symbol, -Chr) %>% 
  inner_join(gene_location %>% select(ensemble_id, gene_symbol, s1), by = c('TargetID'='ensemble_id')) %>% 
  mutate(TargetID=gene_symbol) %>%
  select(-s1, -gene_symbol, -Coord) %>%
  select(one_of(c('TargetID', intersect(colnames(snp),colnames(gene))))) # common cases
if (length(argv) == 1) {
  gene_expression <- gene_expression[1:gene_number, ]
}

# remove column id for gene_location 
gene_location <- gene_location %>% select(-ensemble_id)

# write data 
print('write ../inputs/for_matrix-eqtl/geuvadis_genotype_data.txt')
write.table(snps, file = '../inputs/for_matrix-eqtl/geuvadis_genotype_data.txt', sep = '\t', row.names = F, quote = F)

print('write /inputs/for_matrix-eqtl/geuvadis_snp_location_data.txt')
write.table(snps_location, file = '../inputs/for_matrix-eqtl/geuvadis_snp_location_data.txt', sep = '\t', row.names = F, quote = F)

print('write ../inputs/for_matrix-eqtl/geuvadis_expression_data.txt')
write.table(gene_expression, file = '../inputs/for_matrix-eqtl/geuvadis_expression_data.txt', sep = '\t', row.names = F, quote = F)

print('write ../inputs/for_matrix-eqtl/geuvadis_mRNA_location_data.txt')
write.table(gene_location, file = '../inputs/for_matrix-eqtl/geuvadis_mRNA_location_data.txt', sep = '\t', row.names = F, quote = F)

print('Done!')