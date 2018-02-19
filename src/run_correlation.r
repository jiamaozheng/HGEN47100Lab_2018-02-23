# Correlation Analysis 
#
# Author: Jiamao Zheng (jiamaoz@yahoo.com) 02/18/2018
#
# Rscript run_correlation.r
#
##

# install.packages('data.table')
# install.packages('dplyr')
# install.packages('tidyr)
# install.packages('ggplot2')
# install.packages('ggpmisc)
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggpmisc)))
suppressWarnings(suppressMessages(library(ggplot2)))

# setwd('~/Dropbox/Im_lab/2018/lab_session/HGEN47100Lab_2018-02-23/src/')

# arg variables
argv <- commandArgs(trailingOnly = TRUE)
gene_numbers <- argv[1]
if (length(argv) >1) {
  print('warning: only one numeric argument (e.g. 18 ) is allowed, or no argument by default!')
  break
}

# predicted gene expression 
print('load predicted gene expression..')
predictedExpression <- fread('../results/predixcan_results/geuvadis_predicted_expression.txt') %>% 
  gather(TargetID, predicted_expression, -FID, -IID) %>%
  select(TargetID, FID, predicted_expression)

# observed gene expression 
print('load observed gene expression..')
observedExpression <- fread('zcat < ../datasets/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz') %>% 
  select(-Gene_Symbol, -Chr, -Coord) %>%
  gather(FID, observed_expression, -TargetID)

# join between predicted and observed gene expression 
print('join predicted and observed gene expression..')
joined <- predictedExpression %>%
  inner_join(observedExpression, by=c('TargetID', 'FID'))

gene_anno <- fread('zcat < ../datasets/gencode.v12.annotation.gtf.gz') %>% # get ensemble_id and gene symbol pair 
  filter (V3=='gene') %>% filter(V1 != 'chrX' & V1 != 'chrY' & V1 != 'chrM') %>%
  select(V2, V3, V9) %>% 
  rowwise() %>% 
  mutate(V2=strsplit(gsub("\"", "", strsplit(V9, ';')[[1]][5]), ' ')[[1]][3]) %>% 
  mutate(V3=strsplit(gsub("\"", "", strsplit(V9, ';')[[1]][1]), ' ')[[1]][2]) %>% 
  select(V3, V2) %>% 
  distinct(V2, .keep_all = TRUE)
colnames(gene_anno) = c('TargetID', 'gene_symbol')
joined <- joined %>% inner_join(gene_anno, by=c('TargetID')) %>% select(gene_symbol, FID, predicted_expression, observed_expression)  # replace ensemble_id with gene_symbol 

# correlation plots for all genes 
if (length(argv) == 0) {
    genes <- unique(joined$gene_symbol)
} else if (length(argv) == 1) {
    if (as.numeric(gene_numbers) > length(unique(joined$gene_symbol))) {
       print(paste('warning: please enter a number less than ', length(unique(joined$gene_symbol)), sep=''))
       break 
    } else {
       genes <- unique(joined$gene_symbol)[1:gene_numbers]
    }
} 
summary_table <- data.frame()

for (gene in genes){
  print(paste('plot ../results/correlations/', gene, '.png',sep=''))
  filtered <- joined %>% filter(gene_symbol==gene)
  my.formula <- filtered$predicted_expression ~ filtered$observed_expression
  p <- ggplot(data=filtered, aes(x=predicted_expression, y = observed_expression)) + 
       geom_point() +
       labs(title= gene, x="Predicted expression", y=paste("Observed expression")) +
       geom_smooth(method='lm', se=FALSE, color="red") +
       stat_poly_eq(formula = my.formula, aes(label = paste( ..rr.label.., sep = "~~~")), parse = TRUE, size = 8) + 
       theme(axis.title.x = element_text(size=12, face='bold'), 
             axis.title.y=element_text(size=12, face='bold'), 
             title = element_text(size=12, face='bold')
    )
  # output plots
  ggsave(paste('../results/correlations/', gene, '.png',sep=''), width=8, height=8)

  # r2 and pvalue 
  cortest <- cor.test(filtered$observed_expression, filtered$predicted_expression, method = "kendall", alternative = "greater")
  r2=summary(lm(my.formula, filtered))$r.squared
  summary_table = rbind(summary_table, t(data.frame(c(gene, r2, cortest$p.value))))
}

# output summary table
print('output all genes r2 and pvalue to ../results/correlations/summary_correlation.csv')
colnames(summary_table) = c('gene_symbol', 'r2', 'pvalue')
write.csv(summary_table, file='../results/correlations/summary_correlation.csv', row.names=F)

print('Done!')
