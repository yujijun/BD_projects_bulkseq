# --------------
# Date:  2019-12-23 12:49:30 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:This script is to generate a DESeqDataSeq file by DESeq2 package.
# 
library(DESeq2)
library(tidyverse)
####load the count table and sample info####
load("~/Documents/01-Work/06-BD_project/BD_projects/data/expr.count.rda")
expr.count <- expr.count[,-19]
expr.count <- expr.count[-grep("[.]",rownames(expr.count)),] #delete the gene with the name of "[.]"
sample.info <- data.frame(sample = colnames(expr.count)) #generate sample info by sample information.
sample.info <- sample.info %>% 
  mutate(patient =  str_remove(sample,pattern = '\\d+')) %>% 
  mutate(group =  case_when(patient=='HC' ~ 'HC',
                            TRUE ~ 'BD')) %>%
  column_to_rownames(var = 'sample')
sample.info$patient <- factor(sample.info$patient,levels = c("HC","BDU","BDV")) #the sample.info vector should be factor format.
sample.info$group <- factor(sample.info$group,levels = c("HC","BD"))

dds <- DESeqDataSetFromMatrix(countData =round(as.matrix(expr.count)), #create a DESeqDataSet
                                   colData = sample.info,
                                   design = ~ group) #if you would like to add more contrast, you can add them design = ~patient+group

#pre-filtering: perform a minimal pre-filetering to keep only rows that have at least 10 reads total.
dds <- dds[rowSums(counts(dds)) >= 10,]
dds<-estimateSizeFactors(dds) 
all(rownames(sample.info) == colnames(expr.count)) #check
dds_tmp <- DESeq(dds)








