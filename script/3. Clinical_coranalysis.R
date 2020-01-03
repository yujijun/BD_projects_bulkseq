library(tidyverse)
require(gdata)
library(corrplot)
#preprocessing about clinicale dataset####
#This is the data preprocessing 
patient_info <- read.xls("/Users/yujijun/Documents/01-Work/06-BD_project/BD_projects/data-raw/clinical_info/patient info.xlsx",sheet = 2)
patient_info_important <- patient_info %>% filter(improtant == "Yes")
patient_info_all <- patient_info_important[,-grep("是否异常",colnames(patient_info_important))]
clinical_info <- patient_info_all[,-c(1,3,4,5)]
rownames(clinical_info) <- clinical_info$item
clinical_info <- clinical_info[,-1]
colnames(clinical_info) <- c(paste0("BDV",seq(1,6)),paste0("BDU",seq(1,3)))
#I also generate a clinical_info.rda file under "~/Documents/01-Work/06-BD_project/BD_projects/data/clinical_info.rda"
#so next time, user can load clinical dataset by
load("~/Documents/01-Work/06-BD_project/BD_projects/data/clinical_info.rda")

##Calculating correlation coefficient####
load("~/Documents/01-Work/06-BD_project/BD_projects/data/expr.count.rda")
BDpatient_matrix <- expr.count[,c(c(14:19),c(11:13))]
gene <- c("C1QA","C1QB","C1QC","FCER1A","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DRB1")
cor_all <- c()
pvalue_all <- c()
for (gene_i in gene){
  print(gene_i)
  gene_i_cor <- c()
  gene_i_pvalue <-c()
  for (j in seq(1,nrow(clinical_info))){
    print(j)
    gene_i_value <- as.numeric(BDpatient_matrix[gene_i,])
    j_test <- cor.test(gene_i_value,as.numeric(clinical_info[j,]))
    j_cor <- as.numeric(j_test$estimate)
    gene_i_cor <- c(gene_i_cor,j_cor)
    gene_i_pvalue <- c(gene_i_pvalue,j_test$p.value)
  }
  cor_all <- c(cor_all,gene_i_cor)
  pvalue_all <- c(pvalue_all,gene_i_pvalue)
}

###change shape into df:####
shape <- function(cor_vector,gene=gene,clinical_info=clinical_info){
  #input:
  #cor_vector: correlation vector generated from last step.
  #gene: gene vector is used for colnames.
  #clinical_info: row names of clinical info are used as row names of cor_df.
  #output: 
  #dataframe: rows are gene and columns are clinical index.
  cor_df <- matrix(cor_vector,ncol=8)
  colnames(cor_df) <- gene
  rownames(cor_df) <- rownames(clinical_info)
  cor_df <- t(cor_df)
  return(cor_df)
}
cor_df <- shape(cor_all,gene = gene,clinical_info = clinical_info)
pvalue_df <- shape(pvalue_all,gene = gene,clinical_info = clinical_info)


####visulization####
par(margin(0,0,0,0))
corrplot(cor_df, p.mat=pvalue_df,insig = "label_sig",
         sig.level = c(.001,.01,.05),pch.cex = .9, pch.col = "white")
