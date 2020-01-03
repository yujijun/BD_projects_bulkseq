# --------------
# Date:  2019-12-22 16:59:07 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:This script is for original dataset pre-processing:change FPKM into TPM 

#Input dataset#
expr.data.annot=read.delim('./data-raw/expr.fpkm.company.expression.genename',header = T)

fpkm_to_tpm <- function(x){
  y <- x/sum(x,na.rm = T)*10^6
  return(y)
}

expr.data.annot.tpm <- apply(expr.data.annot[,-1],2,fpkm_to_tpm)
expr.data.annot.tmp.gene <- add_column(as.data.frame(expr.data.annot.tpm),gene_name=expr.data.annot$gene_name,.before = "HC1_fpkm")
colnames(expr.data.annot.tmp.gene) <- gsub("_fpkm","_tpm",colnames(expr.data.annot.tmp.gene))
write.table(expr.data.annot.tmp.gene,'./data-raw/expr.tpm.company.expression.genename',sep = "\t",quote = F,col.names = T,row.names = F)
