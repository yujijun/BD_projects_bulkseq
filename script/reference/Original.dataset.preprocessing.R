# --------------
# Date:  2019-12-22 17:19:27 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: original dataset preprocessing
# 

input_raw <- read.table("./data-raw/all_compare-PBMC bulk-RNA.xls",header = T,sep = "\t",quote = "")

#original company count table#
name.select <- c(paste0("HC",seq(1,10)),paste0("BDU",seq(1,3)),paste0("BDV",seq(1,6)))
name.select <- paste(name.select,"_count",sep = "")
name.select <- c("gene_name",name.select)
expr.count.company <- input_raw[,name.select]
expr.count.company <- expr.count.company[!duplicated(expr.count.company$gene_name),]
write.table(expr.count.company,'./data-raw/expr.count.company.expression',sep = "\t",quote = F,col.names = T,row.names = F)
#original company FPKM table with gene name#
expr.fpkm.company <- input_raw[,c(55,grep("_fpkm",colnames(input_raw)))]
expr.fpkm.company <- expr.fpkm.company[!duplicated(expr.fpkm.company$gene_name),]
write.table(expr.fpkm.company,'./data-raw/expr.fpkm.company.expression.genename',sep = "\t",quote = F,col.names = T,row.names = F)

#original company FPKM table with gene ID#
expr.fpkm.company <- input_raw[,c(1,grep("_fpkm",colnames(input_raw)))]
expr.fpkm.company <- expr.fpkm.company[!duplicated(expr.fpkm.company$gene_id),]
write.table(expr.fpkm.company,'./data-raw/expr.fpkm.company.expression.geneid',sep = "\t",quote = F,col.names = T,row.names = F)
