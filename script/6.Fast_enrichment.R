# --------------
# Date:  2020-01-06 10:50:15 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: Fast enrichment analysis for a specific gene list.
# 
require(clusterProfiler)
#load test geneset and DE gene matrix
load("~/Documents/01-Work/06-BD_project/BD_projects/data/geneset1.rda")
load("~/Documents/01-Work/06-BD_project/BD_projects/data/biggeneset.rda")
load("~/Documents/01-Work/06-BD_project/BD_projects/data/DEmatrix.rda")
#Enrichment analysis
eg = bitr(geneset1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
mkk <- enrichKEGG(gene = eg$ENTREZID,organism = 'hsa')

#output all enrichment pathway information and detail gene list of pathways.
output <- "~"
outpath.Enrichment <- paste(output,"Enrichment",sep = "/")
ifelse(!dir.exists(outpath.Enrichment),dir.create(outpath.Enrichment),FALSE)
for(i in seq(1,length(mkk$Description))){
  pathwayfilename_tmp <- mkk$Description[i]
  sink(paste(outpath.Enrichment,pathwayfilename_tmp,sep = "/"))
  print(mkk$Description[i])
  UpgeneID <- str_split(mkk$geneID[i],"/")[[1]]
  print(bitr(UpgeneID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db"))
  sink()
}
sink(paste(outpath.Enrichment,"All_pathway.txt",sep = "/"))
print(mkk$Description)
sink()

#output enrichment pathway gene list with up and down regulator information.
viral_geneset1 <- read.table(file = "./output_explore/Enrichment/Upenrichment_Viral carcinogenesis",skip =1,header = T,stringsAsFactors = F)
geneset1 <- viral_geneset1$SYMBOL
geneset1.expr <- DEmatrix[geneset1,]
geneset1.expr <- geneset1.expr[order(geneset1.expr$log2FoldChange,decreasing = T),]
geneset1.expr <- geneset1.expr[,-c(1,3,4,5)]
write.table(geneset1.expr,file = "Viral_carcinogenesis.tsv",sep = "\t",col.names = T,row.names = T,quote = F)

#If someone would like to explore more about this gene list,we can output this gene list 
#and upload them into string(http://version10.string-db.org/).
write(geneset1, file = "/Users/yujijun/Desktop/geneset6.txt",
      append = FALSE, sep = "\t")
write(allgene, file = "/Users/yujijun/Desktop/allgene.txt",
      append = FALSE, sep = "\t")

