# --------------
# Date:  2019-12-04 05:37:21 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is a script for RNAseq pipeline
# install package 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("limma")
# BiocManager::install("Glimma")
#BiocManager::install("org.Hs.eg.db")

library(DESeq2)
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(methods)
library(edgeR)
library(limma)
library(gplots)
library(org.Hs.eg.db)
library(RColorBrewer)
library(statmod)
library("BiocParallel")
register(MulticoreParam(4))
options(stringsAsFactors = F)
setwd("/Users/yujijun/Documents/01-Work/06-BD_project/BD_projects")
####Method one####
#merge and annotation:convert ID to genename/merge all same genes by mean/add colnames/output expression matrix#
#this is for 'count expression
geneAnno <- read.table("./data-raw/Gencode.v27.annotation.genes.csv", header = T, sep = ",")
geneCountData <-  read.table("./data-raw/edgeR.genes.rsem.txt", header = FALSE, sep = "\t")

expr.data.annot = merge(geneAnno[,c('gene_name','gene_id')],
                        geneCountData,by.x='gene_id',by.y="V1")
expr.data.annot = expr.data.annot[,-1]
expr.data.annot <- aggregate(.~ gene_name,expr.data.annot,mean) 
header <- c("gene_name","BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6","HC10","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9")
colnames(expr.data.annot) <- header
expr.data.annot <- expr.data.annot[,c("gene_name","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10","BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6")]
write.table(expr.data.annot,'./data-raw/BD_count.expression',sep = "\t",quote = F,col.names = T)

#All for TPM expression
geneTPMData <-  read.table("./data-raw/edgeR.genes.TPM.rsem.txt", header = FALSE, sep = "\t")

expr.data.annot = merge(geneAnno[,c('gene_name','gene_id')],
                        geneTPMData,by.x='gene_id',by.y="V1")
expr.data.annot = expr.data.annot[,-1]
expr.data.annot <- aggregate(.~ gene_name,expr.data.annot,mean) 
header <- c("gene_name","BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6","HC10","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9")
colnames(expr.data.annot) <- header
expr.data.annot <- expr.data.annot[,c("gene_name","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10","BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6")]
write.table(expr.data.annot,'./data-raw/BD_TPM.expression',sep = "\t",quote = F,col.names = T)

#load the count table and sample info#
expr.data.annot=read.delim('./data-raw/BD_count.expression',header = T)
rownames(expr.data.annot) <- expr.data.annot$gene_name
expr.data.annot <- expr.data.annot[,-1]

sample.info <- data.frame(sample = colnames(expr.data.annot))
sample.info <- sample.info %>% 
  mutate(patient =  str_remove(sample,pattern = '\\d+')) %>% 
  mutate(group =  case_when(patient=='HC' ~ 'HC',
                            TRUE ~ 'BD')) %>%
  column_to_rownames(var = 'sample')

#DEgene identification for BDUvsHC and BDVvsHC
dds <- DESeqDataSetFromMatrix(countData =round(as.matrix(expr.data.annot)),
                              colData = sample.info,
                              design = ~ patient)
#pre-filtering: perform a minimal pre-filetering to keep only rows that have at least 10 reads total.
dds <- dds[rowSums(counts(dds)) >= 10,]
dds$patient <- relevel(dds$patient, ref = "HC") #change the foctor levels

dds<-estimateSizeFactors(dds) #
all(rownames(sample.info) == colnames(expr.data.annot)) #check
dds <- DESeq(dds,parallel = T)
res_BDUvsHC <- results(dds, contrast=c("patient","BDU","HC"))
res_BDVvsHC <- results(dds, contrast=c("patient","BDV","HC"))
plotMA(res_BDUvsHC)
plotMA(res_BDVvsHC)
degene=as.data.frame(res_BDUvsHC)
degene= degene[order(degene$padj),]
write.csv(degene,file ='./data_output/DEgene_BDUvsHC.csv')
degene=as.data.frame(res_BDVvsHC)
degene= degene[order(degene$padj),]
write.csv(degene,file ='./data_output/DEgene_BDVvsHC.csv')


#DEgene identification for BDvsHC
dds <- DESeqDataSetFromMatrix(countData =round(as.matrix(expr.data.annot)),
                              colData = sample.info,
                              design = ~ group)
#pre-filtering: perform a minimal pre-filetering to keep only rows that have at least 10 reads total.
dds <- dds[rowSums(counts(dds)) >= 10,]
dds$group <- relevel(dds$group, ref = "HC") #change the foctor levels

dds<-estimateSizeFactors(dds) #reason???
all(rownames(sample.info) == colnames(expr.data.annot)) #check
dds <- DESeq(dds,parallel = T)
res_BDvsHC <- results(dds, contrast=c("group","BD","HC"))

degene=as.data.frame(res_BDvsHC)
degene= degene[order(degene$padj),]
write.csv(degene,file ='./data_output/DEgene_BDvsHC.csv')

#log fold change shrinkage for visualization and ranking 
#visualization for DEgene
plotMA(res_BDvsHC,main="MAplot for BDvsHC") 
plotMA(res_BDUvsHC,main = "MAplot for BDUvsHC")
plotMA(res_BDVvsHC,main = "MAplot for BDVvsHC")

#summary for all three DEgenes
summary(res_BDvsHC)
summary(res_BDUvsHC)
summary(res_BDVvsHC)

library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("patient","group")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,main = "heatmap of top20 count")

library("RColorBrewer")
vsd <- vst(dds, blind = TRUE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- vsd$patient
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,main = "Heatmap of all sample's distance")
#pca
pcaData <- plotPCA(vsd, intgroup = c( "patient"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = patient)) +
  geom_point(size =2) +
  labs(title = "PCA of gene expression in all samples") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  ggrepel::geom_label_repel(aes(label=name),data = pcaData ) + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 10) ) + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) 
  
ggsave("./figure/PCA of gene expression in all samples.png", height = 10,width = 10)

#Enrichment analysis in DEgenes
library(clusterProfiler)
up.gene = degene %>% 
  rownames_to_column(var = 'gene') %>%
  filter(log2FoldChange>1, padj<0.01) %>% 
  pull(gene)
down.gene = degene %>% 
  rownames_to_column(var = 'gene') %>%
  filter(log2FoldChange<(-1), padj<0.01) %>% 
  pull(gene)
eg = bitr(up.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
mkk <- enrichKEGG(gene = eg$ENTREZID,organism = 'hsa')
p1 <- dotplot(mkk,
        title = 'Upregulated KEGG pathway in BD patients')
p1 + theme(plot.title = element_text(face = "bold",hjust = 0.5,size = 15))
eg.d = bitr(down.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
mkk.d <- enrichKEGG(gene = eg.d$ENTREZID,organism = 'hsa')
p2 <- dotplot(mkk.d,title = 'Downregulated KEGG pathway in BD patients')
p2 + theme(plot.title = element_text(face = "bold",hjust = 0.5,size = 15))

#pathway details
library("pathview")
pathview(gene.data  = mkk,
         pathway.id = "hsa04610",
         species    = "hsa")
browseKEGG(mkk, 'hsa04610')

#volcalno plot
degene=as.data.frame(res_BDvsHC)
#degene=degene[rownames(degene)%in% c(up.gene,down.gene),]
degene<- degene %>% 
  rownames_to_column(var = 'Gene') %>%
  mutate(Threshold = case_when(log2FoldChange>1 &padj<0.01 ~ 'Up',
                               log2FoldChange<1 &padj<0.01 ~ 'Down',
                               TRUE ~ 'No_sig')) %>%
  mutate(marked = case_when(abs(log2FoldChange)>10&padj<0.01 ~ 'marked',
                            TRUE ~ 'No_marked'))
library(ggrepel)
p=ggplot(degene,aes(x=log2FoldChange,y=-log10(padj),fill=Threshold,color=Threshold))
p+geom_point(alpha=0.6)+scale_color_manual(values=c("blue", "grey","red"))+
  labs(x="log2 (fold change)",y="-log10 (p.adj)")+
  geom_text_repel(data=degene[degene$marked =='marked',],
                  aes(label=Gene))

degene=degene[order(degene$log2FoldChange),]
degene[1:17,1]

####method two:edgeR filter out and TMM normalization####
#merge and annotation:make DEGLists, filter lowly expressed genes,
# Load gene expression data
geneAnno <- read.table("./data-raw/Gencode.v27.annotation.genes.csv", header = T, sep = ",")
geneCountData = read.table("./data-raw/edgeR.genes.rsem.txt", header = FALSE, sep = "\t", row.names = 1)
header <- c("BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6","HC10","HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9")
names(geneCountData) = header
geneCountData <- geneCountData[,c("HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10","BDU1","BDU2","BDU3","BDV1","BDV2","BDV3","BDV4","BDV5","BDV6")]

#consider patient seperatly: Make DGELists (normalization with TMM)
patient = as.factor(c(rep("HC",10),rep("BDU",3),rep("BDV",6)))
patient = relevel(patient, ref = c("HC"))
geneList = DGEList(counts=round(geneCountData), genes=rownames(geneCountData), group = patient)
# Update DGELists gene annotation
geneList$genes = merge(geneAnno, geneList$genes, by.x = "gene_id", by.y = "genes")
# filter lowly expressed genes/transcripts and recompute the library sizes
geneList = geneList[rowSums(cpm(geneList) > 1) >= 2, , keep.lib.sizes=FALSE] 
# Apply  TMM normalization
geneList = calcNormFactors(geneList, method="TMM")

# Create design matrix
design = model.matrix(~patient)
colnames(design) <- levels(patient)
# Estimating dispersion
geneList = estimateDisp(geneList, design, robust=TRUE)
# Differential expression: quasi-likelihood F-test
geneQLF = glmQLFit(geneList, design, robust=TRUE)
BDUvsHC = glmQLFTest(geneQLF,coef = 2)
BDVvsHC = glmQLFTest(geneQLF,coef = 3)

geneDE = as.data.frame(topTags(BDVvsHC, n = nrow(BDVvsHC)))
write.table(geneDE, file="./data_output/BDVvsHC.DEgene.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
geneDE.small = geneDE[,c(4,13:17)]
write.table(geneDE.small, file="./data_output/BDVvsHC.DEgene.small.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#consider patient together:
patient = as.factor(c(rep("HC",10),rep("BD",9)))
patient = relevel(patient, ref = c("HC"))
geneList = DGEList(counts=round(geneCountData), genes=rownames(geneCountData), group = patient)
geneList$genes = merge(geneAnno, geneList$genes, by.x = "gene_id", by.y = "genes")
# filter lowly expressed genes/transcripts and recompute the library sizes
geneList = geneList[rowSums(cpm(geneList) > 1) >= 2, , keep.lib.sizes=FALSE] #use rowsum as the finil expression value
# Apply  TMM normalization
geneList = calcNormFactors(geneList, method="TMM")
# Create design matrix
design = model.matrix(~patient)
colnames(design) <- levels(patient)
# Estimating dispersion
geneList = estimateDisp(geneList, design, robust=TRUE)
geneQLF = glmQLFit(geneList, design, robust=TRUE)
BDvsHC = glmQLFTest(geneQLF,coef = 2)
# Gnerate tabular output
geneDE = as.data.frame(topTags(BDvsHC, n = nrow(BDvsHC)))
write.table(geneDE, file="./data_output/BDvsHC.DEgene.all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
geneDE.small = geneDE[,c(4,13:17)]
write.table(geneDE.small, file="./data_output/BDvsHC.DEgene.small.txt", sep = "\t", quote = F, row.names = F, col.names = T)



#visualization of DEgene
pdf("./figure/plotsmear.pdf", width=12, height=6, pointsize=12)
plotSmear(BDvsHC, de.tags=rownames(BDvsHC)[is.de.gene != 0], main = "Gene-level")
dev.off()
#Print top 10 DE genes
topTags(BDvsHC)
is.de.gene = decideTestsDGE(BDvsHC, p.value=0.01)
# Print number of DE genes
summary(is.de.gene)

#original data exploreing
# Data exploration, MD plot???
pdf("./figure/mdplot.pdf", width=19, height=19,pointsize = 12)
par(mfrow = c(4, 5),mar=c(rep(2,4))) # 2 rows, 2 columns
for (i in 1:ncol(geneList)) {
  plotMD(cpm(geneList, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(geneList)[i])
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()
# Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples.
colors = c("blue", "red","yellow")
pdf("./figure/mdsplot.pdf", width=10, height=10,pointsize = 12)
par(mfrow = c(1, 1),mar=c(rep(2,4))) # 2 rows, 2 columns
plotMDS(geneList, top =2000, labels = colnames(geneList), col = colors[patient], main = "Gene-level")
legend("topleft", legend = levels(patient), col = colors, pch = c(15,16), ncol = 2)
dev.off()

####Method three####
#Reference:
#1.Xiaoman code
#2.https://ycl6.gitbooks.io/rna-seq-data-analysis/perform_de_analysis.html
#2.1 https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html
#3.https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#4.limma multiple design:https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#(9.3)
#5. usersguide of limma:http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/doc/usersguide.pdf
#6. DEseq: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#7. what is the different: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#8.Usage of DEseq: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
