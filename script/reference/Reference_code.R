#
library(methods)
library(edgeR)
options(digits=4, width=190)

# Load gene-level and transcription-level annotations
geneAnno = read.table("GENOME_data/Homo_sapiens.GRCh38.86.gene.bed", header = FALSE, sep = "\t")
names(geneAnno) = c("chr","left","right","gid","score","strand","name","source","biotype")
isoformAnno = read.table("GENOME_data/Homo_sapiens.GRCh38.86.transcript.bed", header = FALSE, sep = "\t")
names(isoformAnno) = c("chr","left","right","tid","score","strand","name","source","biotype")

# Load gene expression data
header = c("GM12878_1","GM12878_2","K562_1","K562_2")
geneData = read.table("RNASEQ_data/edgeR.genes.rsem.txt", header = FALSE, sep = "\t", row.names = 1)
names(geneData) = header
isoformData = read.table("RNASEQ_data/edgeR.isoforms.rsem.txt", header = FALSE, sep = "\t", row.names = 1)
names(isoformData) = header

# Make DGELists
celllines = as.factor(c("GM12878","GM12878","K562","K562"))
geneList = DGEList(counts=round(geneData), genes=rownames(geneData), group = celllines)
isoformList = DGEList(counts=round(isoformData), genes=rownames(isoformData), group = celllines)

# Update DGELists gene annotation
geneList$genes = merge(geneAnno, geneList$genes, by.x = "gene_id", by.y = "genes")

# filter lowly expressed genes/transcripts and recompute the library sizes
geneList = geneList[rowSums(cpm(geneList) > 1) >= 2, , keep.lib.sizes=FALSE]
isoformList = isoformList[rowSums(cpm(isoformList) > 1) >= 2, , keep.lib.sizes=FALSE]

# Apply  TMM normalization
geneList = calcNormFactors(geneList, method="TMM")
isoformList = calcNormFactors(isoformList, method="TMM")

# Create design matrix
design = model.matrix(~celllines)
rownames(design) = header

# Data exploration, MD plot
pdf("mdplot.pdf", width=10, height=10, pointsize=12)
par(mfrow = c(2, 2)) # 2 rows, 2 columns
for (i in 1:ncol(geneList)) {
  plotMD(cpm(geneList, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(geneList)[i])
  abline(h=0, col="red", lty=2, lwd=2)
}

par(mfrow = c(2, 2)) # 2 rows, 2 columns
for (i in 1:ncol(geneList)) {
  plotMD(cpm(isoformList, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(isoformList)[i])
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

# Data exploration, MDS plot
pdf("mdsplot.pdf", width=12, height=6, pointsize=12)
colors = c("blue", "red")
par(mfrow = c(1, 2))
plotMDS(geneList, cex = 1, labels = c(1,2,1,2), col = colors[celllines], main = "Gene-level")
legend("topleft", legend = levels(celllines), col = colors, pch = c(15,16), ncol = 2)
plotMDS(isoformList, cex = 1, labels = c(1,2,1,2), col = colors[celllines], main = "Transcript-level")
legend("topleft", legend = levels(celllines), col = colors, pch = c(15,16), ncol = 2)
dev.off()

# Estimating dispersion
geneList = estimateDisp(geneList, design, robust=TRUE)
isoformList = estimateDisp(isoformList, design, robust=TRUE)

# Differential expression: quasi-likelihood F-test
geneQLF = glmQLFit(geneList, design, robust=TRUE)
geneQLFT = glmQLFTest(geneQLF)

isoformQLF = glmQLFit(isoformList, design, robust=TRUE)
isoformQLFT = glmQLFTest(isoformQLF)

# Print top 10 DE genes
#topTags(geneQLFT)

# Print top 10 DE transcripts
#topTags(isoformQLFT)

is.de.gene = decideTestsDGE(geneQLFT, p.value=0.01)
is.de.isoform = decideTestsDGE(isoformQLFT, p.value=0.01)

# Print number of DE genes
#summary(is.de.gene)

# Print number of DE transcripts
#summary(is.de.isoform)

# Average logCPM vs. logFC
pdf("plotsmear.pdf", width=12, height=6, pointsize=12)
par(mfrow = c(1, 2))
plotSmear(geneQLFT, de.tags=rownames(geneQLFT)[is.de.gene != 0], main = "Gene-level")
plotSmear(isoformQLFT, de.tags=rownames(isoformQLFT)[is.de.isoform != 0], main = "Transcript-level")
dev.off()

# Gnerate tabular output
geneDE = as.data.frame(topTags(geneQLFT, n = nrow(geneQLFT)))
geneDE = geneDE[order(geneDE$chr, geneDE$left), c(2:4,1,5:14)]
write.table(geneDE, file="DE_analysis.gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

isoformDE = as.data.frame(topTags(isoformQLFT, n = nrow(isoformQLFT)))
isoformDE = isoformDE[order(isoformDE$chr, isoformDE$left), c(2:4,1,5:14)]
write.table(isoformDE, file="DE_analysis.transcript.txt", sep = "\t", quote = F, row.names = F, col.names = T)
