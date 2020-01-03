library("PoiClaClu")

marker_gene <- read.table("./data-raw/human_immune_CIBERSORT.txt",header = T,stringsAsFactors = F)

poisd <- PoissonDistance(t(counts(dds_BD)))
sampleDists <- poisd$dd
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- c(paste0("HC",seq(1,10)),paste0("BDU",seq(1,3)),paste0("BDV",seq(1,6)))
rownames(sampleDistMatrix) <- c(paste0("HC",seq(1,10)),paste0("BDU",seq(1,3)),paste0("BDV",seq(1,6)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = ,
         clustering_distance_cols = ,
         col = colors,main = "Heatmap of sample poisd distance")
