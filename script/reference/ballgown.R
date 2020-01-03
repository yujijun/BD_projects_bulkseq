library(ballgown)
pheno_data <- read.csv("phenodata.csv")
bg <- ballgown(dataDir = "ballgown",
               samplePattern = "sample",
               pData = pheno_data)
samplesNames(bg)
bgfilt <-subset(bg,'rowVars(texpr(bg))>1',genomesubset=TRUE)（过滤掉表达差异较小的基因）
diff_genes <- stattest(bgfilt,feature='gene',covariate=【自变量】,adjustvars=【无关变量】,meas='FPKM')
diff_genes <- arrange(diff_genes,pval)
write.csv(diff_genes,'diff_genes.csv',row.names=FALSE)
# MA plot
library(ggplot2)
library(cowplot)
results_transcripts$mean <- rowMeans(texpr(bg_chrX_filt))

ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)