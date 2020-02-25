library(GSVA)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
expr.tpm.18 <- expr.tpm[,-19]
gene.adaptive.immue <- read.table("./Reference/GO_ADAPTIVE_IMMUNE_RESPONSE.txt",sep = "\t")
gene.adaptive.immue <- as.character(gene.adaptive.immue[,1])
gene.inflammatory.response <- read.table("./Reference/GO_INFLAMMATORY_RESPONSE.txt",sep = "\t")
gene.inflammatory.response <- as.character(gene.inflammatory.response[,1])
gene.inflammatory.response.cellreceptor <- read.table("./Reference/GO_INFLAMMATORY_RESPONSE_Cell_receptor.txt",sep = "\t")
gene.inflammatory.response.cellreceptor <- as.character(gene.inflammatory.response.cellreceptor[,1])
gene.innate.immune.response <- read.table("./Reference/GO_INNATE_IMMUNE_RESPONSE.txt",sep = "\t")
gene.innate.immune.response <- as.character(gene.innate.immune.response[,1])
gene.innate.immune.inmucosa <- read.table("./Reference/GO_INNATE_IMMUNE_RESPONSE_IN_MUCOSA.txt",sep = "\t")
gene.innate.immune.inmucosa <- as.character(gene.innate.immune.inmucosa[,1])
gene.interferon.alpha <- read.table("./Reference/GO_RESPONSE_TO_INTERFERON_ALPHA.txt",sep = "\t")
gene.interferon.alpha <- as.character(gene.interferon.alpha[,1])
gene.interferon.beta <- read.table("./Reference/GO_RESPONSE_TO_INTERFERON_BETA.txt",sep = "\t")
gene.interferon.beta <- as.character(gene.interferon.beta[,1])
gene.interferon.gamma <- read.table("./Reference/GO_RESPONSE_TO_INTERFERON_GAMMA.txt",sep = "\t")
gene.interferon.gamma <- as.character(gene.interferon.gamma[,1])



#Three kinds of violineplot with specific gene position
Violinplot_singlegene <- function(LFC,output.path,figure.name,genename,width,height,violincol="#9F79EE"){
  empty <- data.frame()  #concat all dataset together
  datasets <- colnames(LFC)
  for(i in seq(1, ncol(LFC))){
    Value = LFC[,i]
    Gene = rownames(LFC)
    Dataset = rep(datasets[i],nrow(LFC))
    dataframe_i <- data.frame(Gene,Value,Dataset)
    empty <- rbind(empty, dataframe_i)
    print(dim(empty))
  }
  empty <- empty[!is.na(empty$Value),]
  p <- ggviolin(empty, x = "Dataset", y = "Value", fill = violincol,
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                add = "boxplot", add.params = list(fill = "white"),
                xlab = "Datasets",
                ylab = "LFC_value",
                title = figure.name) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,color = "#2F4F4F")) +
    theme(axis.text.x.bottom = element_text(size = 7,face = "bold")) +
    theme(axis.text.x.bottom = element_text(angle = 90,hjust=1)) +
    #scale_x_discrete(position = "") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 3,byrow = T)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 5,face = "bold")) +
    theme(axis.title.x = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
    theme(axis.title.y = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
    theme(axis.text.x = element_text(face = "bold",size = 20)) +
    theme(legend.position = "none") +
    geom_jitter(data = empty[empty$Gene==genename,],height = 0, width = 0, aes(colour = "#8B0000"))
  output_path = output.path
  figure_name = paste(figure.name,".png",sep = "")
  png(paste(output_path,figure_name,sep = "/"),height = height,width = width,units ="cm",res=150)
  print(p)
  dev.off()
}
geneset <- list(gene.adaptive.immue = gene.adaptive.immue,
                   gene.inflammatory.response = gene.inflammatory.response,
                   gene.inflammatory.response.cellreceptor = gene.inflammatory.response.cellreceptor,
                   gene.innate.immune.inmucosa = gene.innate.immune.inmucosa,
                   gene.innate.immune.response = gene.innate.immune.response,
                   gene.interferon.alpha = gene.interferon.alpha,
                   gene.interferon.beta = gene.interferon.beta,
                   gene.interferon.gamma = gene.interferon.gamma)
gsva.score <- gsva(as.matrix(expr.tpm.18), geneset,
                   kcdf="Gaussian", mx.diff=T,verbose=FALSE,
                   parallel.sz=50,ssgsea.norm=T)
gsva.score <- gsva.score[-c(3,4),]
gsva.score <- gsva.score[c(1,3,2,4,5,6),]
rownames(gsva.score) <- gsub("gene.","",rownames(gsva.score))
rownames(gsva.score)[1] <- "adaptive.immune"
empty <- data.frame()  #concat all dataset together
datasets <- rownames(gsva.score)
for(i in seq(1, nrow(gsva.score))){
  Value = gsva.score[i,]
  Sample = colnames(gsva.score)
  Dataset = rep(datasets[i],ncol(gsva.score))
  dataframe_i <- data.frame(Sample,Value,Dataset)
  empty <- rbind(empty, dataframe_i)
  print(dim(empty))
}
empty <- empty[!is.na(empty$Value),]
empty <- empty %>% mutate(patient = str_remove(empty$Sample,pattern = '\\d+.*'))
empty$patient <- gsub("BDU","BD",empty$patient)
empty$patient <- gsub("BDV","BD",empty$patient)

#comparisons <- list(c("HC","BD"))
#comparisons <- list(c("gene.adaptive.immue","gene.inflammatory.response"))
figure.name <- "GSVA Score of Different Immune Response "
p <- ggviolin(empty, x = "Dataset", y = "Value", fill = "patient",
              #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot",
              add.params = list(fill = "patient"),
              xlab = "Immune Responses",
              ylab = "GSVA_Value",
              title = figure.name) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,color = "#2F4F4F")) +
  theme(axis.text.x.bottom = element_text(size = 12,face = "bold")) +
  theme(axis.text.x.bottom = element_text(angle = 90,hjust=1)) +
  #scale_x_discrete(position = "") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3,byrow = T)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 10,face = "bold")) +
  theme(axis.title.x = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
  theme(axis.title.y = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
  theme(axis.text.x = element_text(face = "bold",size = 20)) +
  theme(legend.position = "right")  
  #stat_compare_means(comparisons = comparisons,label="p.signif")
  #geom_jitter(data = empty[empty$Gene==genename,],height = 0, width = 0, aes(colour = "#8B0000"))

print(p)

p <- ggplot(empty, aes(x = patient, y = Value,fill=patient)) +
  ggtitle("GSVA score of different immune respones") + 
  xlab("Sample_Type") +
  ylab("GSVA_Score") +
  geom_boxplot(alpha=0.7) + 
  facet_wrap(~Dataset,ncol = 3)+
  stat_compare_means(comparisons = list(c("BD", "HC"))) + 
  theme(plot.title = element_text(hjust=0.5,face ="bold",size = 15)) + 
  theme(axis.title.y = element_text(face = "bold")) + 
  theme(axis.title.x = element_text(face = "bold")) + 
  theme(strip.text.x = element_text(face = "bold",size = 10))
show(p)
