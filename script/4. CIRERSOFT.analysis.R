# --------------
# Date:  2019-12-23 17:58:30 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for cibersoft analysis
#Input: sibersoft result from cibersoft website.
#output:significant plot
# 
#CIBERSOFT analysis:
ciber.data <- read.csv("./data_output/CIBERSOFT_output/CIBERSORT.Output_tpm.csv",header = T,row.names = 1)
ciber.data <- ciber.data[,-c(23,24,25)]
ciber.data.filter <- ciber.data[,apply(ciber.data, MARGIN = 2,sum) > 0.1]
ciber.data.filter$class <- c(rep("HC",10),rep("BD",9))
empty <- data.frame()
samples <- colnames(ciber.data.filter[,-12])
for(i in seq(1,ncol(ciber.data.filter)-1)){
  value <- ciber.data.filter[,i]
  celltype <- rep(samples[i],nrow(ciber.data.filter))
  class <- ciber.data.filter$class
  data.frame_i <- data.frame(value,celltype,class)
  empty <- rbind(empty,data.frame_i)
  print(dim(empty))
}
p <- ggplot(empty, aes(x = class, y = value,fill=class)) +
  ggtitle("Proportion of immune cell type") + 
  xlab("Type of patients") +
  ylab("Proportion of immune cell") +
  geom_boxplot(alpha=0.5) + 
  facet_wrap(~celltype,ncol = 6)+
  stat_compare_means(comparisons = list(c("BD", "HC"))) + 
  theme(plot.title = element_text(hjust=0.5,face ="bold",size = 15)) + 
  theme(axis.title.y = element_text(face = "bold")) + 
  theme(axis.title.x = element_text(face = "bold")) + 
  theme(strip.text.x = element_text(face = "bold",size = 10))
show(p)
ggsave("./figure/TPM_Proprotion_immue.png",height = 8,width = 15)
