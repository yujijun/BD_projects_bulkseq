# --------------
# Date:  2019-12-23 17:34:31 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This script is for venndiagram plot
# 

library(VennDiagram)
BDup <- BDup.gene$gene
BDUup <- BDUup.gene$gene
BDVup <- BDVup.gene$gene
venn.diagram(x= list(BDup = BDup,BDUup = BDUup,BDVup = BDVup), filename = "./figure/Vennplot of Upgene.png", height = 450, width = 450,resolution =300, imagetype="png", 
             col ="transparent", fill =c("cornflowerblue","green","darkorchid1"),
             alpha = 0.5, label.col = c("orange", "white","darkorchid4", "white", "white",
                                        "white",   "darkgreen"),
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             cat.col =c("darkblue", "darkgreen", "orange"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 0) 

BDdown <- BDdown.gene$gene
BDUdown <- BDUdown.gene$gene
BDVdown <- BDVdown.gene$gene
venn.diagram(x= list(BDdown = BDdown,BDUdown = BDUdown,BDVdown = BDVdown), filename = "./figure/Vennplot of downgene.png", height = 450, width = 450,resolution =300, imagetype="png", 
             col ="transparent", fill =c("cornflowerblue","green","darkorchid1"),
             alpha = 0.5, label.col = c("orange", "white","darkorchid4", "white", "white",
                                        "white",   "darkgreen"),
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             cat.col =c("darkblue", "darkgreen", "orange"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 0) 
venn.diagram(list(A=1:10,B=3:8,C=6:9), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename="./figure/VennDiagram.tiff")
