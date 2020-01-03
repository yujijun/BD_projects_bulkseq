#' This is a volcano plot function
#' @param data the dataframe you would like to display
#' @param x the colname of data which you would like to use as x value.
#' @param y the colname of data which you would lite to use as y value.
#' @param fill the colname of data you would like to draw color(fill color)
#' @param color the colname of data you would like to draw color(edge color)
#' @param color_palette all colors your would like to use in your figure.
#' @param label the lable name you would like to put
#' @param label_data  A dataframe which you would like to take label
#' @param alpha the transparency of points
#' @param plotLab the logical value, if you would like to display label text.
#' @param xlab the name of x-asix, with the default string of "log2FoldChange".
#' @param ylab the name of y-asix, with the default string of "-log10 Padj".
#' @param title the name of this figure. default name is "Volcano Plot of genes"
#' @param xline the xintercept, default is c(-1,1).
#' @return a volcano plot of ggplot format
#' @export
Volcano <- function(data = data,
                    x= "log2FoldChange",y="Padj",
                    fill, color = "Threshold",
                    color_palette = c("#8DA0CB","grey","#FC8D62"),
                    label="gene",label_data = subset_data,
                    alpha=0.5,plotLab = T,
                    xlab = "log2 Fold Change",ylab="log10 Padj",
                    title="Volcano Plot of genes",xline=c(-1,1)){
  require(ggrepel)
  require(ggplot2)
  p=ggplot(data,aes_string(x=x,y=y,fill=fill,color=color))+
    geom_point(alpha=alpha)+
    scale_color_manual(values=color_palette)+
    labs(x=xlab,y=ylab,title=title)+theme_classic()+
    geom_hline(yintercept=1.3, linetype="dashed", color = "grey")+
    geom_vline(xintercept=xline, linetype="dashed", color = "grey")+
    geom_point(data = label_data, aes_string(x=x,y=y,fill=fill,color=color))+
    #geom_text_repel(label = label_data$Official)
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 15),legend.position = "right") +
    theme(legend.text = element_text(size = 10,face = "bold")) +
    theme(legend.title = element_text(size = 15,face = "bold"))
    #guides(colour = guide_legend(override.aes = list(shape = 15)))
  #geom_text(aes(label=subset_data, color=subset_data), size=10)
  if (plotLab==T) {
    p<- p+geom_text_repel(data=label_data, aes_string(label=label))
  }
  plot(p)
}
