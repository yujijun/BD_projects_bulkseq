# --------------
# Date:  2020-01-03 16:43:28 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:
# 
library(sva)
library(pamr)
library(limma)

#Prepare all dataset
expr.tpm <- read.table("./data-raw/expr.tpm.company.expression.genename",header = T,sep = "\t")
rownames(expr.tpm) <- expr.tpm[,1]
expr.tpm <- expr.tpm[,-1]
expr.tpm <- expr.tpm[rowSums(expr.tpm) >= 10,]
expr.tpm <- log2(expr.tpm+1)

mod = model.matrix(~as.factor(group),data = sample.info)
mod0 = model.matrix(~1,data = sample.info)

#Next we apply function to estimate the surrogate variables:
expr.tpm <- as.matrix(expr.tpm)
svobj = sva(expr.tpm,mod,mod0)

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)

pValuesSv = f.pvalue(expr.tpm,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

fit = lmFit(expr.tpm,modSv)
contrast.matrix <- contrast.matrix <- cbind("C1"=c(-1,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
dim(topTreat(eb, adjust="BH",coef="C1"))
dim(topTable(eb, adjust="BH",coef="C1",number = Inf,p.value = 0.05,lfc=1))

#Ask some question about 
#Reference:
#1. https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#2. http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html
#3. http://127.0.0.1:30707/library/sva/doc/sva.pdf
