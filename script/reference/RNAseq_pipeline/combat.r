#################################################################
# Function: Remove batch effect
# Call: Rscript combat.r -e expression_dat -b batch_dat
# R packages used: sva,optparse
# Authors: Jingxin Fu
# Last update: 2018-02-20, Jingxin Fu
#################################################################

## install necessary libraries
p = c("sva","optparse")
for(el in p){
  if (!is.element(el, installed.packages()[,1]))install.packages(el, dep=TRUE,repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
}
# clean env
rm(list=ls())
# parse parameters
# parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
  make_option(c("-e", "--expression_dat"), type="character",
              help="Input path of expression file. [Required]"),
  make_option(c("-b", "--batch_dat"), type="character",
              help="Input path of corresponding batch file[Required]"),
  make_option(c("-o","--output",type="character", help="Output files [Required]"))
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')

# get all tmp value for rsem result
getExpr <- function(dirs,column = c('expected_count','TPM','FPKM')){
    allFiles = dir(dirs)
    sampleName = sapply(allFiles,function(x){
        return(strsplit(x,split='\\.')[[1]][1])
    })

    for(f in allFiles){
        expr = read.table(paste0(dirs,'/',f),header=T,sep='\t',stringsAsFactors=F)
        rownames(expr) = expr[,1]
        expr = expr[,column,drop=F]
        #print(head(expr))
        if(which(allFiles == f)>1){
            overlapGene = intersect(rownames(allExprs),rownames(expr))
            allExprs = cbind(allExprs[overlapGene,],expr[overlapGene,])
            rownames(allExprs) = overlapGene
        }else{
            allExprs = expr
        }

    }
    colnames(allExprs) = sampleName
    return(allExprs)
}

ssgsvaFormat <- function(dat){
    dat = cbind(gene_id=rownames(dat),dat)
    return(dat)
}

writeDF <- function(dat,path){
    write.table(dat,path,quote=F,sep='\t',row.names=F)
}
# load data
if(file_test('-f',opts$expression_dat)){
    expr.dat = read.table(opts$expression_dat,sep='\t',header=T,stringsAsFactors=F)
    } else {
    expr.dat = getExpr(opts$expression_dat)
}

## do batch or not
if(file_test('-f',opts$batch_dat)){
    batch.dat = read.table(opts$batch_dat,sep='\t',header=T,stringsAsFactors=F)

    # Get overlap samples
    overlap_sample = intersect(colnames(expr.dat),rownames(batch.dat))
    # stop at low sample numbers
    if(length(overlap_sample) < dim(expr.dat)[2]/2) stop("too few samples")
    print('Processing samples:')
    print(overlap_sample)
    # remove batch effect
    errflag = F
    filter = apply(expr.dat, 1, var)
    exprZero = expr.dat[filter==0,]
    expr.dat = expr.dat[filter!=0,]
    expr.combat = tryCatch(ComBat(expr.dat[,overlap_sample],
                                  batch.dat[overlap_sample,1]),
                           error = function(e)errflag <<- T)
    expr.dat = rbind(expr.dat,exprZero)

    if(errflag){
        print('Cannot do Combat!')
        writeDF(ssgsvaFormat(expr.dat),paste0(opts$output,'.expression'))
    } else {
        expr.combat = rbind(expr.combat,exprZero)
        writeDF(ssgsvaFormat(expr.combat),paste0(opts$output,'.expression'))
    }

} else{
    print('Did not do Combat!')
    writeDF(ssgsvaFormat(expr.dat),paste0(opts$output,'.expression'))
}
