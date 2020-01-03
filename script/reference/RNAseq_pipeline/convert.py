import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--expressionfile', type=str, help='expression file')
parser.add_argument('--annotationfile', type=str, help='annotation file')
parser.add_argument('--prefix', type=str, help='output prefix')


args = parser.parse_args()


##input
# filename = "./Miao_Cohort.expression"
# annotation = "./Gencode.v27.annotation.genes.csv"
filename = args.expressionfile
annotation = args.annotationfile
expression = pd.read_csv(filename,sep = "\t")
gencode = pd.read_csv(annotation)

##generate dictionary for ensemble ID and gene symbol
annot_dict = gencode.set_index('gene_id').to_dict()['gene_name']
##convert
expression['gene_name'] = expression['gene_id'].apply(lambda x: annot_dict[x])

##output
cols = expression.shape[1]
transfered_expression = expression[list(expression.columns[(cols-1):(cols)]) + list(expression.columns[1:(cols-1)])]
##mean the rows with same gene
transfered_redup = transfered_expression.groupby('gene_name', as_index=False)[list(expression.columns[1:(cols)])].mean()
transfered_redup.to_csv(args.prefix+'_convertID.txt',index = False,sep = "\t")
