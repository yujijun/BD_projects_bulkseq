# BD_projects
This is a project about Behcet's disease
I also would like to package some normal useful packages into this projects.

#File introduction:
**/script **
#Completed computational propress under this folder.
1. User can use 1.DEseqDataset_and_batcheffect.Rmd to explore batch effect in samples.  

2. When make the batch information clear, User can generate DESeqDataSet by 2.Generate_DESeqDataSet.R  

2.5 Then User can generate normal result or figure by 2.5 Main.R.  

3. User can also calculate correlation between gene and clinical index by 3. Clinical_coranalysis.R  

4. User can generate immune cells proportion from https://cibersort.stanford.edu/ and then analysis more details by 4.CIBERSOFT.analysis.R.  

5. Then, If required, users can draw a Venn diagram of different genes in different situations。

**/data**
User can find all test data under this folder.
Reference:
1. https://ycl6.gitbooks.io/rna-seq-data-analysis/quantification_using_rsem1.html
2. http://tiramisutes.github.io/2018/12/04/ref-RNA-seq.html
3. https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html (downstream analysis）