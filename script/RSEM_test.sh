#!/bin/bash
 
 #SBATCH -t 0-24:00 
 #SBATCH -c 20 #number of cpus required per task
 #SBATCH -n 10 #launch a maximum of number tasks 
 #SBATCH --mem 200G
 #SBATCH --job-name STAR_align_1
 #SBATCH -o %j.out  #file to which standard out will be written
 #SBATCH -e %j.err

 module load RSEM/1.3.1

rsem-calculate-expression --bam --no-bam-output -p 20 --paired-end  /liulab/jjyu/DB_bulkRNAseq/STAR_result/BDU1Aligned.toTranscriptome.out.bam /liulab/jjyu/Reference/RSEM_ref/hg38_gencode_rsem /liulab/jjyu/DB_bulkRNAseq/RSEM_result/BDU1 