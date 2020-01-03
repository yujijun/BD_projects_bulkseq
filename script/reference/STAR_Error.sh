#!/bin/bash
 
 #SBATCH -t 0-24:00 
 #SBATCH -c 40 #number of cpus required per task
 #SBATCH -n 3 #launch a maximum of number tasks 
 #SBATCH --mem 60G
 #SBATCH --job-name STAR_align_1
 #SBATCH -o %j.out  #file to which standard out will be written
 #SBATCH -e %j.err
 
 module load STAR/2.6.1c
 
/data/programs/STAR/2.6.1c_2018-10-30/bin/Linux_x86_64/STAR --runThreadN 40 --readFilesIn /liulab/jjyu/DB_bulkRNAseq/BDV4_1.fq.gz /liulab/jjyu/DB_bulkRNAseq/BDV4_2.fq.gz --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /liulab/jjyu/DB_bulkRNAseq/STAR_result/BDV4 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --readFilesCommand zcat --genomeDir /liulab/jjyu/Reference/hg38_index &
/data/programs/STAR/2.6.1c_2018-10-30/bin/Linux_x86_64/STAR --runThreadN 40 --readFilesIn /liulab/jjyu/DB_bulkRNAseq/HC8_1.fq.gz /liulab/jjyu/DB_bulkRNAseq/HC8_2.fq.gz --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /liulab/jjyu/DB_bulkRNAseq/STAR_result/HC8 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --readFilesCommand zcat --genomeDir /liulab/jjyu/Reference/hg38_index 