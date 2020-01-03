workflow run_DESeq2 {
  call DeSeq2
  }

task DeSeq2 {
  Array[File] files
  String type
  File annot
  String prefix
  
  Int? num_preempt=0

  command <<<
  set -euo pipefail # bash strict mode

  mkdir deseq2_data
  mv ${sep=" " files} ./deseq2_data/ 

  mkdir deseq2_out
  Rscript /Deseq_script.r -d deseq2_data/ -t ${type} -a ${annot} -o ../deseq2_out/ -p ${prefix}

  >>>

  output {
    File DESeq_output='deseq2_out/${prefix}_DESeq2.txt'
    }

  meta {
    maintainer:"Jin Wang"
    citation: "Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8."
    }

  runtime {
    docker:'jinwang0611/deseq2:v3'
    slurm_docker: 'jinwang0611/deseq2:v3'
    disks: "local-disk 100 HDD"
    preemptible:"${num_preempt}"
    }
  
} 
