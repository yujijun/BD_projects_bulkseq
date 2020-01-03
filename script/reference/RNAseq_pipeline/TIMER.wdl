workflow run_TIMER {
  call TIMER
  }

task TIMER {
  File ExpressionFile  #####input expression data(tab split format,samples as colnames,gene as rownames)
  String cancer   #### cancer type(named as TCGA, such as BRCA, LUAD, STAD)
  String prefix   #### output prefix
  String docker

  Int? num_preempt=0

  command <<<
  set -euo pipefail # bash strict mode

  mkdir TIMER_out
  cd TIMER_out
  Rscript /TIMER_ESTIMATE.R -f ${ExpressionFile} -c ${cancer} -p ${prefix} -o ./

  >>>

  output {
    File TIMER_output='TIMER_out/${prefix}_timer.txt'
    }

  meta {
    maintainer:"Jin Wang"
    citation: "Li B, Severson E, Pignon JC, Zhao H, Li T, Novak J, Jiang P, Shen H, Aster JC, Rodig S, Signoretti S, Liu JS*, Liu XS*. Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biol 2016; 17:174."
    }

  runtime {
    docker:"${docker}"
    slurm_docker:"${docker}"
    disks: "local-disk 100 HDD"
    preemptible:"${num_preempt}"
    }
  
} 

