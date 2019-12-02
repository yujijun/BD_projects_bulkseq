workflow run_combat {
  call combat
  }

task combat {
  Array[File] expression
  File batch
  String prefix
  String docker
  Int? num_preempt=0
  command <<<
  set -euo pipefail # bash strict mode

  mkdir combat_out/ 
  mv ${sep=" " expression} combat_out/ 

  Rscript /combat.r -e combat_out/ -b ${batch} -o combat_out/${prefix} 
  >>>
  
  output {
    File expression_combat='combat_out/${prefix}.expression'
    }

  meta {
    maintainer:"Jingxin Fu"
    citation: "Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Storey JD, Zhang Y and Torres LC (2017). sva: Surrogate Variable Analysis. R package version 3.26.0."
    }
  runtime {
    docker:'${docker}'
    slurm_docker:'${docker}'
    preemptible:"${num_preempt}"
    }
}
