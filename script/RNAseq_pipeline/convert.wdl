
workflow run_convert{
  call convert
}

task convert{
  File expressionfile
  File annotationfile
  String prefix

  String? docker = "jinwang0611/convertid:v1"
  Int? num_preempt=0

  command <<<
  python /convert.py --expressionfile ${expressionfile} --annotationfile ${annotationfile} --prefix ${prefix}

  >>>

  output {
  File convert_output='${prefix}_convertID.txt'
  }

  meta {
  maintainer:"Jin Wang"
  }

  runtime {
  docker:'${docker}'
  slurm_docker: "${docker}"
  disks: "local-disk 100 HDD"
  preemptible:"${num_preempt}"
  }
}
