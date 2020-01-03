#!/bin/bash
cfg="/cluster/jwang/RNA-Seq/WDL/test_muptiple/cromwell_slurm-hack.cfg"
#wdl="/cluster/cidc/utils/slurm_demo.wdl"
wdl="/cluster/jwang/RNA-Seq/WDL/test_muptiple/rnaseq_sample.wdl"
jar="/cluster/soft/cromwell/cromwell-30.2.jar"

java -Dconfig.file=${cfg} -jar $jar run $wdl -i rnaseq_sample_input.json
