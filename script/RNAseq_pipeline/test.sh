#!/bin/bash
cfg="/cluster/jwang/RNA-Seq/WDL/test_muptiple/cromwell_slurm-hack.cfg"
#wdl="/cluster/cidc/utils/slurm_demo.wdl"
wdl="/cluster/jwang/RNA-Seq/WDL/test_muptiple/BaiIndex.wdl"
jar="/cluster/soft/cromwell/cromwell-30.2.jar"

java -Dconfig.file=${cfg} -jar $jar run $wdl -i BaiIndex_input.json
