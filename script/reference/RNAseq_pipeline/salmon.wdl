task salmon {
    Array[File] index
    String prefix
    File fastq1
    File? fastq2
    File? geneMap # Annotation is optional but highly recommended

    # options use full but increase run time
    Int? seqBias
    Int? gcBias

    # Required options
    Int? memory = 16
    Int? disk_space = 100
    Int num_cpu
    Int? num_preempt = 0
    Int? boot_disk_gb = 30

    String? docker = "combinelab/salmon:0.9.1"
    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
        slurm_memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${num_preempt}"
        bootDiskSizeGb: "${boot_disk_gb}"
        zones: "us-east1-b us-east1-c us-east1-d"
    }
    parameter_meta {
       index: "list of salmon index files"
       prefix: "prefix for name of output"
       fastq1: "Fastq in gz or raw format"
       fastq2: "Paired Fastq in gz or raw format"
    }
    command {
      INDEX=$(echo $(ls ${index[1]} | head -n 1) | rev |  cut -d '/' -f 2- | rev)
      if [ -f "${fastq2}" ]
      then
        salmon quant --index $INDEX \
                   -1 "${fastq1}" \
                   -2 "${fastq2}" \
                   -o out_quant \
                   -p ${num_cpu} \
                   --libType A \
                   ${"-g " + geneMap} \
                  $(if [ 1 -eq ${seqBias} || 1 -eq ${gcBias} ]; then echo '--auxDir aux'; fi) \
                  $(if [ 1 -eq ${seqBias} ]; then echo '--seqBias'; fi) \
                  $(if [ 1 -eq ${gcBias} ]; then echo '--gcBias'; fi)
      else
        salmon quant --index $INDEX \
                   -r "${fastq1}" \
                   -o out_quant \
                   -p ${num_cpu} \
                   --libType A \
                   ${"-g " + geneMap} \
                  $(if [ 1 -eq ${seqBias} || 1 -eq ${gcBias} ]; then echo '--auxDir aux'; fi) \
                  $(if [ 1 -eq ${seqBias} ]; then echo '--seqBias'; fi) \
                  $(if [ 1 -eq ${gcBias} ]; then echo '--gcBias'; fi)
      fi      
      cp -r out_quant "${prefix}".salmon
      tar -c "${prefix}".salmon | gzip > "${prefix}".salmon.tar.gz
      touch "out_quant/quant.genes.sf" # will be an empty file if geneMap is not specified
      cp -r out_quant/aux_info "${prefix}".salmon.aux_info
      tar -c "${prefix}".salmon.aux_info | gzip > "${prefix}".salmon.aux_info.tar.gz
    }
    output {
        File salmon_tarball = "${prefix}.salmon.tar.gz"
        File aux_info_tarball = "${prefix}.salmon.aux_info.tar.gz"
        File quant = "out_quant/quant.sf"
        File quant_genes = "out_quant/quant.genes.sf"
        String cmd_info = read_string("out_quant/cmd_info.json") #read_json may be better, but cromwell says its not implemented yet when I try to use it.
        String lib_format_counts = read_string("out_quant/lib_format_counts.json")
        String aux_info_meta_info = read_string("out_quant/aux_info/meta_info.json")

        File libParams_flenDist = "out_quant/libParams/flenDist.txt"
        File logs_salmon_quant = "out_quant/logs/salmon_quant.log"

    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods."
    }
}

workflow salmon_workflow {
    call salmon
}
