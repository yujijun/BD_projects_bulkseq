task cibersort {
    File input_tsv
    File? mixture_file
    String prefix

    String? input_type = "gc27"
    Boolean? absolute = true
    Boolean? gencode_id_format = false
    Int? nperm = 100
    
    # Required options
    Int? memory = 16
    Int? disk_space = 50
    Int? num_cpu = 4
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    # Special case where docker MUST be specified from the user because it cannot be distributed.
    String docker
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
    }
    command {
       if [ -f "${mixture_file}" ]
       then
         MYMIX=" --mixture_file ${mixture_file} "
       else
         MYMIX=" "
       fi
       if [ "${input_type}"==" " ]
       then
           MYTYPE=" "
       else
           MYTYPE=" --input_type ${input_type} "
       fi
       CIBERSORT --tsv_in --tsv_out $MYMIX ${input_tsv} \
            -o "${prefix}.cibersort.tsv" $MYTYPE ${true=' --absolute True ' false=' --absolute False ' absolute} \
            ${true=' --gencode_id_format ' false=' ' gencode_id_format} \
            --specific_tempdir mytemp
    }
    output {
       File output_tsv = "${prefix}.cibersort.tsv"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation1: "Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, Hoang CD, Diehn M, Alizadeh AA. Robust enumeration of cell subsets from tissue expression profiles.  Nat Methods. 2015 May;12(5):453-7. doi: 10.1038/nmeth.3337. Epub 2015 Mar 30. PubMed PMID: 25822800; PubMed Central PMCID: PMC4739640."
    }
}

workflow run_cibersort {
    call cibersort {
    }
}
