workflow rsem_workflow {
    call rsem
}


task rsem {

    File transcriptome_bam
    File rsem_reference
    String IndexPrefix
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    Int? max_frag_len=1000
    String? paired_end
    String? rspd

    command {
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        if [[ ${paired_end} == "true" ]];then
            paired="--paired-end"
        else 
            paired=""
        fi

        if [[ ${rspd} == "true" ]];then
            rspd="--estimate-rspd"
        else 
            rspd=""
        fi

        mkdir rsem_out
        echo "rsem-calculate-expression $paired $rspd \
            --num-threads ${num_threads} \
            --bam \
            --no-bam-output \
            --fragment-length-max ${max_frag_len} \
            ${transcriptome_bam} rsem_reference/${IndexPrefix} ${prefix}"

        rsem-calculate-expression $paired $rspd \
            --num-threads ${num_threads} \
            --bam  \
            --no-bam-output \
            --fragment-length-max ${max_frag_len} \
            ${transcriptome_bam} rsem_reference/${IndexPrefix}  rsem_out/${prefix}

    }

    output {
        File genes="rsem_out/${prefix}.genes.results"
        File isoforms="rsem_out/${prefix}.isoforms.results"
    }

    runtime {
        docker: "jinwang0611/rsem:v1"
        slurm_docker: "jinwang0611/rsem:v1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        citation:"Bo Li.RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011,12:323"
    }
}


