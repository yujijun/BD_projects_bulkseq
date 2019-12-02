task trust {
    String prefix
    File bam_file
    File bai_file
    String genome
    String? overlap_length=10
    String? insert_threshold=10
    Int? memory=8
    Int? num_threads=1
    String? docker="gcr.io/cidc-dfci/trust:v3.0.2"

    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
        cpu: "${num_threads}"
        zones: "us-east1-b us-east1-c us-east1-d"
    }
    command {

        work=$PWD
        fname="$(basename ${bam_file})"
        echo $work
        echo $fname
        trust -f ${bam_file} -o $work/ -l ${overlap_length} -I ${insert_threshold} -g ${genome} -n ${num_threads} -c -E
        mv "$fname-TCR-ALL-coverage.txt" $work/${prefix}.tcr_coverage.txt
        mv "$fname-TCR-ALL.txt" $work/${prefix}.tcr_clones.txt
        mv "$fname-TCR-ALL.fa" $work/${prefix}.tcr_clones.fa
        mv "$fname-TCR-ALL-extended.txt" $work/${prefix}.tcr_extended.txt

    }
    output {
        File tcr_coverage = "${prefix}.tcr_coverage.txt"
        File tcr_clones_tab = "${prefix}.tcr_clones.txt"
        File tcr_clones_fa = "${prefix}.tcr_clones.fa"
        File tcr_clones_ext = "${prefix}.tcr_extended.txt"
    }
    meta {
        maintainer: "Xihao Hu"
        citation: "Li B#, Li T#, Pignon JC, Wang B, Wang J, Shukla SA, Dou R, Chen Q, Hodi FS, Choueiri TK, Wu C, Hacohen N, Signoretti S, Liu JS*, Liu XS*. Landscape of tumor-infiltrating T cell repertoire of human cancers. Nat Genet 2016; 48(7):725-32."
    }
}

workflow TRUST {
    String prefix
    File bam_file
    File bai_file
    String genome
    Int? memory=8
    Int? num_threads=1

    call trust {
       input: prefix=prefix,
              bam_file=bam_file,
              bai_file=bai_file,
              genome=genome,
              memory=memory,
              num_threads=num_threads
    }
}

