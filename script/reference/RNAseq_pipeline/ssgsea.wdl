# Executes the ssgsea method and normalization from Barbie et. al.
#     via the GSVA bioconductor package
task ssgsea {
    File input_tsv
    File input_gmt
    String? prefix = "output"

    String? method = "ssgsea"
    String? kcdf = "Gaussian"
    Boolean? abs_ranking = false
    Int? min_sz
    Int? max_sz
    Int? parallel_sz = 0
    String? parallel_type = "SOCK"
    Boolean? mx_diff = true
    Float? tau
    Boolean? ssgsea_norm = true
    
    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    String? docker = "vacation/gsva:1.0.5"
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
       GSVA --tsv_in --tsv_out --gmt ${input_gmt} ${input_tsv} \
            -o "${prefix}.ssgsea.tsv" --method ${method} --kcdf ${kcdf} \
            ${true=' --abs_ranking ' false=' ' abs_ranking} \
            ${"--min_sz " + min_sz} ${"--max_sz " + max_sz} \
            --parallel_sz ${num_cpu} \
            ${true=' --mx_diff True ' false=' --mx_diff False ' mx_diff} \
            ${" --tau " + tau} \
            ${true=' --ssgsea_norm True ' false=' --ssgsea_norm False ' ssgsea_norm} \
            --specific_tempdir mytemp --verbose
       rm -r mytemp/
    }
    output {
       File output_tsv = "${prefix}.ssgsea.tsv"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation1: "Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462(7269):108-112. doi:10.1038/nature08460."
        citation2: "GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013 Jan 16;14:7. doi:10.1186/1471-2105-14-7. PubMed PMID: 23323831; PubMed Central PMCID: PMC3618321."
    }
}

workflow run_ssgsea {
    call ssgsea {
    }
}
