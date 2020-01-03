task fastqc {
    String prefix
    String? type = "unspecified" # expect "left" "right" or "unspecified"
    File fastq

    # Required options
    Int? memory = 16
    Int? disk_space = 100
    Int? num_cpu = 4
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    String? docker = "biocontainers/fastqc:0.11.5"

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
       fastq: "fastq file to analyze"
    }
    command {
       mkdir output
       mkdir ready
       fastqc "${fastq}" \
              --extract \
              -o output
       HTML=$(ls ./output/*.html)
       ZIP=$(ls ./output/*.zip)
       DATA=$(ls ./output/*/fastqc_data.txt)
       SUMMARY=$(ls ./output/*/summary.txt)
       mv "$HTML" "${prefix}".fastqc.html
       mv "$ZIP" "${prefix}".fastqc.zip
       mv "$DATA" "${prefix}".fastqc.data.txt
       mv "$SUMMARY" "${prefix}".fastqc.summary.txt
       cat "${prefix}".fastqc.data.txt | \
           grep 'Total Sequences' | \
           head -n 1 | \
           cut -f 2 > total.txt
       rm -r output/
    }
    output {
        File html = "${prefix}.fastqc.html"
        File zip = "${prefix}.fastqc.zip"
        File data = "${prefix}.fastqc.data.txt"
        File summary_raw = "${prefix}.fastqc.summary.txt"
        Int sequences = read_int("total.txt")
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc"
    }
}
task post_fastqc {
   File summary_raw

   String? docker = "ubuntu:16.04"

   command <<<
       cat "${summary_raw}" | \
           awk -F "\t" '{print $2 "\t" $1}' > summary_map.txt
    >>>
   output {
        Map[String, String] summary = read_map("summary_map.txt")
    }
    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "2 GB"
        slurm_memory: "2G"
        disks: "local-disk 20 HDD"
        cpu: "1"
        preemptible: "0"
        bootDiskSizeGb: "10"
    }
}
workflow run_fastqc {
    String prefix
    String? type = "unspecified" # expect "left" "right" or "unspecified"
    File fastq

    call fastqc {
       input:
          prefix = prefix,
          type = type,
          fastq = fastq
    }
    call post_fastqc {
       input: summary_raw=fastqc.summary_raw
    }
    output {
        File html = fastqc.html
        File zip = fastqc.zip
        File data = fastqc.data
        File summary_raw = fastqc.summary_raw
        Int sequences = fastqc.sequences
        Map[String, String] summary = post_fastqc.summary
    }
}
