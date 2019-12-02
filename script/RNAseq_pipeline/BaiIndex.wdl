workflow BaiIndex{
    call trust_index
}

task trust_index {
    File bam_file
    String prefix
    Int? memory=8
    String? docker="biocontainers/samtools:v1.3_cv2"

    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
    }
    command {
	samtools sort ${bam_file} > ${prefix}.sorted.bam
	samtools index ${prefix}.sorted.bam > ${prefix}.sorted.bam.bai
    }

    output {
	File sorted_bam = "${prefix}.sorted.bam"
        File bai_file = "${prefix}.sorted.bam.bai"
    }

    meta {
        maintainer: "Jin Wang"
        citation: "Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup; The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352"
    }
}
