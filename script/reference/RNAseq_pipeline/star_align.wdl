workflow star_workflow {
  call star
  }

task star {
    # star required parameters
    File fastq1
    File? fastq2
    File starIndex
    String prefix
    
    ## star optional parameters
    String? annotation_gtf
    Int? outFilterMultimapNmax=20
    Int? alignSJoverhangMin=8
    Int? alignSJDBoverhangMin=1
    Int? outFilterMismatchNmax=999
    Float? outFilterMismatchNoverLmax=0.1
    Int? alignIntronMin=20
    Int? alignIntronMax=1000000
    Int? alignMatesGapMax=1000000
    String? outFilterType="BySJout"
    Float? outFilterScoreMinOverLread=0.33
    Float? outFilterMatchNminOverLread=0.33
    String? outSAMstrandField='intronMotif'
    String? outFilterIntronMotifs='None'
    String? quantMode='TranscriptomeSAM GeneCounts'
    String? outSAMtype='BAM Unsorted'
    String? outSAMunmapped='None'
    Int? chimSegmentMin=15
    Int? chimJunctionOverhangMin=15
    String? chimOutType='WithinBAM SoftClip'
    Int? chimMainSegmentMultNmax=1
    String? genomeLoad='NoSharedMemory'
    File? sjdbFileChrStartEnd

    # system parameters

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command <<<
      set -euo pipefail # bash strict mode,avoid hidden errors 

      if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]];then
        tar -xvvf ${fastq1}
        fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done|paste -s -d ',') # return the exact match formula when cannot find any matching.
        fast2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done| paste -s -d ',')
        
        if [[ $fastq1_abs == *"*_1.fastq*" ]];then # no parted-end FASTQs found; then check single-end FASTQ
          fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
          fastq2_abs=''
        fi

      else
        #make sure paths are absolute
        fastq1_abs=${fastq1}
        fastq2_abs=${fastq2}
        if [[ $fastq1_abs != /* ]]; then
          fastq1_abs=$PWD/$fastq1_abs
          fastq2_abs=$PWD/$fastq2_abs
        fi

      fi
      echo 'FASTQs:'
      echo $fastq1_abs
      echo $fastq2_abs
      # filecmd
      case $fastq1_abs in 
		*.gz) fileCmd='zcat';;
		*.bz) fileCmd='bzcat';;
		*.fastq) fileCmd='cat';;
		*) exit 1
      esac
      
      mkdir starIndex
      tar -xvvf ${starIndex} -C starIndex --strip-components=1
      echo $(date) "--Begin running STAR"
      
      mkdir star_out
      STAR --runMode alignReads --readFilesIn $fastq1_abs $fastq2_abs --genomeDir starIndex --runThreadN ${num_threads} --outFileNamePrefix star_out/${prefix}. --readFilesCommand $fileCmd --outSAMtype ${outSAMtype} \
--chimSegmentMin ${chimSegmentMin} \
--chimJunctionOverhangMin ${chimJunctionOverhangMin} \
--chimOutType ${chimOutType} \
--chimMainSegmentMultNmax ${chimMainSegmentMultNmax} \
--outFilterMultimapNmax ${outFilterMultimapNmax} \
--outFilterMismatchNmax ${outFilterMismatchNmax} \
--outFilterMismatchNoverLmax ${outFilterMismatchNoverLmax} \
--outFilterType ${outFilterType} \
--outFilterScoreMinOverLread ${outFilterScoreMinOverLread} \
--outFilterMatchNminOverLread ${outFilterMatchNminOverLread} \
--outFilterIntronMotifs ${outFilterIntronMotifs} \
--alignSJoverhangMin ${alignSJoverhangMin} \
--alignSJDBoverhangMin ${alignSJDBoverhangMin} \
--alignIntronMin ${alignIntronMin} \
--alignIntronMax ${alignIntronMax} \
--alignMatesGapMax ${alignMatesGapMax} \
--outSAMstrandField ${outSAMstrandField} \
--quantMode ${quantMode} \
--outSAMunmapped ${outSAMunmapped} \
--genomeLoad ${genomeLoad} \ 
${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
${"--sjdbGTFfile " + annotation_gtf}
    >>>

    output {
      File bamFile='star_out/${prefix}.Aligned.out.bam'
      File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
      File chimeric_junctions = "star_out/${prefix}.Chimeric.out.junction"
      File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sam"
      File read_counts = "star_out/${prefix}.ReadsPerGene.out.tab"
      File junctions = "star_out/${prefix}.SJ.out.tab"
      Array[File] logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
      
      }
    meta {
	maintainer:"Jingxin Fu"
	citation:"Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21."
     }

    runtime{
      docker: "jinwang0611/star:v1"
      slurm_docker: "jinwang0611/star:v1"
      memory: "${memory}GB"
      disks: "local-disk ${disk_space} HDD"
      cpu: "${num_threads}"
      preemptible: "${num_preempt}"
     }

  }
