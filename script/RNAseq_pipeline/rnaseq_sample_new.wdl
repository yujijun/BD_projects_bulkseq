import "/cluster/asahu/project/cidc/mda_broad_comparison/broad_run/samtofastq.wdl" as samtofastq_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/fastqc.wdl" as fastqc_wdl
import "/cluster/cidc/cidc-pipelines/wdl/rnaseq/star_align.wdl" as star_wdl
import "/cluster/cidc/cidc-pipelines/wdl/rnaseq/rsem_workflow.wdl" as rsem_wdl
import "/cluster/xmwang/CDK4/trust_clonality.wdl" as trust_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/BaiIndex.wdl" as index_wdl

import "/cluster/cidc/cidc-pipelines/wdl/rnaseq/DESeq2.wdl" as DESeq2_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/ssgsea.wdl" as ssgsea_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/cibersort.wdl" as cibersort_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/convert.wdl" as convert_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/TIMER.wdl" as timer_wdl

import "/cluster/jwang/PathSeq/picard.wdl" as picard_wdl
import "/cluster/jwang/PathSeq/pathseq.wdl" as pathseq_wdl



workflow RNA_Seq{
	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	String prefix

	##SAMTOFASTQ requirements
    Int samtofastq_memory
    Int samtofastq_disk_space
    Int samtofastq_num_threads
    Int samtofastq_num_preempt

	##STAR requirements
	File starIndex
	Int star_memory
	Int star_disk_space
	Int star_num_threads
	Int star_num_preempt

	##RSEM requirements
	File rsem_reference
	String IndexPrefix
	Int rsem_memory
	Int rsem_disk_space
	Int rsem_num_threads
	Int rsem_num_preempt

	##DESeq2 requirement
	String type
	File annot
	Int DESeq2_num_preempt

	#ssgsea requirement
	String input_gmt
	#convert requirement
	File gencode_annotation
	#TIMER requirement
	String cancer

	
	###running separate samples
	scatter (sample in inputSamples) {

		call samtofastq_wdl.samtofastq {
            input: 
            	input_bam=sample[1], 
            	prefix=sample[0],
            	memory=samtofastq_memory,
            	disk_space=samtofastq_disk_space,
            	num_threads=samtofastq_num_threads,
            	num_preempt=samtofastq_num_preempt

        }


		call fastqc_wdl.run_fastqc as tumor_fastqc1 {
		 	input:
		 		fastq = samtofastq.fastq1,
		 		prefix = sample[0]+'_1',
		 		type = 'left'
		 }
		
		call fastqc_wdl.run_fastqc as tumor_fastqc2 {
		 	input:
		 		fastq = samtofastq.fastq2,
		 		prefix = sample[0]+'_2',
		 		type = 'right'
		 }
		
		call star_wdl.star as star_align {
			input:
				fastq1 = samtofastq.fastq1,
				fastq2 = samtofastq.fastq2,
				starIndex = starIndex,
				prefix = sample[0],
				memory = star_memory,
				disk_space = star_disk_space,
				num_threads = star_num_threads,
				num_preempt = star_num_preempt

		}

		call rsem_wdl.rsem as rsem_expression {
			input:
				transcriptome_bam = star_align.transcriptome_bam,
				rsem_reference = rsem_reference,
				IndexPrefix = IndexPrefix,
				prefix = sample[0],
				memory = rsem_memory,
				disk_space = rsem_disk_space,
				num_threads = rsem_num_threads,
				num_preempt = rsem_num_preempt,
				paired_end = "true",
				rspd = "true"
		}

		call index_wdl.trust_index as index{
		input:
			bam_file = star_align.bamFile,
			prefix = sample[0]	
		}

		call trust_wdl.trust as trust{
		input:
			prefix = sample[0],
			bam_file = index.sorted_bam,
			bai_file = index.bai_file
		} 

		call picard_wdl.picard as PICARD{
    		input:
    			fastq1 = samtofastq.fastq1,
    			fastq2 = samtofastq.fastq2,
    			sample_name = sample[0]
    	}

    	call pathseq_wdl.PathseqPipeline as PATHSEQ{
    		input:
    			input_bam = PICARD.PathseqBAM,
				sample_name = sample[0]
    	}


	}

#run Mulitple samples

	call DESeq2_wdl.DeSeq2{
		input:
			files = rsem_expression.genes,
			type = type,
			annot = annot,
			prefix = prefix,
			num_preempt = DESeq2_num_preempt
	    }

	#call combat_wdl.combat as Combat{

    #            input:
    #                    expression = rsem_genes_results,
    #                    prefix = prefix,
    #                   docker = "jingxin/combat:3.0"
    #   }
	
	call convert_wdl.convert as Convert{
		input:
			expressionfile = rsem_expression.genes,
			annotationfile = gencode_annotation,
			prefix = prefix
	}

	call cibersort_wdl.cibersort as Cibersort{
        input:
			input_tsv = Convert.convert_output,
			prefix = prefix,
			docker = "gcr.io/cidc-dfci/cibersort:1.06.1"
        }
	
	call timer_wdl.TIMER as Timer{
		input:
			ExpressionFile = Convert.convert_output,
			cancer = cancer,
			prefix = prefix,
			docker = "jinwang0611/timer2:v2"
	}

	call ssgsea_wdl.ssgsea as Ssgsea{
		input:
			input_tsv = Convert.convert_output,
			input_gmt = input_gmt,
			min_sz = 10,
			max_sz = 500,
			memory = 8,
			disk_space = 20,
			num_cpu = 4,
			boot_disk_gb = 10,
			num_preempt = 0	
	}
}
