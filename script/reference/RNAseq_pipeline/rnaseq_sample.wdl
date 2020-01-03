import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/star_align.wdl" as star_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/rsem_workflow.wdl" as rsem_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/fastqc.wdl" as fastqc_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/salmon.wdl" as salmon_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/trust.wdl" as trust_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/BaiIndex.wdl" as index_wdl

workflow RNA_Seq{
	Map[String,File] fastqs


	###running separate samples
	call fastqc_wdl.run_fastqc as fastqc1{
		input:
			fastq = fastqs['fastq1'],
			prefix = fastqs['prefix']+'_1',
			type = 'left'
	}

	call fastqc_wdl.run_fastqc as fastqc2{
                input:
                        fastq = fastqs['fastq2'],
                        prefix = fastqs['prefix']+'_2',
                        type = 'right'
        }

	call star_wdl.star as star {
		input:
			fastq1 = fastqs['fastq1'],
			fastq2 = fastqs['fastq2'],
			prefix = fastqs['prefix'],
	}

	call rsem_wdl.rsem as rsem {
		input:
			transcriptome_bam = star.transcriptome_bam,
			prefix = fastqs['prefix'],

			## Default
			IndexPrefix="rsem_reference",
			paired_end = "true",
			rspd = "true"
	}

	call salmon_wdl.salmon as salmon{
		input:
			prefix = fastqs['prefix'],
			fastq1 = fastqs['fastq1'],
			fastq2 = fastqs['fastq2']
	}
 	
	call index_wdl.trust_index as index{
		input:
			bam_file = star.bamFile,
			prefix = fastqs['prefix']	
	}

	call trust_wdl.trust as trust{
		input:
			prefix = fastqs['prefix'],
			bam_file = index.sorted_bam,
			bai_file = index.bai_file
	} 

	meta {
		maintainer: "Jin Wang, Jingxin Fu"
		citation1:"Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21."
		citation2:"Bo Li.RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011,12:323"
		citation3: "Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8."
	}
}
