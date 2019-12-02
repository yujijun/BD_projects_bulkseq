import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/DESeq2.wdl" as DESeq2_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/combat.wdl" as combat_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/ssgsea.wdl" as ssgsea_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/cibersort.wdl" as cibersort_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/convert.wdl" as convert_wdl
import "/cluster/jwang/RNA-Seq/WDL/test_muptiple/TIMER.wdl" as timer_wdl

workflow RNA_Seq{
	Array[File]  rsem_genes_results
	String prefix
	String input_gmt
	File gencode_annotation
	String cancer

	call DESeq2_wdl.DeSeq2 as Deseq2{
		input:
			files = rsem_genes_results,
			type = 'rsem',
			prefix = prefix,
			num_preempt = 0 
	}
	call combat_wdl.combat as Combat{

                input:
                        expression = rsem_genes_results,
                        prefix = prefix,
                        docker = "jingxin/combat:3.0"
        }
	
	call convert_wdl.convert as Convert{
		input:
			expressionfile = Combat.expression_combat,
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

	meta{
		maintainer: "Jingxin Fu, Jin Wang"
		citation1:"Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21."
		citation2:"Bo Li.RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011,12:323"
		citation3: "Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8."
	}


}
