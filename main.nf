#! /usr/bin/env nextflow

/* Name: RNASeq pipeline
 * Auth: Sharu Paul
 * Desc: NextFlow dsl2 pipeline for RNASeq analysis
 */

// Help statement
def helpMsg() {
  log.info """\
          RNA-SEQ PIPELINE
     ============================
   Usage:
      nextflow run main.nf --mode <paired_end|quantseq> --indir <input data directory> -profile <nextflow profile(s)>

   Mandatory Arguments:         
    --mode                  Analysis mode: paired_end for conventional paired-end RNA-seq or quantseq for single-end QuantSeq
    --indir                 Path to directory containing input data 

   Input data:      [Will look for data in directory specified in --indir by default, one or more of following 
                    need to be specified if in a different directory, a subdirectory, or in case of error in 
                    finding the data (glob pattern mismatch)]

    --reads                 Read files. Paired-end glob for paired_end mode, e.g. "rawReads/*_{R1,R2}.fastq.gz";
                            single-end glob for quantseq mode, e.g. "rawReads/*.fastq.gz"
    --cdna                  Reference cDNA file. Required when Salmon quantification or strandedness checking is enabled
    --fasta                 Reference genome fasta file
    --gff                   Reference genome GFF file
   
   Required Analysis Arguments:
    --do_trim              Run trim_galore on reads: true or false
    --check_strandedness   Run Salmon automatic library-type inference and compare with --featurecounts_strand: true or false
    --strandedness_fail_on_mismatch
                            Fail before FeatureCounts if Salmon inference disagrees with --featurecounts_strand: true or false
    --do_salmon            Run Salmon quantification: true or false
    --do_deseq2            Run DESeq2 differential expression analysis: true or false
    --hisat_args           HISAT2 alignment arguments. Pass "" if no extra options are intended
    --featurecounts_strand FeatureCounts strandedness: 0 unstranded, 1 stranded, 2 reversely stranded
    --featurecounts_gene_type
                            Annotation feature type for gene counts, e.g. gene or exon
    --featurecounts_mrna_type
                            Annotation feature type for mRNA counts, e.g. mRNA or transcript
    --featurecounts_attribute
                            Annotation attribute used as the count ID, e.g. ID or gene_id
    --featurecounts_extra_args
                            Extra FeatureCounts arguments. Pass "" if no extra options are intended

   Conditional Arguments:
    --trim_args            Required when --do_trim true in paired_end mode
    --quant_trim_args      Required when --do_trim true in quantseq mode unless --trim_args is provided
    --strandedness_salmon_args
                            Required when --check_strandedness true. Extra Salmon args for inference, e.g. "--validateMappings"
    --hisat_index_args     Required when --hisatindex is not provided. Pass "" if no extra options are intended
    --salmon_index_args    Required when Salmon is enabled and --salmonindex is not provided. Pass "" if no extra options are intended
    --sal_quant_args       Required when Salmon is enabled. Include library type and other Salmon quant options
    --metadata             Required when --do_deseq2 true. Sample metadata CSV/TSV with a sample column
    --deseq_design         Required when --do_deseq2 true. DESeq2 design formula, e.g. "~ batch + treatment"
    --deseq_contrast       Required when --do_deseq2 true. DESeq2 contrast, e.g. "treatment,treated,control"

   Optional Arguments:    [default value]
    --threads               Number of threads [16]
    --outdir                Output directory name [RNAseq_Results]
    --salmonindex           Path to salmon index. Provide directory containing prebuilt salmon index files 
                            [If not provided, index is built by default]
    --hisatindex            Path to hisat index. Provide directory containing prebuilt Hisat2 index files 
                            [If not provided, Hisat will build an index by default] 
    
   Nextflow Arguments: (notice single "-" instead of double "--") 
    -profile                Nextflow profiles available: singularity, docker, slurm
    -resume                 Resume last run

    --help                  Print this help statement      
           """
     .stripIndent()
}

if(params.help){
  helpMsg()
  exit 0
}

include { rnaseq } from './Workflows/rnaseq.nf'

// Workflow
workflow {

 rnaseq()

}
