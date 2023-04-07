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
      nextflow run main.nf --indir <input data directory> -profile <nextflow profile(s)>

   Mandatory Arguments:         
    --indir                 Path to directory containing input data 

   Input data:   [Will look for data in directory specified in --indir by default, one or more of following 
               need to be specified by glob pattern (e.g. --reads "rawReads/*_{R1,R2}.fastq.gz") if in a 
               different folder, a subfolder, or in case of error in finding the data (glob pattern mismatch)]

    --reads                 Paired-end reads
    --cdna                  Reference cDNA file
    --fasta                 Reference genome fasta file
    --gff                   Reference genome GFF file
   
   Optional Arguments:    [default value]
    --threads               Number of threads [36]
    --outdir                Output directory name [RNAseq_Results]
    --trim_args             Additional arguments for trim_galore ["--fastqc"]
    --salmonindex           Path to salmon index. Provide directory containing prebuilt salmon index files 
                            [If not provided, index is built by default]
    --sal_quant_args        Additional arguments for salmon quant ["--libType=A --validateMappings"]
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
