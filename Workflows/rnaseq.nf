#! /usr/bin/env nextflow

include { fastqc } from '../modules/fastqc.nf'
include { salmon_index; salmon_quant } from '../modules/salmon.nf'
include { multiqc } from '../modules/multiqc.nf'
include { trim_galore } from '../modules/trim_galore.nf' 
include { hisat_index; hisat; samtools } from '../modules/hisat.nf'
include { featureCounts_gene; featureCounts_mRNA; featureCounts_geneMult } from '../modules/counts.nf'
include { deseq2_from_featurecounts } from '../modules/deseq2.nf'


// Complete RNASeq Workflow

workflow rnaseq {

if (!params.reads) { reads = "${params.indir}/*_{1,2}.{fq,fastq,fq.gz,fastq.gz}" 
  } else { reads = "${params.reads}" } 
if (!params.cdna) { cdna = "${params.indir}/*rna.{fna,fna.gz}" 
  } else { cdna = "${params.cdna}" }
if (!params.fasta) { fasta = "${params.indir}/*genomic.{fa,fasta,fna,fna.gz}" 
  } else { fasta = "${params.fasta}" }
if (!params.gff) { gff = "${params.indir}/*.{gff,gff.gz}"
  } else { gff = "${params.gff}" }

// Channels  
    reads_ch = channel.fromPath(reads, checkIfExists:true)
    readpairs_ch = channel.fromFilePairs(reads, checkIfExists:true)
    cdna_ch = channel.fromPath(cdna, checkIfExists:true)
    genome_ch = channel.fromPath(fasta, checkIfExists:true)
    gff_ch = channel.fromPath(gff, checkIfExists:true) 
  
// Processes

  // Fastqc
    reads_ch.collate(4) | fastqc

  // Trimming
   if (params.do_trim) {
     trim_galore(readpairs_ch)
     workflow_reads_ch = trim_galore.out.trim_reads
   } else {
     workflow_reads_ch = readpairs_ch
   }

  // Salmon quantification 
   if (params.salmonindex) {
     pre_salmon_ch = channel.fromPath(params.salmonindex, checkIfExists:true)
     pre_salmon_ch | combine(workflow_reads_ch) | salmon_quant
   } else {
     salmon_index(cdna_ch) | combine(workflow_reads_ch) | salmon_quant 
   }

  // Hisat2 Alignment 
   if (params.hisatindex) {
    pre_hisat_ch = channel.fromPath(params.hisatindex, checkIfExists:true)
    pre_hisat_ch | combine(workflow_reads_ch) | hisat
   } else {
     hisat_index(genome_ch)  
     hisat_index.out | combine(workflow_reads_ch) | hisat 
   }  
  
  // Sam to bam and stats
    samtools(hisat.out.read_sam)
 
  // Counts   
    samtools.out.read_bam | combine(gff_ch) | (featureCounts_gene & featureCounts_mRNA & featureCounts_geneMult)

  // Differential expression
   if (params.do_deseq2) {
     if (!params.metadata) {
       error "DESeq2 requires --metadata when --do_deseq2 is true"
     }
     if (!params.deseq_design) {
       error "DESeq2 requires --deseq_design when --do_deseq2 is true"
     }
     if (!params.deseq_contrast) {
       error "DESeq2 requires --deseq_contrast when --do_deseq2 is true"
     }

     metadata_ch = channel.fromPath(params.metadata, checkIfExists:true)
     gene_counts_ch = featureCounts_gene.out.flatten().filter { count_file -> count_file.name.endsWith('_genecounts.txt') }
     deseq2_from_featurecounts(gene_counts_ch.collect(), metadata_ch, params.deseq_design, params.deseq_contrast)
   }

  // Multiqc    
   if (params.do_trim) {
     salmon_quant.out | concat(fastqc.out) | concat(trim_galore.out.trim_fqc) | concat(featureCounts_gene.out) | concat(featureCounts_mRNA.out) | concat(featureCounts_geneMult.out) | collect | multiqc
   } else {
     salmon_quant.out | concat(fastqc.out) | concat(featureCounts_gene.out) | concat(featureCounts_mRNA.out) | concat(featureCounts_geneMult.out) | collect | multiqc
   }

}

workflow.onComplete { 
    println ( workflow.success ? "Workflow Completed Successfully!" : "Oops .. something went wrong" )
    
}
