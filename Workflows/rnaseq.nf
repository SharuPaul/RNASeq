#! /usr/bin/env nextflow

include { fastqc } from '../modules/fastqc.nf'
include { salmon_index; salmon_quant } from '../modules/salmon.nf'
include { multiqc } from '../modules/multiqc.nf'
include { trim_galore } from '../modules/trim_galore.nf' 
include { hisat_index; hisat; samtools } from '../modules/hisat.nf'
include { featureCounts_gene; featureCounts_mRNA; featureCounts_geneMult } from '../modules/counts.nf'


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
    trim_galore(readpairs_ch)

  // Salmon quantification 
   if (params.salmonindex) {
     pre_salmon_ch = channel.fromPath(params.salmonindex, checkIfExists:true)
     pre_salmon_ch | combine(trim_galore.out.trim_reads) | salmon_quant
   } else {
     salmon_index(cdna_ch) | combine(trim_galore.out.trim_reads) | salmon_quant 
   }

  // Hisat2 Alignment 
   if (params.hisatindex) {
    pre_hisat_ch = channel.fromPath(params.hisatindex, checkIfExists:true)
    pre_hisat_ch | combine(trim_galore.out.trim_reads) | hisat
   } else {
     hisat_index(genome_ch)  
     hisat_index.out | combine(trim_galore.out.trim_reads) | hisat 
   }  
  
  // Sam to bam and stats
    samtools(hisat.out.read_sam)
 
  // Counts   
    samtools.out.read_bam | combine(gff_ch) | (featureCounts_gene & featureCounts_mRNA & featureCounts_geneMult)

  // Multiqc    
    salmon_quant.out | concat(fastqc.out) | concat(trim_galore.out.trim_fqc) | concat(featureCounts_gene.out) | concat(featureCounts_mRNA.out) | concat(featureCounts_geneMult.out) | collect | multiqc

}

workflow.onComplete { 
    println ( workflow.success ? "Workflow Completed Successfully!" : "Oops .. something went wrong" )
    
}
