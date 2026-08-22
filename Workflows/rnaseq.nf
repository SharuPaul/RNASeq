#! /usr/bin/env nextflow

include { fastqc } from '../modules/fastqc.nf'
include { salmon_index; salmon_quant; salmon_infer_strandedness } from '../modules/salmon.nf'
include { multiqc } from '../modules/multiqc.nf'
include { trim_galore } from '../modules/trim_galore.nf' 
include { hisat_index; hisat; samtools } from '../modules/hisat.nf'
include { featureCounts_gene; featureCounts_mRNA; featureCounts_geneMult } from '../modules/counts.nf'
include { deseq2_from_featurecounts } from '../modules/deseq2.nf'


// Complete RNASeq Workflow

workflow rnaseq {

mode = params.mode ? params.mode.toString().toLowerCase().replace('-', '_') : null
if (!mode) {
  error "Missing required argument --mode. Use --mode paired_end or --mode quantseq"
}
if (!['paired_end', 'quantseq'].contains(mode)) {
  error "Invalid --mode '${params.mode}'. Use --mode paired_end or --mode quantseq"
}

is_paired_end = mode == 'paired_end'
require_param = { value, name, help_text = null ->
  if (value == null) {
    error "Missing required argument --${name}" + (help_text ? ". ${help_text}" : "")
  }
}

to_boolean = { value, name ->
  require_param(value, name, "Use true or false")
  if (value instanceof Boolean) {
    return value
  }
  def normalized_value = value.toString().toLowerCase()
  if (!['true', 'false', '1', '0', 'yes', 'no', 'y', 'n'].contains(normalized_value)) {
    error "Invalid --${name} '${value}'. Use true or false"
  }
  return normalized_value in ['true', '1', 'yes', 'y']
}

do_trim = to_boolean(params.do_trim, 'do_trim')
do_salmon = to_boolean(params.do_salmon, 'do_salmon')
do_deseq2 = to_boolean(params.do_deseq2, 'do_deseq2')
check_strandedness = to_boolean(params.check_strandedness, 'check_strandedness')
strandedness_fail_on_mismatch = to_boolean(params.strandedness_fail_on_mismatch, 'strandedness_fail_on_mismatch')

if (!params.indir) {
  error "Missing required argument --indir"
}

if (do_trim) {
  if (is_paired_end && params.trim_args == null) {
    error "Missing required argument --trim_args when --do_trim true. Set this from the read/QC requirements"
  }
  if (!is_paired_end && params.trim_args == null && params.quant_trim_args == null) {
    error "Missing required trimming arguments when --do_trim true in quantseq mode. Provide --quant_trim_args or --trim_args"
  }
}

if (!params.hisatindex) {
  require_param(params.hisat_index_args, 'hisat_index_args', 'Pass an empty quoted string if no HISAT2 index options are intended')
}
require_param(params.hisat_args, 'hisat_args', 'Pass an empty quoted string if no HISAT2 alignment options are intended')

if (do_salmon) {
  if (!params.salmonindex) {
    require_param(params.salmon_index_args, 'salmon_index_args', 'Pass an empty quoted string if no Salmon index options are intended')
  }
  require_param(params.sal_quant_args, 'sal_quant_args', 'Set library type and other Salmon quant options explicitly')
}

if (check_strandedness) {
  if (!params.salmonindex) {
    require_param(params.salmon_index_args, 'salmon_index_args', 'Pass an empty quoted string if no Salmon index options are intended')
  }
  require_param(params.strandedness_salmon_args, 'strandedness_salmon_args', 'Set Salmon options for library-type inference explicitly')
}

require_param(params.featurecounts_strand, 'featurecounts_strand', 'Use 0 unstranded, 1 stranded, or 2 reversely stranded')
require_param(params.featurecounts_gene_type, 'featurecounts_gene_type', 'Set the annotation feature type for gene counts')
require_param(params.featurecounts_mrna_type, 'featurecounts_mrna_type', 'Set the annotation feature type for mRNA counts')
require_param(params.featurecounts_attribute, 'featurecounts_attribute', 'Set the annotation attribute used as the count ID')
require_param(params.featurecounts_extra_args, 'featurecounts_extra_args', 'Pass an empty quoted string if no extra FeatureCounts options are intended')

featurecounts_strand = params.featurecounts_strand.toString()

if (!['0', '1', '2'].contains(featurecounts_strand)) {
  error "Invalid --featurecounts_strand '${params.featurecounts_strand}'. Use 0, 1, or 2"
}

if (!params.reads) {
  reads = is_paired_end ? "${params.indir}/*_{1,2}.{fq,fastq,fq.gz,fastq.gz}" : "${params.indir}/*.{fq,fastq,fq.gz,fastq.gz}"
  } else { reads = "${params.reads}" }
if (!params.fasta) { fasta = "${params.indir}/*genomic.{fa,fasta,fna,fna.gz}"
  } else { fasta = "${params.fasta}" }
if (!params.gff) { gff = "${params.indir}/*.{gff,gff.gz}"
  } else { gff = "${params.gff}" }

// Channels  
    reads_ch = channel.fromPath(reads, checkIfExists:true)
    if (is_paired_end) {
      workflow_input_reads_ch = channel.fromFilePairs(reads, checkIfExists:true)
    } else {
      workflow_input_reads_ch = channel.fromPath(reads, checkIfExists:true).map { read -> tuple(read.name.replaceFirst(/(\.fq|\.fastq)(\.gz)?$/, ''), read) }
    }
    genome_ch = channel.fromPath(fasta, checkIfExists:true)
    gff_ch = channel.fromPath(gff, checkIfExists:true) 
  
// Processes

  // Fastqc
    reads_ch.collate(4) | fastqc

  // Trimming
   if (do_trim) {
     trim_galore(workflow_input_reads_ch)
     workflow_reads_ch = trim_galore.out.trim_reads
   } else {
     workflow_reads_ch = workflow_input_reads_ch
   }

  // Salmon index for quantification and strandedness inference
   if (do_salmon || check_strandedness) {
     if (!params.cdna) { cdna = "${params.indir}/*rna.{fna,fna.gz}"
       } else { cdna = "${params.cdna}" }
     cdna_ch = channel.fromPath(cdna, checkIfExists:true)

     if (params.salmonindex) {
       salmon_index_ch = channel.fromPath(params.salmonindex, checkIfExists:true)
     } else {
       salmon_index(cdna_ch)
       salmon_index_ch = salmon_index.out
     }
   }

  // Salmon strandedness inference
   if (check_strandedness) {
     salmon_index_ch | combine(workflow_reads_ch) | salmon_infer_strandedness
   }

  // Salmon quantification 
   if (do_salmon) {
     salmon_index_ch | combine(workflow_reads_ch) | salmon_quant
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
   counts_input_ch = samtools.out.read_bam.combine(gff_ch)
   if (check_strandedness) {
     strandedness_gate_ch = salmon_infer_strandedness.out.report.collect()
     counts_input_ch = counts_input_ch.combine(strandedness_gate_ch).map { read_bam, genome_gff, strandedness_gate -> tuple(read_bam, genome_gff, 'checked') }
   } else {
     counts_input_ch = counts_input_ch.map { read_bam, genome_gff -> tuple(read_bam, genome_gff, 'not_checked') }
   }
   counts_input_ch | (featureCounts_gene & featureCounts_mRNA & featureCounts_geneMult)

  // Differential expression
   if (do_deseq2) {
     if (!params.metadata) {
       error "DESeq2 requires --metadata when --do_deseq2 is true"
     }
     if (!params.deseq_design) {
       error "DESeq2 requires --deseq_design when --do_deseq2 is true"
     }
     has_inline_contrast = params.deseq_contrast != null && params.deseq_contrast.toString().trim()
     has_contrast_file = params.deseq_contrast_file != null && params.deseq_contrast_file.toString().trim()

     if (has_inline_contrast && has_contrast_file) {
       error "DESeq2 accepts either --deseq_contrast or --deseq_contrast_file, not both"
     }
     if (!has_inline_contrast && !has_contrast_file) {
       error "DESeq2 requires --deseq_contrast or --deseq_contrast_file when --do_deseq2 is true"
     }

     deseq_contrasts = params.deseq_contrast
     if (has_contrast_file) {
       contrast_file = file(params.deseq_contrast_file)
       if (!contrast_file.exists()) {
         error "DESeq2 contrast file not found: ${params.deseq_contrast_file}"
       }

       raw_contrast_lines = contrast_file.text.readLines()
       contrast_lines = raw_contrast_lines
         .withIndex()
         .collect { line, index -> [number: index + 1, text: line.trim()] }
         .findAll { line -> line.text && !line.text.startsWith('#') }

       if (contrast_lines.isEmpty()) {
         error "DESeq2 contrast file has no contrast lines: ${params.deseq_contrast_file}"
       }

       contrast_lines.each { line ->
         contrast_parts = line.text.split(',').collect { part -> part.trim() }
         if (contrast_parts.size() != 3 || contrast_parts.any { part -> !part }) {
           error "Invalid DESeq2 contrast at ${params.deseq_contrast_file}:${line.number}. Use variable,numerator,denominator"
         }
       }

       deseq_contrasts = contrast_lines.collect { line -> line.text }.join(';')
     }

     metadata_ch = channel.fromPath(params.metadata, checkIfExists:true)
     gene_counts_ch = featureCounts_gene.out.flatten().filter { count_file -> count_file.name.endsWith('_genecounts.txt') }
     deseq2_from_featurecounts(gene_counts_ch.collect(), metadata_ch, params.deseq_design, deseq_contrasts)
   }

  // Multiqc    
   multiqc_ch = fastqc.out
   if (do_salmon) {
     multiqc_ch = multiqc_ch.concat(salmon_quant.out)
   }
   if (check_strandedness) {
     multiqc_ch = multiqc_ch.concat(salmon_infer_strandedness.out.report)
   }
   if (do_trim) {
     multiqc_ch = multiqc_ch.concat(trim_galore.out.trim_fqc)
   }

   multiqc_ch.concat(featureCounts_gene.out).concat(featureCounts_mRNA.out).concat(featureCounts_geneMult.out).collect() | multiqc

}

workflow.onComplete { 
    println ( workflow.success ? "Workflow Completed Successfully!" : "Oops .. something went wrong" )
    
}
