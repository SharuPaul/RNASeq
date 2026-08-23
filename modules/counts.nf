// FeatureCounts

process featureCounts_gene {
    tag "${read_bam.simpleName}"
    label 'counts_gene'
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'
 
    input: 
     tuple path(read_bam), path(genome_gff), val(strandedness_gate)

    output: 
     path("*")

    publishDir "${params.outdir}/05_Counts/Gene", mode: 'copy'

    script:
     def paired_flag = params.mode.toString().toLowerCase().replace('-', '_') == 'paired_end' ? '-p' : ''
     """
     featureCounts -T ${params.threads} ${params.featurecounts_extra_args} ${paired_flag} -s ${params.featurecounts_strand} -t ${params.featurecounts_gene_type} -g ${params.featurecounts_attribute} -a ${genome_gff} \
      -o ${read_bam.simpleName}_genecounts.txt \
      ${read_bam}
     """
}

process featureCounts_mRNA {
    tag "${read_bam.simpleName}"
    label 'counts_mRNA'
    container 'quay.io/biocontainers/subread'
 
    input:
     tuple path(read_bam), path(genome_gff), val(strandedness_gate)

    output:
     path("*")

    publishDir "${params.outdir}/05_Counts/mRNA", mode: 'copy'

    script:
     def paired_flag = params.mode.toString().toLowerCase().replace('-', '_') == 'paired_end' ? '-p' : ''
     """
     featureCounts -T ${params.threads} ${params.featurecounts_extra_args} ${paired_flag} -s ${params.featurecounts_strand} -t ${params.featurecounts_mrna_type} -g ${params.featurecounts_attribute} -a ${genome_gff} \
      -o ${read_bam.simpleName}_mRNAcounts.txt \
      ${read_bam}
     """
}

process featureCounts_geneMult {
    tag "${read_bam.simpleName}"
    label 'counts_Multimap'
    container 'quay.io/biocontainers/subread'
 
    input:
     tuple path(read_bam), path(genome_gff), val(strandedness_gate)

    output:
     path("*")

    publishDir "${params.outdir}/05_Counts/MultiMapping", mode: 'copy'

    script:
     def paired_flag = params.mode.toString().toLowerCase().replace('-', '_') == 'paired_end' ? '-p' : ''
     """
     featureCounts -T ${params.threads} ${params.featurecounts_extra_args} -M ${paired_flag} -s ${params.featurecounts_strand} -t ${params.featurecounts_gene_type} -g ${params.featurecounts_attribute} \
      -a ${genome_gff} -o ${read_bam.simpleName}_geneMultcounts.txt \
      ${read_bam}
     """
}
