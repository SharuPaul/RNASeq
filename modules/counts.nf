// FeatureCounts

process featureCounts_gene {
    tag "${read_bam.simpleName}"
    label 'counts_gene'
 
    input: 
     tuple path(read_bam), path(genome_gff)

    output: 
     path("*")

    publishDir "${params.outdir}/05_Counts/Gene", mode: 'copy'

    script:
     """
     featureCounts -T ${params.threads} -t gene -g ID -p -a ${genome_gff} \
      -o ${read_bam.simpleName}_genecounts.txt \
      ${read_bam}
     """
}

process featureCounts_mRNA {
    tag "${read_bam.simpleName}"
    label 'counts_mRNA'
 
    input:
     tuple path(read_bam), path(genome_gff)

    output:
     path("*")

    publishDir "${params.outdir}/05_Counts/mRNA", mode: 'copy'

    script:
     """
     featureCounts -T ${params.threads} -t mRNA -g ID -p -a ${genome_gff} \
      -o ${read_bam.simpleName}_mRNAcounts.txt \
      ${read_bam}
     """
}

process featureCounts_geneMult {
    tag "${read_bam.simpleName}"
    label 'counts_Multimatch'
 
    input:
     tuple path(read_bam), path(genome_gff)

    output:
     path("*")

    publishDir "${params.outdir}/05_Counts/MultiMapping", mode: 'copy'

    script:
     """
     featureCounts -T ${params.threads} -M -p -t gene -g ID \
      -a ${genome_gff} -o ${read_bam.simpleName}_geneMultcounts.txt \
      ${read_bam}
     """
}
