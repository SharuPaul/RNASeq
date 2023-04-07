// Hisat2 index and alignment, Samtools sam to bam, stats

process hisat_index {
    label "Hisat_index"
    container = 'nfcore/rnaseq'

   input:
    path(genome)

   output:
    path("hisat2_index"), emit: index
   
   publishDir "${params.outdir}/04_Hisat2", mode: 'copy'

   script:
    """
    mkdir hisat2_index
    proc=\$(((`nproc`)))
    hisat2-build -p \$proc ${genome} hisat2_index/${genome.simpleName}
    """

}

process hisat {
   label "Hisat2_align"
   
   input:
    tuple path(index), val(readname), path(read_pairs)

   output:
    path("*.sam"), emit: read_sam
    path("*.log"), emit: hisat_logs

    publishDir "${params.outdir}/04_Hisat2", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('.sam'))
         "sam_files/${filename}"
      else if (filename.endsWith('.log'))
         "hisat_logs/${filename}"
      else
         filename
     }

   script:
    """
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
    proc=\$(((`nproc`)))
    hisat2 -p \$proc -x \$INDEX -1 ${read_pairs[0]} -2 ${read_pairs[1]} -S ${readname}.sam &> ${readname}.log
    """
 
}


process samtools {
   label "Samtools"
   
   input:
    path(samfile)

   output:
    path("*_sorted.bam"), emit: read_bam
    path("*_stats.txt"), emit: stats

   publishDir "${params.outdir}/04_Hisat2", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('.bam'))
         "bam_files/${filename}"
      else if (filename.endsWith('_stats.txt'))
         "Samtools_flagstats/${filename}"
      else
         filename
     }

   script:
    """
    proc=\$(((`nproc`)))
    samtools view --threads \$proc -uS ${samfile} | samtools sort -o ${samfile.simpleName}_sorted.bam --threads \$proc -
    samtools flagstat ${samfile.simpleName}_sorted.bam > ${samfile.simpleName}_stats.txt
    """

}
