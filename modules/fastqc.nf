// Quality Control using fastqc

process fastqc {
   label 'fastqc'
   tag "FASTQC on $reads"
    
   input:
    path(reads)

   output:
    path("*_fastqc*")

    publishDir "${params.outdir}/01_Fastqc", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith('.html'))
       "html/${filename}"
    else if (filename.endsWith('.zip'))
       "zips/${filename}"
    else
       filename
   }

   script:
    """
    fastqc -t ${params.threads} -f fastq -q ${reads}
    """

}
