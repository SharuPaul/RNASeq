// Generate multiqc report

process multiqc {
   label 'multiqc'
   container 'quay.io/biocontainers/multiqc'
    
   input:
    path('*')

   output:
    path('multiqc_report.html')

    publishDir "${params.outdir}", mode: 'copy'

   script:
    """
    multiqc .
    """

}
