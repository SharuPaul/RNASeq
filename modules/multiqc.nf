// Generate multiqc report

process multiqc {
   label 'multiqc'
   container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    
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
