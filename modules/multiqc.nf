// Generate multiqc report

process multiqc {
   label 'multiqc'
    
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
