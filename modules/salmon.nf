// Create binary index and quantify with Salmon

process salmon_index { 
    label 'salmon_index'
    
   input:
    path(cdna)

   output:
    path("*")

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

   script:
    """
    proc=\$(((`nproc`)))
    salmon index --threads \$proc -t $cdna -i index 
    """
}


process salmon_quant {
    label 'salmon_quantification'
    tag "$pair_id"
     
   input:
    tuple path(index), val(pair_id), path(reads)

   output:
    path("*")

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

   script:
    """
    salmon quant $params.sal_quant_args -i $index -1 ${reads[0]} -2 ${reads[1]} -o ${pair_id}_quant 
    """
}
