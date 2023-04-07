// Trimming 

process trim_galore {
    tag "$name"
    label "trim"
        
   input:
    tuple val(name), path(reads)
  
   output:
    tuple val(name), path ("*.fq"), emit: trim_reads
    path ("*"), emit: trim_fqc
   
   publishDir "${params.outdir}/02_Trim", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith('.fq'))
       "trim_reads/${filename}"
    else if (filename.endsWith('.zip'))
       "zips/${filename}"
    else if (filename.contains('report'))
       "reports/${filename}"
    else if (filename.endsWith('fastqc.html'))
       "fastqc/${filename}" 
    else
       filename
   }

   script:
    """
    proc=\$(((`nproc`)))
    if ((\$proc > 14)); then ncores=4
    elif ((\$proc > 11)); then ncores=3
    elif ((\$proc > 8)); then ncores=2
    else ncores=1
    fi
    trim_galore --cores \$ncores $params.trim_args --paired ${reads}
    """
}
