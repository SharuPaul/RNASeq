// Trimming 

process trim_galore {
    tag "$name"
    label "trim"
    container 'quay.io/biocontainers/trim-galore'
        
   input:
    tuple val(name), path(reads)
  
   output:
    tuple val(name), path ("trimmed_reads/*"), emit: trim_reads
    path ("trim_qc/*"), emit: trim_fqc
   
   publishDir "${params.outdir}/02_Trim", mode: 'copy', saveAs: { filename ->
    def basename = filename.tokenize('/').last()
    if (basename.endsWith('.fq') || basename.endsWith('.fq.gz') || basename.endsWith('.fastq') || basename.endsWith('.fastq.gz'))
       "trim_reads/${basename}"
    else if (basename.endsWith('.zip'))
       "zips/${basename}"
    else if (basename.contains('report'))
       "reports/${basename}"
    else if (basename.endsWith('fastqc.html'))
       "fastqc/${basename}"
    else
       basename
   }

   script:
    def mode = params.mode.toString().toLowerCase().replace('-', '_')
    def trim_mode = mode == 'paired_end' ? '--paired' : ''
    def trim_args = mode == 'quantseq' && params.quant_trim_args != null ? params.quant_trim_args : params.trim_args
    def read_args = reads instanceof List ? reads.join(' ') : reads
    """
    mkdir trimmed_reads trim_qc
    proc=\$(((`nproc`)))
    if ((\$proc > 14)); then ncores=4
    elif ((\$proc > 11)); then ncores=3
    elif ((\$proc > 8)); then ncores=2
    else ncores=1
    fi
    trim_galore --cores \$ncores ${trim_args} ${trim_mode} ${read_args}
    find . -maxdepth 1 -type f \\( -name '*_val_*.fq*' -o -name '*_trimmed.fq*' \\) -exec mv {} trimmed_reads/ \\;
    find . -maxdepth 1 -type f \\( -name '*trimming_report.txt' -o -name '*_fastqc.zip' -o -name '*_fastqc.html' \\) -exec mv {} trim_qc/ \\;
    """
}
