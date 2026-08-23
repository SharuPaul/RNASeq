// Create binary index and quantify with Salmon

process salmon_index { 
    label 'salmon_index'
    container 'quay.io/biocontainers/salmon'
    
   input:
    path(cdna)

   output:
    path("index")

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

   script:
    """
    proc=\$(((`nproc`)))
    salmon index --threads \$proc ${params.salmon_index_args} -t $cdna -i index
    """
}


process salmon_quant {
    label 'salmon_quantification'
    container 'quay.io/biocontainers/salmon'
    tag "$pair_id"
     
   input:
    tuple path(index), val(pair_id), path(reads)

   output:
    path("*")

    publishDir "${params.outdir}/03_Salmon", mode: 'copy'

   script:
    def mode = params.mode.toString().toLowerCase().replace('-', '_')
    def single_read = reads instanceof List ? reads[0] : reads
    def read_args = mode == 'paired_end' ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${single_read}"
    """
    salmon quant ${params.sal_quant_args} -i $index ${read_args} -o ${pair_id}_quant
    """
}

process salmon_infer_strandedness {
    label 'salmon_strandedness'
    container 'quay.io/biocontainers/salmon'
    tag "$sample_id"

   input:
    tuple path(index), val(sample_id), path(reads)

   output:
    path("${sample_id}_salmon_strandedness.tsv"), emit: report
    path("${sample_id}_salmon_lib_format_counts.json"), emit: lib_format_counts
    path("${sample_id}_salmon_meta_info.json"), emit: meta_info

    publishDir "${params.outdir}/03_Salmon/Strandedness", mode: 'copy'

   script:
    def mode = params.mode.toString().toLowerCase().replace('-', '_')
    def single_read = reads instanceof List ? reads[0] : reads
    def read_args = mode == 'paired_end' ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${single_read}"
    def fail_on_mismatch = params.strandedness_fail_on_mismatch.toString().toLowerCase() in ['true', '1', 'yes', 'y'] ? 'true' : 'false'
    """
    proc=\$(((`nproc`)))
    salmon quant -p \$proc -i $index -l A ${params.strandedness_salmon_args} ${read_args} -o ${sample_id}_salmon_auto

    meta_info="${sample_id}_salmon_auto/aux_info/meta_info.json"
    lib_counts="${sample_id}_salmon_auto/lib_format_counts.json"
    if [ ! -f "\$lib_counts" ]; then
      lib_counts="${sample_id}_salmon_auto/aux_info/lib_format_counts.json"
    fi

    cp "\$meta_info" "${sample_id}_salmon_meta_info.json"
    if [ -f "\$lib_counts" ]; then
      cp "\$lib_counts" "${sample_id}_salmon_lib_format_counts.json"
    else
      touch "${sample_id}_salmon_lib_format_counts.json"
    fi

    inferred=\$(sed -n 's/.*"expected_format"[[:space:]]*:[[:space:]]*"\\([^"]*\\)".*/\\1/p' "\$meta_info" | head -n 1)
    if [ -z "\$inferred" ] && [ -f "\$lib_counts" ]; then
      inferred=\$(sed -n 's/.*"expected_format"[[:space:]]*:[[:space:]]*"\\([^"]*\\)".*/\\1/p' "\$lib_counts" | head -n 1)
    fi

    case "\$inferred" in
      U|IU|OU|MU) inferred_featurecounts_strand="0" ;;
      SF|ISF|OSF|MSF) inferred_featurecounts_strand="1" ;;
      SR|ISR|OSR|MSR) inferred_featurecounts_strand="2" ;;
      *) inferred_featurecounts_strand="unknown" ;;
    esac

    if [ "\$inferred_featurecounts_strand" = "${params.featurecounts_strand}" ]; then
      status="PASS"
    else
      status="FAIL"
    fi

    {
      printf 'sample\\tsalmon_expected_format\\tinferred_featurecounts_strand\\tconfigured_featurecounts_strand\\tstatus\\n'
      printf '${sample_id}\\t%s\\t%s\\t${params.featurecounts_strand}\\t%s\\n' "\$inferred" "\$inferred_featurecounts_strand" "\$status"
    } > "${sample_id}_salmon_strandedness.tsv"

    if [ "${fail_on_mismatch}" = "true" ] && [ "\$status" = "FAIL" ]; then
      echo "Salmon inferred FeatureCounts strandedness \$inferred_featurecounts_strand from library type \$inferred, but --featurecounts_strand is ${params.featurecounts_strand}" >&2
      exit 1
    fi
    """
}
