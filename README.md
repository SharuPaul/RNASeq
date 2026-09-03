# RNASeq
A Nextflow DSL2 pipeline for conventional paired-end RNA-seq or single-end QuantSeq reads.


## Workflow
1. FastQC - Quality check
2. Trim Galore - Adapter trimming and FastQC
3. Salmon - Optional index building and quantification
4. HISAT2 - Index building and alignment
5. SAMtools - SAM to BAM conversion and flagstat report
6. featureCounts - Gene, mRNA, and multi-mapping counts with mode-aware paired-end and strandedness settings
7. DESeq2 - Optional differential expression analysis from raw gene counts
8. MultiQC - Generate a MultiQC report


## Requirements
1. Nextflow
2. Either Singularity or Docker to use containers. Without containers, these tools are needed: fastqc, trim-galore, salmon, hisat2, samtools, subread, multiqc, and the R packages DESeq2 and ggplot2 when running differential expression.
3. Git


## Usage
Clone the repo using this code:

```
git clone git@github.com:SharuPaul/RNASeq.git
```

And run this command to get help statement:

```
nextflow run main.nf --help
```

```
Usage:
   nextflow run main.nf --mode <paired_end|quantseq> --indir <input data directory> -profile <nextflow profile(s)>

   Mandatory Arguments:         
    --mode                  Analysis mode: paired_end for conventional paired-end RNA-seq or quantseq for single-end QuantSeq
    --indir                 Path to directory containing input data 

   Input data:      [Will look for data in directory specified in --indir by default, one or more of following 
                    need to be specified if in a different directory, a subdirectory, or in case of error in 
                    finding the data (glob pattern mismatch)]

    --reads                 Read files. Paired-end glob for paired_end mode, e.g. "rawReads/*_{R1,R2}.fastq.gz";
                            single-end glob for quantseq mode, e.g. "rawReads/*.fastq.gz"
    --cdna                  Reference cDNA file. Required when Salmon quantification or strandedness checking is enabled
    --fasta                 Reference genome fasta file
    --gff                   Reference genome GFF file
   
   Required Analysis Arguments:
    --do_trim              Run trim_galore on reads: true or false
    --check_strandedness   Run Salmon automatic library-type inference and compare with --featurecounts_strand: true or false
    --strandedness_fail_on_mismatch
                            Fail before FeatureCounts if Salmon inference disagrees with --featurecounts_strand: true or false
    --do_salmon            Run Salmon quantification: true or false
    --do_deseq2            Run DESeq2 differential expression analysis: true or false
    --hisat_args           HISAT2 alignment arguments.
    --featurecounts_strand FeatureCounts strandedness: 0 unstranded, 1 stranded, 2 reversely stranded
    --featurecounts_gene_type
                            Annotation feature type for gene counts, e.g. gene or exon
    --featurecounts_mrna_type
                            Annotation feature type for mRNA counts, e.g. mRNA or transcript
    --featurecounts_attribute
                            Annotation attribute used as the count ID, e.g. ID or gene_id
    --featurecounts_extra_args
                            Extra FeatureCounts arguments.

   Conditional Arguments:
    --trim_args            Required when --do_trim true in paired_end mode
    --quant_trim_args      Required when --do_trim true in quantseq mode unless --trim_args is provided
    --strandedness_salmon_args
                            Required when --check_strandedness true. Extra Salmon args for inference, e.g. "--validateMappings"
    --hisat_index_args     Required when --hisatindex is not provided.
    --salmon_index_args    Required when Salmon is enabled and --salmonindex is not provided.
    --sal_quant_args       Required when Salmon is enabled. Include library type and other Salmon quant options
    --metadata             Required when --do_deseq2 true. Sample metadata CSV/TSV with a sample column
    --deseq_design         Required when --do_deseq2 true. DESeq2 design formula, e.g. "~ batch + treatment"
    --deseq_contrast       Required when --do_deseq2 true unless --deseq_contrast_file is provided. One or more
                            DESeq2 contrasts separated by semicolons, e.g. "treatment,treated,control;treatment,dose2,control"
    --deseq_contrast_file  Required when --do_deseq2 true unless --deseq_contrast is provided. One contrast per line

   Optional Arguments:    [default value]
    --threads               Number of threads [16]
    --outdir                Output directory name [RNAseq_Results]
    --salmonindex           Path to salmon index. Provide directory containing prebuilt salmon index files 
                            [If not provided, index is built by default]
    --hisatindex            Path to hisat index. Provide directory containing prebuilt Hisat2 index files 
                            [If not provided, Hisat will build an index by default] 
    
   Nextflow Arguments: (notice single "-" instead of double "--") 
    -profile                Nextflow profiles available: singularity, docker, slurm
    -resume                 Resume last run

    --help                  Print this help statement   
```

Run the pipeline using this command:

```
nextflow run main.nf --mode paired_end --indir <input data directory> -profile <nextflow profile(s)> \
  --do_trim true \
  --trim_args "--fastqc" \
  --check_strandedness false \
  --strandedness_fail_on_mismatch true \
  --do_salmon true \
  --salmon_index_args "" \
  --sal_quant_args "--libType=A --validateMappings" \
  --hisat_index_args "" \
  --hisat_args "" \
  --featurecounts_strand 0 \
  --featurecounts_gene_type gene \
  --featurecounts_mrna_type mRNA \
  --featurecounts_attribute ID \
  --featurecounts_extra_args "" \
  --do_deseq2 false
```

Prebuilt indexes for Salmon and HISAT2 can be supplied, and additional Nextflow arguments can also be used. The program will look for input data in the directory specified by `--indir` by default. If some data is in a different folder or a subfolder, specify it with the appropriate argument, such as `--reads` or `--cdna`.

Tool arguments that affect biological interpretation are intentionally required.

## Output layout
- `01_Fastqc`: FastQC reports for input reads
- `02_Trim`: trimmed reads and Trim Galore/FastQC reports when `--do_trim true`
- `03_Salmon`: Salmon index, quantification, and strandedness-check files when enabled
- `04_Hisat2`: HISAT2 index, SAM/BAM files, alignment logs, and SAMtools flagstat reports
- `05_Counts`: gene, mRNA, and multi-mapping featureCounts results
- `06_DifferentialExpression`: DESeq2 count matrix and results when `--do_deseq2 true`
- `multiqc_report.html`, `timeline.html`, and `report.html`: run summaries

## Annotation style
Annotation style means the feature names and attributes used inside the genome annotation file. FeatureCounts uses `-t` to choose which feature rows to count and `-g` to choose which attribute groups those rows into a count ID.

GFF3-style annotations often use feature types like `gene` and `mRNA` with attributes like `ID` or `Parent`. GTF-style annotations often use feature types like `exon` with attributes like `gene_id` or `transcript_id`.

For example:

```
--featurecounts_gene_type gene --featurecounts_attribute ID
```

may fit a GFF3 file, while:

```
--featurecounts_gene_type exon --featurecounts_attribute gene_id
```

is a common GTF-style gene-count setup.

## Modes
Use `--mode paired_end` for conventional paired-end RNA-seq. This uses paired-end trimming, HISAT2 alignment, and FeatureCounts counting when those steps are explicitly enabled.

Use `--mode quantseq` for single-end QuantSeq data. This uses single-end trimming, HISAT2 `-U` alignment, and featureCounts without paired-end counting. If `--reads` is not set, the default QuantSeq read glob is:

```
<indir>/*.{fq,fastq,fq.gz,fastq.gz}
```

QuantSeq sample names come from the FASTQ basename after removing `.fq`, `.fastq`, `.fq.gz`, or `.fastq.gz`.

Use `--quant_trim_args` for QuantSeq trimming. If it is not set, `--trim_args` is used. Salmon quantification can run in QuantSeq mode with single-end reads, but it is usually disabled by setting `--do_salmon false`.

For QuantSeq, set `--featurecounts_strand` to match the library strandedness:

```
0 = unstranded
1 = stranded
2 = reversely stranded
```

You can ask Salmon to check strandedness by setting `--check_strandedness true`. This requires `--cdna` or a matching cDNA file under `--indir`. The pipeline runs Salmon with `-l A`, writes strandedness reports under `03_Salmon/Strandedness`, and compares Salmon's inferred library type with `--featurecounts_strand`.

The Salmon-to-featureCounts mapping is:

```
U, IU, OU, MU       -> featureCounts -s 0
SF, ISF, OSF, MSF  -> featureCounts -s 1
SR, ISR, OSR, MSR  -> featureCounts -s 2
```

Example QuantSeq run:

```
nextflow run main.nf --mode quantseq --indir <input data directory> -profile <nextflow profile(s)> \
  --do_trim true \
  --quant_trim_args "--fastqc" \
  --check_strandedness true \
  --strandedness_salmon_args "--validateMappings" \
  --strandedness_fail_on_mismatch true \
  --do_salmon false \
  --featurecounts_strand 2 \
  --featurecounts_gene_type gene \
  --featurecounts_mrna_type mRNA \
  --featurecounts_attribute ID \
  --do_deseq2 false
```

## Differential expression
DESeq2 can be enabled after gene counting by providing sample metadata, a design formula, and one or more contrasts. It uses the gene-level featureCounts files ending in `_genecounts.txt`.

```
nextflow run main.nf --mode paired_end --indir <input data directory> -profile <nextflow profile(s)> \
  --do_trim true \
  --trim_args "--fastqc" \
  --check_strandedness false \
  --strandedness_fail_on_mismatch true \
  --do_salmon true \
  --sal_quant_args "--libType=A --validateMappings" \
  --featurecounts_strand 0 \
  --featurecounts_gene_type gene \
  --featurecounts_mrna_type mRNA \
  --featurecounts_attribute ID \
  --do_deseq2 true \
  --metadata samples.csv \
  --deseq_design "~ batch + treatment" \
  --deseq_contrast "treatment,treated,control;treatment,dose2,control"
```

The metadata file may be CSV or TSV and must include a `sample` column. Additional columns can be used in the design formula as experimental conditions or covariates:

```
sample,treatment,batch,timepoint
sample1,control,batch1,day0
sample2,control,batch2,day0
sample3,treated,batch1,day0
sample4,treated,batch2,day0
```

Contrasts can also be provided in a file with one contrast per line. Blank lines and lines beginning with `#` are ignored:

```
# variable,numerator,denominator
treatment,treated,control
treatment,dose2,control
```

Then use:

```
--deseq_contrast_file contrasts.txt
```

Sample names must match the names produced from the read input. Metadata samples must be unique. Use either `--deseq_contrast` or `--deseq_contrast_file`, not both.

The merged raw count matrix is written to `06_DifferentialExpression/gene_count_matrix.tsv`. DESeq2 result files are written in a design-specific subfolder, such as `06_DifferentialExpression/batch_treatment`.

DESeq2 outputs:

```
gene_count_matrix.tsv
<design>/deseq2_results_<variable>_<numerator>_<denominator>.tsv
<design>/deseq2_summary_<variable>_<numerator>_<denominator>.tsv
<design>/deseq2_normalized_counts.tsv
<design>/deseq2_vst_counts.tsv
<design>/deseq2_pca.pdf
<design>/deseq2_pca_ellipse.pdf
<design>/deseq2_pca_centroids.tsv
<design>/deseq2_pca_centroid_distances.tsv
<design>/deseq2_ma_plot_<variable>_<numerator>_<denominator>.pdf
<design>/deseq2_volcano_<variable>_<numerator>_<denominator>.pdf
<design>/deseq2_volcano_<variable>_<numerator>_<denominator>.png
<design>/deseq2_session_info.txt
```

## Citation
If you use this pipeline, please cite:

Sharu Paul. RNASeq v1.0.0. GitHub: https://github.com/SharuPaul/RNASeq

Citation metadata is also available in `CITATION.cff`.
