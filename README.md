# RNASeq
A pipeline for RNASeq analysis on paired-end reads implemented with NextFlow dsl2.


## Workflow
1. Fastqc - Quality Check
2. Trim_galore - Adapter trimming and fastqc
3. Salmon - Index building and quantification
4. Hisat2 - Index building and Alignment
5. Samtools - sam to bam conversion, generate stats report with flagstat
6. FeatureCounts - Count genes, mRNAs, and genes with multi-mapping reads
7. DESeq2 - Differential expression analysis from raw gene counts
8. Multiqc - Generate a multiqc report


## Requirements
1. Nextflow
2. Either Singularity or Docker to use containers. If not using containers, these software/modules are needed: fastqc, trimgalore, salmon, hisat2, samtools, subread, multiqc, and bioconductor-deseq2 when running differential expression.
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
   nextflow run main.nf --indir <input data directory> -profile <nextflow profile(s)>

   Mandatory Arguments:         
    --indir                 Path to directory containing input data 

   Input data:      [Will look for data in directory specified in --indir by default, one or more of following 
                    need to be specified if in a different directory, a subdirectory, or in case of error in 
                    finding the data (glob pattern mismatch)]

    --reads                 Paired-end reads (glob pattern, e.g. "rawReads/*_{R1,R2}.fastq.gz")
    --cdna                  Reference cDNA file
    --fasta                 Reference genome fasta file
    --gff                   Reference genome GFF file
   
   Optional Arguments:    [default value]
    --threads               Number of threads [16]
    --outdir                Output directory name [RNAseq_Results]
    --trim_args             Additional arguments for trim_galore ["--fastqc"]
    --salmonindex           Path to salmon index. Provide directory containing prebuilt salmon index files 
                            [If not provided, index is built by default]
    --sal_quant_args        Additional arguments for salmon quant ["--libType=A --validateMappings"]
    --hisatindex            Path to hisat index. Provide directory containing prebuilt Hisat2 index files 
                            [If not provided, Hisat will build an index by default] 
    --do_deseq2             Run DESeq2 differential expression analysis [false]
    --metadata              Sample metadata CSV/TSV for DESeq2. Must include a sample column matching sample names
    --deseq_design          DESeq2 design formula, e.g. "~ batch + treatment"
    --deseq_contrast        DESeq2 contrast, e.g. "treatment,treated,control"
    
   Nextflow Arguments: (notice single "-" instead of double "--") 
    -profile                Nextflow profiles available: singularity, docker, slurm
    -resume                 Resume last run

    --help                  Print this help statement   
```

Run the pipeline using this command:

```
nextflow run main.nf --indir <input data directory> -profile <nextflow profile(s)>
```

Prebuilt indexes for salmon and hisat can be supplied, and addtitional nextflow arguments can also be used. The program will look for input data in directory specified by `--indir` by default. If some data is in a different folder or a subfolder, and it cannot be located automatically, then you can specify that using the appropriate arguments (e.g. `--reads` or `--cdna` ).

## Differential expression
DESeq2 can be enabled after gene counting by providing sample metadata, a design formula, and a contrast:

```
nextflow run main.nf --indir <input data directory> -profile <nextflow profile(s)> \
  --do_deseq2 true \
  --metadata samples.csv \
  --deseq_design "~ batch + treatment" \
  --deseq_contrast "treatment,treated,control"
```

The metadata file may be CSV or TSV and must include a `sample` column. Additional columns can be used in the design formula as experimental conditions or covariates:

```
sample,treatment,batch,timepoint
sample1,control,batch1,day0
sample2,control,batch2,day0
sample3,treated,batch1,day0
sample4,treated,batch2,day0
```

Sample names must match the names produced by the paired-read input. DESeq2 results are written to `06_DifferentialExpression` inside the output directory.
