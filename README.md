# RNASeq
A pipeline for RNASeq analysis on paired-end reads implemented with NextFlow dsl2.


## Workflow
1. Fastqc - Quality Check
2. Trim_galore - Adapter trimming and fastqc - trimmed reads are used for the rest of the workflow
3. Salmon - Index building and quantification
4. Hisat2 - Index building and Alignment
5. Samtools - sam to bam conversion, generate stats report with flagstat
6. FeatureCounts - Count genes, mRNAs, and genes with multi-mapping reads
7. Multiqc - Generate a multiqc report


## Requirements
1. Nextflow
2. Either Singularity or Docker to use containers. If not using containers, these software/modules are needed: fastqc, trimgalore, salmon, hisat2, samtools, subread, and multiqc.
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

   Input data:   [Will look for data in directory specified in --indir by default, one or more of following 
               need to be specified by glob pattern (e.g. --reads "rawReads/*_{R1,R2}.fastq.gz") if in a 
               different directory, a subdirectory, or in case of error in finding the data (glob pattern mismatch)]

    --reads                 Paired-end reads
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
    
   Nextflow Arguments: (notice single "-" instead of double "--") 
    -profile                Nextflow profiles available: singularity, docker, slurm
    -resume                 Resume last run

    --help                  Print this help statement 
```

Run the pipeline using this command:

```
nextflow run main.nf --indir <input data directory> -profile <nextflow profile(s)>
```

You can also supply prebuilt indexes for salmon and hisat, and use any nextflow arguments. The program will look for input data in directory specified by --indir by default. If some data is in a different folder or a subfolder, and it cannot be located automatically, then you can specify that using the appropriate arguments (reads, cdna, fasta, gff).
