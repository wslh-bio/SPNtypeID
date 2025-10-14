# SPNtypeID

SPNtypeID is a [Nextflow](https://www.nextflow.io/) pipeline used for genome assembly and serotyping of *Streptococcus pneumoniae*.

![SPNtypeID](https://github.com/wslh-bio/SPNtypeID/actions/workflows/workflow_test.yml/badge.svg)
![GPL-3.0](https://img.shields.io/github/license/wslh-bio/SPNtypeID)
![Static Badge](https://img.shields.io/badge/release-v1.9.0-%2300FFFF)

### Table of Contents:
[Usage](#using-the-workflow)  
[Input](#input)  
[Parameters](#parameters)  
[Workflow outline](#workflow-outline)  
[Read trimming and quality assessment](#read-trimming-and-quality-assessment)  
[Genome assembly](#genome-assembly)  
[Assembly quality assessment](#genome-assembly-quality-assessment)  
[Genome coverage](#genome-coverage-quality-assessment)  
[Contamination detection](#contamination-detection)  
[Serotyping](#serotyping)                                                                                                                                  
[Output](#output-files)  
[Results file explanation](#results-file-explanation)             
[Credits](#credits)  
[Contributions and Support](#contributions-and-support)  
[Citations](#citations)  

### Using the workflow
The pipeline is designed to start from raw, paired-end Illumina reads. Start the pipeline using:
```
nextflow run SPNtypeID/main.nf --input [path-to-samplesheet] --outdir [path-to-outdir] --ntc_regex [what-expression-to-look-for-in-ntcs] --runname [what-to-call-run] -profile [docker,singularity,aws]
```

or from github using:
```
nextflow run wslh-bio/SPNtypeID -r [version] --input [path-to-samplesheet] --outdir [path-to-outdir] -ntc_regex [what-expression-to-look-for-in-ntcs] --runname [what-to-call-run] -profile [docker,singularity,aws]
```

You can also test the pipeline with example data using `-profile test` or `-profile test_full`:
```
nextflow run SPNtypeID/main.nf --outdir [path-to-outdir] -profile test[_full],[docker,singularity]
```
**As of version 1.8.0, the workflow validation step is not functional and the pipeline will fail if run with `-profile docker,test_full`.**


### Input
SPNTypeID's inputs are paired Illumina FASTQ files for each sample and a comma separated sample sheet containing the sample name, the path to the forward reads file, and the path to the reverse reads file for each sample. A sample sheet can be created using the [fastq_dir_to_samplesheet.py](https://github.com/wslh-bio/SPNtypeID/blob/main/bin/fastq_dir_to_samplesheet.py) script or by hand.  An example of the sample sheet's format can be seen in the table below and found [here](https://raw.githubusercontent.com/wslh-bio/SPNtypeID/main/samplesheets/workflow_test.csv).

| sample  | fastq_1 | fastq_2 |
| ------------- | ------------- | ------------- |
| sample_name  | /path/to/sample_name_R1.fastq.gz | /path/to/sample_name_R2.fastq.gz |

### Parameters
SPNTypeID's main parameters and their defaults are shown in the table below:
| Parameter  | Parameter description and default |
| ------------- | ------------- |
| contaminants  | Path to fasta of contaminants for removal, defaults to BBDuk's adapters fasta |
| maxpctother  | Sets the maximum percentage of reads from other organisms (default: 1.0) |
| minavgreadq | Sets the minimum average read quality score (default: 30) |
| mincoverage | Sets the minimum coverage (default: 40) |
| minlength | Minimum read length for trimming (default: 10) |
| minpctspn | Sets the minimum percentage of reads that must be S. pneumoniae (default: 60.0) |
| minpctstrep | Sets the minimum percentage of reads that must be Streptococcus (default: 80.0) |
| ntc_read_limit | Sets the maximum number of reads allowed in the no template control (default: 10000) |
| ntc_spn_read_limit | Sets the BBDuk trimming quality score value (default: 500) |
| qualitytrimscore | Sets the BBDuk trimming quality score value (default: 10) |
| trimdirection | Sets the BBDuk trimming direction (default: 'lr') |
| workflow_test | Run the workflow test (default: false) |
| ntc_regex | Regex pattern for identifying no template control (NTC) files. This is a mandatory parameter if a run has an NTC. (default: null) |
| split_regex | Regex pattern to split sample ID from rest of file name for fastq files with no run name (default: '_S\\\d+') |
| run_name_regex | Regex pattern for run names found in fastq file (default: 'WI-\\\D\\\d{4}-\\\d{6}[A-Za-z]*') |
| maxcontigs | Set the maximum number of contigs allowed in an assembly (default: 300) |

### Workflow outline

<img src ='/assets/SPNtypeID.png'>

#### Read trimming and quality assessment
Read trimming and cleaning is performed using [BBtools v38.76](https://jgi.doe.gov/data-and-tools/bbtools/) to trim reads of low quality bases and remove PhiX contamination. Then [FastQC v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used assess the quality of the raw and cleaned reads. [Bioawk v1.0](https://github.com/lh3/bioawk) is used to calculate the mean and median quality of the cleaned reads.

#### Genome assembly
Assembly of the cleaned and trimmed reads is performed using [Shovill v1.1.0](https://github.com/tseemann/shovill).

#### Genome assembly quality assessment
Quality assessment of the assemblies is performed using [QUAST v5.0.2](http://bioinf.spbau.ru/quast).

#### Genome coverage quality assessment
Mean and median genome coverage is determined by mapping the cleaned reads back their the assembly using [BWA v0.7.17-r1188](http://bio-bwa.sourceforge.net/) and calculating depth using [Samtools v1.10](http://www.htslib.org/).

#### Genome length assessment
Genome length is assessed by comparing the expected *S. pneumoniae* genome length to the observed genome length and calculating a Z score. These statistics (which can be found [here](/assets/databases/NCBI_Assembly_stats_20240124.txt) were obtained from the [PHoeNIx](https://github.com/CDCgov/phoenix) pipeline, which calculated them from 9266 publicly available *S. pneumoniae* genomes.

#### Contamination detection
Contamination is detected by classifying reads using [Kraken v1.0.0](https://ccb.jhu.edu/software/kraken2/).

#### Serotyping
Serotyping is performed using [SeroBA v2.0.4](https://github.com/GlobalPneumoSeq/seroba).

### Output files
Example of pipeline output:
```
outdir
├── assembly_stats_summary
│   └── assembly_stats_results_summary.tsv
├── bbduk
│   ├── *_1.fastq.gz
│   ├── *_2.fastq.gz
│   ├── *.adapter.stats.txt
│   ├── *.bbduk.log
│   └── *.trim.txt
├── bbduk_summary
│   └── bbduk_results.tsv
├── bioawk
│   └── *.qual.tsv
├── calculate_assembly_stats
│   └── *_Assembly_ratio_20240124.tsv
├── coverage_stats
│   └── coverage_stats.tsv
├── create_report
│   └── *_spntypeid_report.csv
├── fastqc
│   ├── *_fastqc.html
│   └── *_fastqc.zip
├── fastqc_summary
│   └── fastqc_summary.tsv
├── kraken_ntc
│   └── *.kraken.txt ***
├── kraken_sample
│   └── *.kraken.txt
├── kraken_summary
│   └── kraken_results.tsv
├── multiqc
│   ├── multiqc_data
│   ├── multiqc_plots
│   └── multiqc_report.html
├── percent_strep_summary
│   └── percent_strep_results.tsv
├── pipeline_info
├── quality_stats
│   └── quality_stats.tsv
├── quast
│   ├── *.quast.report.tsv
│   └── *.transposed.quast.report.tsv
├── quast_summary
│   └── quast_results.tsv
├── rejected_samples
│   └── Empty_samples.csv ***
├── samtools
│   ├── *.bam
│   ├── *.depth.tsv
│   └── *.stats.txt
├── seroba
│   ├── seroba.log
│   └── *.pred.csv
├── seroba_summary
│   └── seroba_results.tsv
└── shovill
    ├── *.contigs.fa
    ├── *.sam
    └── *_shovill_output
        ├── contigs.gfa
        ├── shovill.corrections
        ├── shovill.log
        └── spades.fasta
```
 *** = Optional output

**Notable result files:**  
**`<runname>`_spntypeid_report.csv** - Summary table of each step in SPNtypeID  
**multiqc_report.html** - HTML report generated by MultiQC  
**Empty_samples.csv** - Lists any samples that are empty and were removed from the pipeline. If no samples were empty, file will be absent from output directory.

### Results file explanation
| Output header | Purpose |
| ------------- | ------------- |
|Sample| Unique sample identifier|
|Run| Which run the sample is on, dictated by the `--runname` param |
|Total Reads| How many reads identified in the sample |
|Reads Removed| How many reads removed from the total reads |
|Median Read Quality| The median value of Phred scores |
|Average Read Quality| The mean value of Phred scores | 
|Contigs (#)| How many contigs present in the sample |
|N50| Length of the shortest contig where contigs of greater than this length contain 50% of the total bases |
|Assembly Length (bp)| Total size of the assembly |
|Ratio of Actual:Expected Genome Length| Relationship between the actual genome length to the expected S. pneumoniae genome length |
|z-score| The isolate's relationship compared to the mean S. pneumo genome length |
|Median Coverage| Median amount of times each base was sequenced  |
|Average Coverage| Mean amount of times each base was sequenced |
|Percent Strep| Percentage of reads from Streptococcus |
|Percent SPN| Percentage of reads specifically from S. pneumoniae |
|SecondGenus| Other genus present, if detected |
|Percent SecondGenus| Percentage of other genus, if detected |
|Serotype| Serotype determined for each sample |
|Kraken Database Version| Version of the Kraken database utilized |
|All NTC reads| List of all reads detected in all no template controls, if present. If '999999' in column, no NTC was provided |
|All NTC SPN reads| List of all S. pneumoniae reads in all no template controls, if present. If '999999' in column, no NTC was provided |
|Max NTC read| Highest amount of reads found in all no template controls. If '999999' in column, no NTC was provided |
|Max NTC SPN read| Highest amount of S. pneumoniae reads found in all no template controls. If '999999' in column, no NTC was provided |
|SPNtypeID Version| Version of the SPNTypeID pipeline used for analysis |
|Comments| Any notable comments from all steps throughout the pipeline |


## Credits
SPNTypeID was written by [Dr. Kelsey Florek](https://github.com/k-florek), [Dr. Abigail C. Shockey](https://github.com/AbigailShockey), and [Eva Gunawan, MS](https://github.com/evagunawan).

We thank the bioinformatics group at the Wisconsin State Laboratory of Hygiene for all of their contributions. 

## Contributions and Support
If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
If you use SPNtypeID for your analysis, please cite it using the following:

`K. Florek, E. Gunawan, & A.C. Shockey (2025). SPNtypeID (Version 1.9.2) [https://github.com/wslh-bio/SPNtypeID/tree/main].`

An extensive list of references for the tools used by Dryad can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
