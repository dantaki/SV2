# SV<sup>2</sup> User Guide

## Table of Contents
1. [Introduction](#introduction)
   * [Preprocessing](#preprocessing)
   * [Feature Extraction](#feature-extraction)
   * [Genotyping](#genotyping)
2. [Installation](#installation)
   * [Prerequisites](#prerequisites)
   * [`pip` install](#install-with-pip)
   * [Manual Install](#manual-install)
   * [Configure](#configure)
3. [Input](#input)
   * [Sample Information](#sample-information)
      * [BAM](#bam)
      * [SNV VCF](#snv-vcf)
   * [SV Input](#sv-input)
      * [BED Input](#bed-input)
      * [VCF Input](#vcf-input)
   * [FASTA](#fasta)
4. [Output](#output)
   * [Preprocessing Output](#preprocessing-output)
   * [Feature Output](#feature-output)
   * [Genotyping Output](#genotyping-output)
      * [VCF Output](#vcf-output)
      * [BED Output](#bed-output)
5. [Performance](#performance)
   * [Genotyping Accuracy](#genotyping-accuracy)
   * [*De Novo* Mutations](#de-novo-mutations)
6. [Tutorial](#tutorial)
   * [Overview](#overview)
   * [Download](#download-tutorial-files)
   * [Generate Input](#generate-input)
   * [Genotype](#genotype-sv)
7. [Troubleshooting](#troubleshooting)

## Introduction

SV<sup>2</sup> (support-vector structural-variant genotyper) is a machine learning algorithm for genotyping deletions and duplications from paired-end whole genome sequencing data. SV<sup>2</sup> can rapidly integrate variant calls from multiple SV discovery algorithms into a unified callset with high genotyping accuracy and detection of *de novo* mutations. 

SV<sup>2</sup> is an open source software written in Python/Cython that exploits four features of SV: read depth, discordant paired-ends, split-reads, and heterozygous allele depth in a supervised support vector machine classifier. 

Required inputs include a BAM file with supplementary alignment tags (SA), a SNV VCF file with allele depth, and either a BED or VCF of SV to genotype. SV<sup>2</sup> operates in three stages: preprocessing, feature extraction , and genotyping. SV<sup>2</sup> outputs a BED and VCF file of genotypes along with annotations for genes, repeats, and other befitting statistics for SV analysis. 

For more information or citing SV<sup>2</sup> please refer to the [bioRxiv preprint](http://biorxiv.org/content/early/2017/03/17/113498.article-metrics).

### Preprocessing

SV<sup>2</sup> preprocessing records the median coverage, insert size, read length for each chromosome for downstream normalization of features. Preprocessing statistics are obtained for each chromosome in 100 random nonoverlapping regions 100kb in length. The default random seed is 42, but can be altered with the `[-s|-seed] INT` option.

[Preprocessing output](#preprocessing-output) can be found in `sv2_preprocessing/` in the current working directory. If preprocessing has been completed, supplying `-pre sv2_preprocessing/` to the SV<sup>2</sup> command will bypass this stage.
 
### Feature Extraction

Before feature extraction, a mask is applied to SV regions. The mask includes segmental duplications, short tandem-repeats, centromeres, telomeres, and unmappable regions and is included here `$SV2_INSTALL_DIR/sv2/src/resources/annotation_files/hg*flags.bed.gz`. The genome mask with merged positions is located here `$SV2_INSTALL_DIR/sv2/src/resources/hg*_unmapped.bed.gz`. SV calls that completely overlap masked elements cannot be genotyped and are represented as `./.` in the output VCF.

SV<sup>2</sup> genotypes using four features of SV: depth of coverage, discordant paired-end, split-reads, and heterozygous allele depth. [Feature extraction output](#feature-output) is located in `sv2_features/` in the current working directory. If feature extraction has completed and a VCF of multiple samples is desired, supplying the option `-feats sv2_features/` to the SV<sup>2</sup> command with skip feature extraction. 

#### Depth of Coverage

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/doc.png "Depth of Coverage")

Depth of coverage is estimated via the number of reads spanning a locus for SV greater than 1000bp. For smaller SVs, depth of coverage was recorded as the median per-base-pair read depth. 

Coverage features are first normalized according to the median chromosome coverage. For SVs overlapping pseudoautosomal regions on male sex chromosomes, normalization implements the median genome coverage. Normalized coverage is then corrected for GC content, adapted from [CNVator](http://genome.cshlp.org/content/21/6/974.long), for either PCR or PCR-free libraries. 

For PCR-free libraries supply the `-pcrfree` flag.

SV<sup>2</sup> cannot genotype SVs when normalized coverage exceeds 5.0 (10 autosomal copies).

#### Discordant Paired-Ends

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/dpe.png "Discordant Paired-Ends")

Discordant paired-ends contain insert sizes greater than the chromosome median plus five times the median absolute deviation. SV<sup>2</sup> only considers discordant paired-ends if both mates bridge the putative breakpoint by +/- 500bp. Likewise, SV<sup>2</sup> requires that both mates rest on opposite sides of the breakpoint. The resulting number of discordant paired-ends that meet these criteria are then normalized by the number of concordant paired-ends that span 500bp windows of the start and end positions of the SV.
   
#### Split-Reads

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sr.png "Split-Reads")

Split-reads are those with supplementary alignments. To reduce noise, SV<sup>2</sup> only considers split-reads if the primary and supplementary alignments bridge the breakpoint by +/- 500bp. Likewise, both alignments must map to opposite sides of the breakpoint. The resulting number of split-reads are normalized to the number of concordant reads that span 500bp windows of the start and end positions of the SV. 

#### Heterozygous Allele Depth

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/had.png "Heterozygous Allele Depth")

Akin to B-allele frequency in microarrays, heterozygous allele depth is defined as the median ratio of minor allele reads to major allele reads for every heterozygous SNV within the SV. This feature is parsed from a SNV VCF that contains either `AD` or `DPR` in the format column.  

### Genotyping

SV<sup>2</sup> genotypes SV with six [support vector machine classifiers](http://scikit-learn.org) that are trained with respect to SV type and length.

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/clf.png "SVM Classifiers")

Each classifier, with the exception of Duplication SNV implements depth of coverage, discordant paired-ends, and split-reads as features. The Duplication SNV classifier employs depth of coverage and heterozygous allele depth as features. 

[Genotyping output](#genotyping-output) in BED and VCF format are located in `sv2_genotypes/` in the current working directory.

## Installation

### Prerequisites
* [python 2.7](https://www.python.org/)
  * [cython](https://github.com/cython/cython)
  * [numpy](http://www.numpy.org/)
  * [pandas](http://pandas.pydata.org/)
  * [pybedtools](https://daler.github.io/pybedtools/)
  * [pysam 0.9+](https://github.com/pysam-developers/pysam)
  * [scikit-learn 0.17](http://scikit-learn.org)
* [bedtools 2.25.0](https://github.com/arq5x/bedtools2/releases) or later


### Install with `pip`

*Recommended*

```
pip install http://downloads.sourceforge.net/project/sv2/sv2-1.2.tar.gz
```

### Manual Install

> [Source Files :floppy_disk:](#source-files)

```
wget http://downloads.sourceforge.net/project/sv2/sv2-1.2.zip # sv2-1.2.tar.gz also available
unzip sv2-1.2.zip
cd sv2-1.2/

python setup.py install [--prefix <PYTHONPATH>]
```
* define `--prefix <PYTHONPATH>` for local installation
* ignore numpy compilation warnings

### Configure

SV<sup>2</sup> requires one FASTA file to run.

* [hg19 FASTA](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)

* [hg38 FASTA](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)

```
sv2 -hg19 <hg19.fasta> -hg38 <hg38.fasta>
```

### Check Installation


```
usage: sv2 [-h] [-i I] [-r R] [-c C] [-g G] [-pcrfree] [-s S] [-o O]
           [-pre PRE] [-feats FEATS] [-hg19 HG19] [-hg38 HG38]

                       ____
  _____________   ___ |___ \
 /   _____/\   \ /   // ___/
 \_____  \  \   Y   //_____)
 /        \  \     /
/_________/   \___/
Support Vector Structural Variation Genotyper
Version 1.2        Author: Danny Antaki <dantaki at ucsd dot edu>

optional arguments:
  -h, --help       show this help message and exit

genotype arguments:
  -i I, -in I      Tab delimited input [ ID, BAM-PATH, VCF-PATH, M/F ]
  -r R, -cnv R     SV to genotype. Either in BED or VCF format
  -c C, -cpu C     Parallelize sample-wise. 1 per cpu
  -g G, -genome G  Reference genome build [ hg19, hg38 ]
  -pcrfree         GC content normalization for PCR free libraries
  -s S, -seed S    Preprocessing: integer seed for genome shuffling
  -o O, -out O     output
  -pre PRE         Preprocessing output directory
  -feats FEATS     Feature extraction output directory

configure arguments:
  -hg19 HG19       hg19 FASTA
  -hg38 HG38       hg38 FASTA
```

If you get this error: `Error detail: Resource temporarily unavailable` there is not enough memory to run SV<sup>2</sup> in your shell.

## Input

The SV<sup>2</sup> command requires two input files: sample information and a list of SVs to genotype. The sample information input file contains paths to BAM and VCF files. 

A FASTA file for hg19 or hg38 is required for SV<sup>2</sup> to run. Refer to the [config](#configure) section for details.  

### Sample Information

`[-i|-in] SAMPLE_INFORMATION.txt` 

ID | BAM PATH | VCF PATH | Gender [M/F]
--- | --- | --- | ---
NA12878 | /bam/NA12878.bam | /vcf/NA12878_SNVs.vcf.gz | F
HG00096 | /bam/HG00096.bam | /vcf/HG00096_SNVs.vcf.gz | M

The sample information file must be tab-delimited with four columns. Failure to format this file correctly will result in a `Segmentation Fault`.

Each row of the sample information file corresponds to a sample. The sample information file can contain multiple samples and can be run in parallel given `[-c|-cpu NUMBER_OF_THREADS]`. 

The first column is the sample identifier. The second column contains the full path to the BAM file for that sample. Likewise, the third column contains the full path to the SNV VCF file. The SNV VCF can contain multiple samples, however the sample identifier in the first column must uniquely match one of the sample identifiers in the VCF header. The final column has the gender. Either a character `[M|F]` or an integer `[1|2]` is recognized. Integer encoding for gender is the same as PLINK fam files: 1 for male, 2 for female. Unknown genders are not supported. The purpose of including gender is to properly normalize coverage for males on sex chromosomes. 
 
#### BAM

BAM files must contain supplementary alignment tags (SA) to recognize split-reads. Without SA tags, SV<sup>2</sup> may not genotype variants well. SV<sup>2</sup> supports [BWA-MEM](http://github.com/lh3/bwa) alignments.

BAM files must also be indexed with the `.bai` file in the same directory as the BAM file. 

#### VCF

SNV VCF files must contain allele depth encoded as either `AD` or `DPR` in the format column. Additionally, VCF files must be compressed with [bgzip](http://www.htslib.org/doc/tabix.html) and indexed with [tabix](http://www.htslib.org/doc/tabix.html) with the index file in the same path as the VCF. 

When defining a VCF for a given sample, the sample identifier in the VCF header must uniquely match to the sample identifier defined in the SV<sup>2</sup> input.

### SV Input

`[-r|-cnv] SV_LIST` 

SV<sup>2</sup> requires either a BED or VCF file of deletions and duplications.

#### BED Input

BED files must be tab-delimited with at least four columns, extra columns are ignored. SV type must be contain either `DEL` or `DUP`.

The first four columns must be formatted in this order below,

CHROM | START | END | SV TYPE
---- | ---- | ----- | -----
chr1 | 100 | 200 | DEL
chr1 | 500 | 1000 | DUP

Failure to format this file correctly will result in a `Segmentation Fault`.

#### VCF Input

SV<sup>2</sup> only supports uncompressed VCF files. VCF files must contain both `END=` and `SVTYPE=` in the INFO column. `SVTYPE=` must contain either `DEL` or `DUP`.

### FASTA

SV<sup>2</sup> requires a FASTA file in either hg19 (GRCh37) or hg38 for GC content correction. 
 
* [hg38 FASTA](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)

* [hg19 FASTA](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)

Configure SV<sup>2</sup> with `sv2 -hg19 <hg19.fasta> -hg38 <hg38.fasta>`

## Output

All output is generated in the current working directory. SV<sup>2</sup> generates output files for each stage of genotyping. 

### Preprocessing Output

Preprocessing output is located in `sv2_preprocessing/` in the current working directory. Preprocessing can be skipped if a directory containing SV<sup>2</sup> preprocessing output via the `-pre sv2_preprocessing/` argument. 

##### Example Preprocessing Output

ID | CHROM | COVERAGE_MEDIAN | READ_LENGTH_MEDIAN | INSERT_SIZE_MEDIAN | INSERT_SIZE_MAD | BAM_BP_PARSED | SNP_COVERAGE_MEDIAN | SNP_PARSED
---- | ---- | -------------- | ------------------ | ------------------ | --------------- | ------------- | ------------------- | ----------
HG00096 | chr21 | 5.2 | 250.0 | 449.0 | 109.711 | 6323985 | 45.0 | 61259

### Feature Output

Feature extraction output is located in the `sv2_features/` directory in the current working directory. Feature extraction can be skipped if a directory containing SV<sup>2</sup> feature extraction output via the `-feats sv2_features/` argument. This is useful for merging genotypes for many samples. 

### Genotyping Output

SV<sup>2</sup> outputs genotypes in BED format and VCF format in the `sv2_genotypes/` directory in the current working directory. Output file name is defined with the `[-o|-out] OUTPUT_NAME` option, the default output name is `sv2_genotypes`.

VCF output contains recommended standard filters and stringent filters for *de novo* mutation discovery. 

#### VCF Output

Median Phred-adjusted ALT genotype likelihood scores are located in the `QUAL` column; median REF genotype likelihood scores are located in the `INFO` column prefixed by `REF_GTL=`.

`PASS` in the `FILTER` column represent SV<sup>2</sup> standard filters along with additional filters typically used in SV analysis, such as greater than 50% overlap to segmental duplications.

Stringent filters for *de novo* mutation discovery can be found in the `INFO` column prefixed with `DENOVO_FILTER=`.

Overlap to 1000 Genomes Phase 3 deletions and duplications can be found in the `INFO` column prefixed with `1000G_ID=` and `1000G_OVERLAP=` 

Genic overlap for each SV is recorded in the `INFO` column prefixed with `GENES=`. Transcripts are pipe-delimited with details such as the refGene name delimited by commas. 

Variants that cannot be genotyped are reported as  `./.`.
 
#### BED Output

BED output does not contain filtering information or annotations that are found in the VCF output. Both BED and VCF output are given the same file name, but with different extensions. The BED extension is `.txt`.

## Performance

For details of SV<sup>2</sup> genotyping performance, please refer to the [preprint](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)

### Genotyping Accuracy

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png "SV2 Genotyping Accuracy")

SV<sup>2</sup> has superior genotyping accuracy when genotyping the union of SVTyper and Manta SVs. Truth sets were generated from SNV arrays. Please refer to the [preprint](http://biorxiv.org/content/early/2017/03/17/113498) for more details.

### *De Novo* Mutations

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_fdr.png "SV2 De Novo Filters")

SV<sup>2</sup> provides stringent filters for *de novo* mutation discovery in the output VCF. The filters are located in the `INFO` column as `DENOVO_FILTER=`. SV<sup>2</sup> does not call variants as *de novo*, rather SV<sup>2</sup> `DENOVO_FILTER=` is a guide for filtering putative *de novo* mutations. For more details on how *de novo* filters were constructed please refer to the [preprint](http://biorxiv.org/content/early/2017/03/17/113498).

## Tutorial

### Overview
Run SV<sup>2</sup> to genotype SV called on chromosome 21 in two samples from the [1000Genomes Project](http://www.internationalgenome.org/). 

Included in the tutorial are down-sampled (10%) BAM files aligned using [BWA-MEM](https://github.com/lh3/bwa) to hg19. 

**Note:** since the BAM files are down-sampled, this tutorial is merely an example on how to execute SV<sup>2</sup>

Also included are SNV calls produced with [GATK Haplotype Caller](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php). The VCFs are [bgzipped and tabix indexed](http://www.htslib.org/doc/tabix.html). SNV calls were generated using high coverage alignments not the subsampled alignments.
 
[ForestSV](https://sites.google.com/site/sebatlab/software-data) deletion and duplication predictions are provided in BED format. 

[1000 Genomes Phase 3 integrated SV](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map) calls are included in VCF format. 

A hg19 FASTA file is required to complete this tutorial.

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz

# index FASTA
samtools faidx human_g1k_v37.fasta

# configure SV<sup>2</sup>
sv2 -hg19 human_g1k_v37.fasta
``` 

**Ensure you have properly [installed](#installation) and [tested](#check-installation) SV<sup>2</sup> before starting the tutorial.**

### Download Tutorial Files

```
wget http://downloads.sourceforge.net/project/sv2/sv2_tutorial.zip
unzip sv2_tutorial.zip
cd sv2_tutorial/
```

#### Tutorial Files

##### Alignments 
* HG00096_chr21_sub.bam 
* HG00096_chr21_sub.bam.bai
* HG01051_chr21_sub.bam
* HG01051_chr21_sub.bam.bai

##### SNV VCFs
* HG00096_chr21.vcf.gz
* HG00096_chr21.vcf.gz.tbi
* HG01051_chr21.vcf.gz
* HG01051_chr21.vcf.gz.tbi

##### SV Predictions to Genotype
* chr21_forestSV.bed
* ALL.wgs.integrated_sv_map_v1.20130502.chr21.sv.genotypes.vcf

##### Input Script
* make_sv2_input.py

### Generate Input

`make_sv2_input.py` is executed in the `sv2_tutorial/` directory. It will generate the sample input files used in this tutorial.

```
python make_sv2_input.py
ls *sv2_input.txt

# HG00096_sv2_input.txt
# HG01051_sv2_input.txt
# sv2_input.txt
```
### Genotype SV

#### Standard Usage

##### Run SV<sup>2</sup> for individual samples

```
sv2 -i HG00096_sv2_input.txt -r chr21_forestSV.bed -o HG00096_sv2
sv2 -i HG01051_sv2_input.txt -r chr21_forestSV.bed -o HG01051_sv2
```

##### Run SV<sup>2</sup> for multiple samples

```
sv2 -i sv2_input.txt -r chr21_forestSV.bed -o sv2_forestSV
```

##### Parallelize by Sample

When given multiple samples, SV<sup>2</sup> can run samples in parallel reducing run time.

```
sv2 -i sv2_input.txt -r chr21_forestSV.bed -c 2 -o sv2_forestSV
```

#### Genotype SV: Skip Preprocessing

Skipping preprocessing allows the user to genotype another set of SVs with the same samples.

**Note:** Unless you move the features files to another location they will be replaced!

```
# move ForestSV features to another location, if you generated them

mv sv2_features/ sv2_forestSV_features/

sv2 -i sv2_input.txt -r ALL.wgs.integrated_sv_map_v1.20130502.chr21.sv.genotypes.vcf -pre sv2_preprocessing/ -o sv2_tut_1KGP

```

#### Genotype SV: Skip Feature Extraction

Skipping feature extraction allows the user to combine or subset samples in the final output.

```
head -n 1 sv2_input.txt >sv2_subset.txt

sv2 -i sv2_input.txt -r chr21_forestSV.bed -pre sv2_preprocessing/ -feats sv2_forestSV_features/ -o sv2_tut_forestSV_subset

```

## Troubleshooting

Since SV<sup>2</sup> uses Cython, error messages may be cryptic. A `Segmentation Fault` typically indicates the input files are not formatted correctly. 

If you encounter this error: `UserWarning: Trying to unpickle estimator SVC from version pre-0.18 when using version 0.18.1. This might lead to breaking code or invalid results. Use at your own risk.`, install scikit-learn v0.17 `pip install scikit-learn==0.17`

For any bugs or errors in SV<sup>2</sup> please contact Danny Antaki <dantaki@ucsd.edu>

