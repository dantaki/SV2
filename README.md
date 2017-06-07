	
![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2.png "Support Vector Structural Variation Genotyper")

Support Vector Structural Variation Genotyper

*A genotyper for the rest of us*

[![DOI](https://zenodo.org/badge/80166279.svg)](https://zenodo.org/badge/latestdoi/80166279) [![Docker Automated buil](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg)](https://hub.docker.com/r/dantaki/sv2/)

## Table of Contents
* [Preprint](#preprint)
* [User Guide](#user-guide)
   * [Tutorial](#tutorial)
* [Installation](#getting-started)
* [Options](#options)
* [Input](#input)
* [Output](#output)
* [Performance](#performance)
* [Usage](#usage)
* [Requirements](#requirements)
* [Source Files](#source-files)
* [Credits](#credits)
* [Citing SV<sup>2</sup>](#citing-sv2)
* [History](#history)
* [License](#license)
* [Contact](#contact)

## Preprint

[bioRxiv](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)

SV<sup>2</sup> (support-vector structural-variant genotyper) is a machine learning algorithm for genotyping deletions and duplications from paired-end whole genome sequencing data. SV<sup>2</sup> can rapidly integrate variant calls from multiple SV discovery algorithms into a unified callset with [high genotyping accuracy](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png) and detection of [*de novo* mutations](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_fdr.png).

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_flowchart.png "Support Vector Structural Variation Genotyper Work Flow")

## [User Guide][UserGuide]

[click here :notebook:][UserGuide]

### [Tutorial](https://github.com/dantaki/SV2/blob/master/doc/README.md#tutorial)

[UserGuide]:doc/README.md

## Getting Started
#### 1: Download Source Files :floppy_disk:
```
wget http://downloads.sourceforge.net/project/sv2/sv2-1.2.zip # sv2-1.2.tar.gz also available
unzip sv2-1.2.zip
```

> [Source Files :floppy_disk:](#source-files)

#### 2: Compile and Install from Source 
```
cd sv2-1.2/
python setup.py install # ignore numpy warnings
```

#### 3: Configure SV<sup>2</sup>

```
sv2 -hg19 <hg19.fasta> [-hg38 <hg38.fasta>] 
```
SV<sup>2</sup> requires one FASTA file to run.

## Options
`sv2 --help`

Genotyping Option | Description
----------------- | -------------
-i \| -in | Tab-delimited input [ID, BAM path, VCF path, Gender]
-r \| -cnv | SV to genotype. BED or VCF
-c \| -cpu | Parallelize sample-wise. 1 per CPU 
-g \| -genome | Reference genome build [hg19, hg38]. Default: hg19
-pcrfree | GC content normalization for PCR-free libraries
-s \| -seed | Random seed for genome shuffling in preprocessing. Default: 42
-o \| -out | Output name
-pre | Preprocessing output directory. *Skips preprocessing*
-feats | Feature output directory. *Skips feature extraction*

| Configure Option | Description |
| ---------------- | ---------- |
| -hg19 | hg19 FASTA file | 
| -hg38 | hg38 FASTA file |

## Input
### Sample information < -i >
Tab-delimited file containing sample information. Gender can also be encoded as 1 for M and 2 for F

ID | BAM PATH |  VCF PATH | Gender [M/F]
--- | --- | --- | ---
NA12878 | /bam/NA12878.bam | /vcf/NA12878_SNVs.vcf.gz | F 
HG00096 | /bam/HG00096.bam | /vcf/HG00096_SNVs.vcf.gz | M

* BAM format
  * Supplementary alignment tags (SA) are required for split-read analysis
* VCF format
  * Allele Depth (AD) is required 
  * bgzip and tabix indexed VCF

Refer to the [User Guide](https://github.com/dantaki/SV2/tree/master/doc#sample-information) for more details. 

### Variants to genotype < -r >
* BED format
  * Tab-delimited: first four columns 
    * Chromosome
    * Start
    * End
    * Type: DEL | DUP
* VCF format
  * SVTYPE= DEL | DUP
  * Must have END=

Refer to the [User Guide](https://github.com/dantaki/SV2/tree/master/doc#sv-input) for more details. 

 ## Output
 
 Output is generated in the current working directory. 
 
 `sv2_preprocessing/` contains preprocessing output. `sv2_features/` contains feature extraction output. 
 
 `sv2_genotypes/` contains output in tab-delimited BED format and VCF format.
 
 *Output VCF comes with gene annotations and other useful statistics*

 For more detail on SV<sup>2</sup> output, please refer to the [User Guide](https://github.com/dantaki/SV2/tree/master/doc#output)
 
## Performance

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png "Genotyping ROC curve")

##### [Performance of *de novo* mutations](https://github.com/dantaki/SV2/blob/master/doc/README.md#de-novo-mutations)

Please refer to the [preprint](#preprint) for performance details. 

## Usage

* SV<sup>2</sup> is designed for human whole genome short-read sequencing libraries. Given deletion and duplication positions, SV<sup>2</sup> returns a VCF with predicted copy number genotypes.
* Whole genome alignments from the [1000 Genomes Project](http://www.1000genomes.org/) were used for training. Validated genotypes were obtained from the phase 3 integrated structural variation call set ([DOI:10.1038/nature15394](http://dx.doi.org/10.1038%2Fnature15394); PMID:    26432246).
* Features for genotyping include coverage, discordant paired-ends, split-reads, and heterozygous allele depth ratio.
   * BAM files must have supplementary alignment tags (SA).
   * SNV VCF must contain Allele Depth (AD/DPR). SV<sup>2</sup> can accommodate [GATK Haplotype Caller](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) and [FreeBayes](https://github.com/ekg/freebayes) VCFs.
      * SNV VCF must be compressed and indexed with [bgzip and tabix](http://www.htslib.org/doc/tabix.html)
* SV<sup>2</sup> operates with a bi-allelic model with a copy number range of 0-4
* Output is in VCF format.
   * Median Phred-adjusted ALT likelihoods are reported in the QUAL column
   * SV<sup>2</sup> standard filters are reported in the FILTER column
   * SV<sup>2</sup> stringent filters for *de novo* discovery are located in the INFO column as `DENOVO_FILTER=`
   * Positions are annotated based on their overlap to genes, RepeatMasker, segmental duplications, 1000 Genomes phase 3 CNV, and more
* SVs with estimated autosome copy number >10 cannot be genotyped. 

---

## Requirements
* [python 2.7](https://www.python.org/)
  * [cython](https://github.com/cython/cython)
  * [numpy](http://www.numpy.org/)
  * [pandas](http://pandas.pydata.org/)
  * [pybedtools](https://daler.github.io/pybedtools/)
  * [pysam 0.9+](https://github.com/pysam-developers/pysam)
  * [scikit-learn v0.17](http://scikit-learn.org/)

* [bedtools 2.25.0](https://github.com/arq5x/bedtools2/releases) or later

*SV<sup>2</sup> requires python 2.7*

*SV<sup>2</sup> has been tested on Linux and MacOS with [bioconda](https://bioconda.github.io/)*

## Source Files

[Source Forge](https://sourceforge.net/projects/sv2/files/)

[GitHub](https://github.com/dantaki/SV2/releases/tag/SV2v1.2)

*Please do not use `git clone` on this repository*

## Credits

### Author:

* Danny Antaki
    * dantaki@ucsd.edu

### Acknowledgements:
* William Brandler
* Jonathan Sebat
    * Sebat Lab http://sebatlab.ucsd.edu/index.php/software-data

## Citing SV<sup>2</sup>

For citing SV<sup>2</sup> please refer to the preprint: [bioRxiv](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)


## History

[SV<sup>2</sup> version 1.1](https://github.com/dantaki/SV2/tree/v1.1) used in Brander, Antaki, Gujral,  et al. *bioRxiv* 2017: [DOI](http://biorxiv.org/content/early/2017/04/04/102327)

[gtCNV version 0.1](https://github.com/dantaki/gtCNV/tree/Version-0.1) used in Brander, Antaki, Gujral,  et al. *AJHG* 2016: [DOI](http://dx.doi.org/10.1016/j.ajhg.2016.02.018) PMID:    27018473

## License 
MIT License

Copyright (c) 2017 Danny Antaki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact
:mailbox:
dantaki@ucsd.edu
:metal:
