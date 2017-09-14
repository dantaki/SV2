	
![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2.png "Support Vector Structural Variation Genotyper")

Support Vector Structural Variation Genotyper

*A genotyper for the rest of us*

[![DOI](https://zenodo.org/badge/80166279.svg)](https://zenodo.org/badge/latestdoi/80166279)

## Table of Contents

* [Preprint](#preprint)
* [User Guide](#user-guide)
   * [Tutorial](#tutorial)
* [Getting Started](#getting-started)
   * [Installation](#installation)
   * [Requirements](#requirements)
   * [Configure SV<sup>2</sup>](#configure-sv2)
   * [Options](#options)
* [Input](#input)
* [Output](#output)
* [Training](#training)
* [Performance](#performance)
* [Usage](#usage)
* [Credits](#credits)
* [Citing SV<sup>2</sup>](#citing-sv2)
* [History](#history)
* [License](#license)
* [Contact](#contact)

## [Preprint](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)

SV<sup>2</sup> (support-vector structural-variant genotyper) is a machine learning algorithm for genotyping deletions and duplications from paired-end whole genome sequencing data. SV<sup>2</sup> can rapidly integrate variant calls from multiple SV discovery algorithms into a unified callset with [high genotyping accuracy](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png) and detection of [*de novo* mutations](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_fdr.png).

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_flowchart.png "Support Vector Structural Variation Genotyper Work Flow")

## [User Guide](https://github.com/dantaki/SV2/wiki)

[click here :notebook:](https://github.com/dantaki/SV2/wiki)

### [Tutorial](https://github.com/dantaki/SV2/wiki/tutorial)

## Getting Started

### [Installation](https://github.com/dantaki/SV2/wiki/installation)

[Install with `pip`](https://github.com/dantaki/SV2/wiki/installation#install-with-pip) *Recommended* 

```
pip install https://github.com/dantaki/SV2/releases/download/sv2v1.3.2/sv2-1.3.2.tar.gz 
```
*Advanced users:* SV<sup>2</sup> can be [manually installed from source](https://github.com/dantaki/SV2/wiki/installation#manual-install). 

### Requirements
* [python 2.7](https://www.python.org/)
  * [cython](https://github.com/cython/cython)
  * [numpy](http://www.numpy.org/)
  * [pandas](http://pandas.pydata.org/)
  * [pybedtools](https://daler.github.io/pybedtools/)
  * [pysam 0.9+](https://github.com/pysam-developers/pysam)
  * [scikit-learn v0.17](http://scikit-learn.org/)

* [bedtools 2.25.0](https://github.com/arq5x/bedtools2/releases) or later

### [Configure SV<sup>2</sup>](https://github.com/dantaki/SV2/wiki/installation#configure)

Before running SV<sup>2</sup>, define your fasta files. 

```
sv2 -hg19 <hg19.fasta> [-hg38 <hg38.fasta>] 
```

### [Options](https://github.com/dantaki/SV2/wiki/options#)

Please refer to the [User Guide](https://github.com/dantaki/SV2/wiki/options) for details on SV<sup>2</sup> options.

## [Input](https://github.com/dantaki/SV2/wiki/input)

### Input Arguments

| Input Arguments | Description 
| ----------------| -----------
-i \| -in       | Sample information [ID, BAM-PATH,VCF-PATH, M/F]
-b \| -bed      | Bed file(s) of SVs
-v \| -vcf      | VCF file(s) of SVs

### [Sample information < -i >](https://github.com/dantaki/SV2/wiki/input#sample-information)

Tab or space delimited file containing sample information. Gender can also be encoded as 1 for M and 2 for F

ID | BAM PATH |  VCF PATH | Gender [M/F]
--- | --- | --- | ---
NA12878 | /bam/NA12878.bam | /vcf/NA12878_SNVs.vcf.gz | F 
HG00096 | /bam/HG00096.bam | /vcf/HG00096_SNVs.vcf.gz | M

* BAM format
  * Supplementary alignment tags (SA) are required for split-read analysis
* VCF format
  * Allele Depth is required 
  * bgzip and tabix indexed VCF

Refer to the [User Guide](https://github.com/dantaki/SV2/wiki/input#sample-information) for more details. 

### [Variants to genotype <-b ... > <-v ... >](https://github.com/dantaki/SV2/wiki/input#sv-input)

SV<sup>2</sup> can accept multiple BED and VCF files.

```
sv2 -i in.txt -b del.bed dup.bed ... -v del.vcf ...
```

* [BED format](https://github.com/dantaki/SV2/wiki/input#bed-input)
  BED files are either tab or space delimited. The first four columns must be [CHROM  START  END  TYPE]

* [VCF format](https://github.com/dantaki/SV2/wiki/input#vcf-input)
  VCF files must contain `SVTYPE=` and `END=` in the INFO column 

* SV Type: must contain either `DEL` or `DUP`

Refer to the [User Guide](https://github.com/dantaki/SV2/wiki/input#sv-input) for more details. 

## [Output](https://github.com/dantaki/SV2/wiki/Output)
 
 Output is generated in the current working directory. 
 
 * `sv2_preprocessing/` contains preprocessing output. 

 * `sv2_features/` contains feature extraction output. 
 
 * `sv2_genotypes/` contains output in tab-delimited BED format and VCF format.
 
*Output VCF comes with gene annotations and other useful statistics*

### [Merging SVs](https://github.com/dantaki/SV2/wiki/Output#merging-svs)

SV<sup>2</sup> provides the option to merge SV calls. By default this option is off. 

Details found in the [User Guide](https://github.com/dantaki/SV2/wiki/Output#merging-svs)

```
# merge SV after genotyping
sv2 -i in.txt [-b ...] [-v ...] -merge

# merge SV with >50% reciprocal overlap
sv2 -i in.txt [-b ...] [-v ...] -min-ovr 0.5
```

For more detail on SV<sup>2</sup> output, please refer to the [User Guide](https://github.com/dantaki/SV2/wiki/output)
 

## [Training](https://github.com/dantaki/SV2/wiki/Training)

Advanced users can retrain SV<sup>2</sup> genotyping classifiers with the original or custom training set. 

SV<sup>2</sup> includes a script for [generating training features](https://github.com/dantaki/SV2/wiki/Training#custom-feature-extraction)

```
sv2train -i <in.txt> [-b ...] [-v ...] ...
```

Included is a [jupyter notebook guide for training classifiers](https://github.com/dantaki/SV2/blob/master/sv2/training/sv2_training.ipynb). This guide will produce a JSON file containing paths of the new classifiers. 

Users can then [load the classifiers](https://github.com/dantaki/SV2/wiki/Training#adding-new-classifiers-to-sv2) into SV<sup>2</sup> with this command:

```
sv2 -load-clf myclf.json
```

To [genotype with new classifiers](https://github.com/dantaki/SV2/wiki/Training#genotyping-with-new-classifiers), pass the name of the classifier to the SV<sup>2</sup> command

```
sv2 -i <in.txt> [-b ...] [-v ...] -clf myclf
```

For more details please refer to the [User Guide](https://github.com/dantaki/SV2/wiki/Training)

## [Performance](https://github.com/dantaki/SV2/wiki/Performance)

SV<sup>2</sup> performance was measured with independent cohorts using Illumina 2.5M arrays and PacBio SMRT sequencing. 

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png "Genotyping ROC curve")

##### [Performance of *de novo* mutations](https://github.com/dantaki/SV2/wiki/performance#de-novo-mutations)

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
