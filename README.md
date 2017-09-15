	
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
* [Input](#input)
* [Output](#output)
* [Training](#training)
* [Performance](#performance)
* [Requirements](#requirements)
* [Citing SV<sup>2</sup>](#citing-sv2)
* [License](#license)
* [Contact](#contact)

## [Preprint](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)

SV<sup>2</sup> (support-vector structural-variant genotyper) is a machine learning algorithm for genotyping deletions and duplications from paired-end whole genome sequencing data. SV<sup>2</sup> can rapidly integrate variant calls from multiple SV discovery algorithms into a unified callset with [high genotyping accuracy](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png) and detection of [*de novo* mutations](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_fdr.png).

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_flowchart.png "Support Vector Structural Variation Genotyper Work Flow")

## [User Guide](https://github.com/dantaki/SV2/wiki)

### [Tutorial](https://github.com/dantaki/SV2/wiki/tutorial)

---

## Getting Started

1. [Installation](https://github.com/dantaki/SV2/wiki/installation)

[Install with `pip`](https://github.com/dantaki/SV2/wiki/installation#install-with-pip) *Recommended* 

```
$ pip install https://github.com/dantaki/SV2/releases/download/sv2v1.3.2/sv2-1.3.2.tar.gz 
```

2. [Configure SV<sup>2</sup>](https://github.com/dantaki/SV2/wiki/installation#configure)

```
$ sv2 -hg19 /full/path/to/hg19.fasta [-hg38 /full/path/to/hg38.fasta] 
```

3. Run SV<sup>2</sup>

```
$ sv2 -i samples.txt -b sv.bed -v sv.vcf -o my_genotypes
```

[:notebook: SV<sup>2</sup> Options Documentation](https://github.com/dantaki/SV2/wiki/options#)

---

## [Input](https://github.com/dantaki/SV2/wiki/input)

SV<sup>2</sup> *requires* a sample information file and either a BED or VCF file of SVs to genotype.

1. Sample Information [:notebook: Documentation](https://github.com/dantaki/SV2/wiki/input#sample-information)
 
```
-i | -in    PATH        Sample information file
```
Tab or space delimited file containing sample information.

ID | BAM PATH |  VCF PATH | Gender [M/F]
--- | --- | --- | ---
NA12878 | /bam/NA12878.bam | /vcf/NA12878_SNVs.vcf.gz | F 
HG00096 | /bam/HG00096.bam | /vcf/HG00096_SNVs.vcf.gz | M


2. SVs to Genotype [:notebook: Documentation](https://github.com/dantaki/SV2/wiki/input#sv-input)

SV<sup>2</sup> can accept multiple BED and VCF files.

```
-b | -bed    ...       BED file(s) of SVs
```

```
-v | -vcf    ...       VCF file(s) of SVs
```

Documentation on required SV file formats:
  * [BED format](https://github.com/dantaki/SV2/wiki/input#bed-input)

  * [VCF format](https://github.com/dantaki/SV2/wiki/input#vcf-input)

---

## [Output](https://github.com/dantaki/SV2/wiki/Output)
 
 Output is generated in the current working directory. 
 
 * `sv2_preprocessing/` contains preprocessing output 

 * `sv2_features/` feature extraction output
 
 * `sv2_genotypes/` genotype output
 
*Output VCF comes with gene annotations and other useful annotations*

* SV<sup>2</sup> can merge divergent breakpoints. By default this option is off. 

* [:notebook: Merging Documentation](https://github.com/dantaki/SV2/wiki/Output#merging-svs)

---

## [Training](https://github.com/dantaki/SV2/wiki/Training)

Advanced users can retrain SV<sup>2</sup> genotyping classifiers with the original or custom training set. 

[:notebook: Training Documentation](https://github.com/dantaki/SV2/wiki/Training)

---

## [Performance](https://github.com/dantaki/SV2/wiki/Performance)

SV<sup>2</sup> performance was measured with independent cohorts using Illumina 2.5M arrays and PacBio SMRT sequencing. 

![alt text](https://raw.githubusercontent.com/dantaki/SV2/master/png/sv2_auc.png "Genotyping ROC curve")

##### [Performance of *de novo* mutations](https://github.com/dantaki/SV2/wiki/performance#de-novo-mutations)

Please refer to the [preprint](#preprint) for performance details. 

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

[:notebook: Usage Documentation](https://github.com/dantaki/SV2/wiki/usage)

---

## Citing SV<sup>2</sup>

For citing SV<sup>2</sup> please refer to the preprint: [bioRxiv](http://biorxiv.org/content/early/2017/03/17/113498) : [doi](https://doi.org/10.1101/113498)

---

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

---

## Contact
:mailbox:
dantaki@ucsd.edu
:metal: