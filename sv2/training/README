**Warning: for advanced users. Knowledge of python and scikit-learn is recommended**

This is a guide for users to retrain the default SV<sup>2</sup> classifiers and to [train custom supervised classifiers](Training#custom-feature-extraction).

Included in the SV<sup>2</sup> source package are the original training set and a jupyter notebook containing instructions for (re)training. 

---

## Methodology

1. Get a training set
    * [SV<sup>2</sup> Training Set](Training#sv2-training-set)
    * [Generate your own training features](Training#custom-feature-extraction)
2. [Train the models](Training#training-svm-classifiers)
3. [Add your models to SV<sup>2</sup>](Training#adding-new-classifiers-to-sv2)
4. [Genotype with your model](Training#genotyping-with-new-classifiers)

---

## SV<sup>2</sup> Training Set

The default training set is packaged with the source package:

```
$ ls sv2-VERSION/sv2/training/1kgp_training_data

1kgp_highcov_del_gt1kb.txt
1kgp_highcov_del_lt1kb.txt         
1kgp_highcov_del_malesexchrom.txt           
1kgp_highcov_dup_har.txt
1kgp_lowcov_dup_breakpoint.txt
1kgp_lowcov_dup_malesexchrom.txt
```
These files can be used for retraining in the [training SVM classifiers section](Training#training-svm-classifiers)

---

## Custom Feature Extraction

`sv2train` is a script designed for advanced users that wish to train genotyping classifiers with their own data. 

Given SV input, SV<sup>2</sup> will generate features for training a new classifier, given user-defined genotype labels.

```
$ sv2train -i <in.bam> [-b ...] [-v ...] -snv in.vcf.gz -p in.ped -o <sv2>

$ ls sv2_training_features/
    
    sv2_training_features_deletion_gt1kb.txt
    sv2_training_features_deletion_lt1kb.txt
    sv2_training_features_deletion_male_sex_chrom.txt
    sv2_training_features_duplication_breakpoint.txt
    ...
```

The header is formatted for the companion [jupyter notebook](https://github.com/dantaki/SV2/blob/master/sv2/training/sv2_training.ipynb), please do not alter it.

**VERY IMPORTANT:bangbang:**

before training, users have to populate the values in `copy_number`. The expected values for the companion [jupyter notebook](https://github.com/dantaki/SV2/blob/master/sv2/training/sv2_training.ipynb) is the following:

* Biallelic SVs

| copy_number | VCF genotype |
| ----------- | ------------ |
| 0           | 1/1 (DEL:HOM)|
| 1           | 0/1 (DEL:HET)|
| 2           | 0/0 (REF)    |
| 3           | 0/1 (DUP:HET)|
| 4           | 1/1 (DUP:HOM)|

* SVs on Male Sex Chromosomes

| copy_number | VCF genotype |
| ----------- | ------------ |
| 0           | 1 (DEL:ALT)  |
| 1           | 0 (REF)      |
| 2           | 1 (DUP:ALT)  |

The companion [jupyter notebook](https://github.com/dantaki/SV2/blob/master/sv2/training/sv2_training.ipynb) encodes genotype labels as copy number for simplicity. This is useful for users that wish to include variants with multiple alleles.

#### Examples of multiallelic SVs
| REF | ALT | Genotype | copy_number |
| ----| --- | -------- | ----------- | 
| \<CN1\> | \<CN0\>,\<CN2\>  | 2/2 | 4        |
| \<CN1\> | \<CN0\>,\<CN2\>  | 1/2 | 2        |
| \<CN1\> | \<CN2\>,\<CN3\>  | 0/2 | 4        |

---

## Training SVM Classifiers

The jupyter notebook is located in the source package here: `sv2-VERSION/sv2/training/sv2_training.ipynb`

A copy is also available on [github](https://github.com/dantaki/SV2/blob/master/sv2/training/sv2_training.ipynb)

This notebook is designed for retraining the default training set, but users can supply their own data given correct file names. 

It is important to chose a name for your classifier, this name will be later [loaded into SV<sup>2</sup>](Training#adding-new-classifiers-to-sv2)

The output of the jupyter notebook is a JSON file containing the paths to the trained classifiers. The models are saved in pickle `.pkl` files. 

**VERY IMPORTANT:bangbang:** do not alter the paths in the JSON file or the pickle files themselves.

---

## Adding New Classifiers to SV<sup>2</sup>

A JSON file containing paths to classifier models is required to add new classifiers. 

Pass the JSON file to the SV<sup>2</sup> `-load-clf` command

```
$ sv2 -load-clf myclf.json
```

This command appends new classifiers to the SV<sup>2</sup> classifier JSON file located here: `$SV2_INSTALL_LOCATION/sv2/config/sv2_clf.json`

---

## Genotyping with New Classifiers

After loading the classifiers with the `-load-clf` command, users can specify which model to genotype with the `-clf <classifier-name>` option. 

genotype with default classifiers
```
$ sv2 ... -clf default
```
genotype with a classifier named `myclf`
```
$ sv2 ...  -clf myclf
```
