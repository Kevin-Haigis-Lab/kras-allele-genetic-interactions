# The origin, distribution, and genetic interactions of *KRAS* alleles across cancer types

Joshua H. Cook, Giorgio E. M. Melloni, Doga C. Gulhan, Peter J. Park, Kevin M. Haigis

## Abstract

Mutational activation of *KRAS* promotes the initiation and progression of cancers, especially in the colorectum, pancreas, lung, and blood plasma, with varying prevalence of specific activating missense mutations. 
Although epidemiological studies connect specific alleles to clinical outcomes, the mechanisms underlying the distinct clinical characteristics of mutant *KRAS* alleles are unclear. 
Here, we analyzed 13,492 samples from these four tumor types to examine allele- and tissue-specific genetic properties associated with oncogenic *KRAS* mutations. 
The prevalence of known mutagenic mechanisms partially explained the observed spectrum of *KRAS* activating mutations, but the prevalence of most alleles across the different cancers could not be explained in this way, suggesting that biological selection underlies the tissue-specific frequencies of mutant alleles. 
Consistent with experimental studies that have identified distinct signaling properties associated with each mutant form of K-RAS, a genetic analysis revealed that each *KRAS* allele was associated with a distinct and tissue-specific comutation network. 
Moreover, we identified genetic dependencies associated with specific mutant *KRAS* alleles. 
Overall, this analysis demonstrates that the genetic interactions associated with oncogenic *KRAS* mutations are allele- and tissue-specific, underscoring the complexity that drives their clinical behaviors.

---

## Software

#### Machine Requirements

All compute was performed using the [Harvard Medical School Research Computing](https://rc.hms.harvard.edu) [High-Performance Computing environment](https://wiki.rc.hms.harvard.edu/display/O2).

- system: Linux
- release: 3.10.0-1062.el7.x86_64
- machine: x86_64
- processor: x86_64
- CPU cores: 4
- max. required RAM: 75GB
- interpreter: 64bit

#### Required Software

- gcc v.6.2.0
- R v.4.0.1 ([libraries](config/R-libraries.txt))
- python v3.7.4 ([libraries](config/python-env.yaml))
- git v.2.9.5
- conda2 v.4.2.13
- imageMagick v.6.9.1.10
- java v.jdk-1.8u112
- gsea v.3.0
- qpdf v.8.4.0
- libxml v.2.9.4
- openblas v.0.2.19
- freetype v.2.7
- pandoc v.2.1.1
- texlive v.2007
- jpeg v.9b
- gsl v.2.3
- libpng v.1.6.26
- hdf5 v.1.10.1
- gdal v.3.0.2
- tiff v.4.0.7
- glib v.2.50.2
- cairo v.1.14.6

## Reproduce the results

A [guide](reproduction-guide.md) has been created to aid in reproducing the results of this paper.
Please open an [issue](https://github.com/jhrcook/comutation/issues) if you run into any problems.

---

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-4.0-276DC3.svg?style=flat&logo=R)](https://cran.r-project.org)
[![python](https://img.shields.io/badge/Python-3.7-3776AB.svg?style=flat&logo=python)](https://www.python.org)

