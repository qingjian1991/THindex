# THindex
Tumor Heterogeneity Index Estimating

1. Concepts

In this package, we use the ecology methods to estimate the Tumor Heterogeneity(TH) based on their mutated loci of variant allele frequencis(VAF).

The inferHeterogeneityPlus estimate the TH based on two different methods in the Package [vegan](https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf): 1) **Diveristy indices** and 2) **Taxonomic indices**. See also <http://www.coastalwiki.org/wiki/Measurements_of_biodiversity> for the concepts.

2. Software

Required R packages including:

(1)maftools 

(2)mclust 

(3)tidyverse 

(4)data.table 

Installing the package by the `devtools`.
```R
devtools::install_github("qingjian1991/THindex")
```

