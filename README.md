# GEEMiRKAT

Type: Package

Title: An adaptive microbiome regression-based kernel association test using generalised estimating equations for longitudinal microbiome studies

Version: 1.0

Author: Han Sun

Maintainer: Jianjun Zhang; Han Sun sunh529@henau.edu.cn

Imports: phyloseq, vegan, GUniFrac, PGEE, MiRKAT, devtools, ACAT, GEEMiRKAT

Description: GEEMiRKAT is used for testing associations between the microbiome and the host phenotypes from longitudinal microbiome data.

License: GPL-2

Encoding: UTF-8

LazyData: true

URL: https://github.com/SunHan5/GEEMiRKAT



## Introduction

This R package, **GEEMiRKAT**, can be used for testing associations between the microbiome and the host phenotypes adaptively from longitudinal microbiome data. It can be applied to datasets with diverse types of outcomes to study the association between diverse types of host phenotype and microbiome, such BMI (Gaussian distribution), disease status (Binomial distribution) or number of tumors (Poisson distribution). In addition, it can be also applied to cross-sectional data.



## Installation 

phyloseq:
```
BiocManager::install("phyloseq")
```

vegan:
```
install.packages("vegan")
```

GUniFrac:
```
install.packages("GUniFrac")
```

PGEE:
```
install.packages("PGEE")
```

MiRKAT:
```
install.packages("MiRKAT")
```

devtools:
```
install.packages("devtools")
```

ACAT:
```
devtools::install_github("yaowuliu/ACAT")
```

You may install `GEEMiRKAT` from GitHub using the following code: 

```
devtools::install_github("SunHan5/GEEMiRKAT", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------



## Usage
```
GEEMiRKAT(y, id, covs = NULL, Ks, model, n.perm = 5000)
```



## Arguments
* _y_ - response variable (i.e., host phenotype of interest). Exponential family of distributions (e.g., Gaussian, Binomial, Poisson) outcomes.
* _id_ - A vector for identifying the sequence of subjects/clusters of longitudinal data.
* _covs_ - covariate (e.g., age, gender). Default is covs = NULL.
* _Ks_ - Kernel matrix based on diverse distances.
* _model_ - "gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
* _n.perm_ - A number of permutations. Default is n.perm=5000. 



## Values
_$Ind.pvs_: The p-values for the individual tests based on different dissimilarity metrics and working correlation matrices.

_$GEEMiRKAT.pvs_: The p-values for the adaptive test and omnibus tests.



## Example

Import requisite R packages: 

```
library(phyloseq)
library(vegan)
library(GUniFrac)
library(PGEE)
library(MiRKAT)
library(ACAT)
library(GEEMiRKAT)
```


Import example microbiome data:

```
data(CD)
y <- CD@sam_data$CD
covs <- as.factor(CD@sam_data$smoker)
covs <- as.data.frame(covs)
id <- CD@sam_data$id
otu.tab <- CD@otu_table
tree <- CD@phy_tree
```


Calculate the similarity matrix between samples
```
unifracs <- GUniFrac(otu.tab, tree, alpha = c(0.5, 1))$unifracs
D.weighted = unifracs[, , "d_1"]
D.weighted50 = unifracs[, , paste("d_", 0.5, sep = "")]
D.unweighted = unifracs[, , "d_UW"]
D.BC = as.matrix(vegdist(otu.tab , method="bray"))

K.weighted = D2K(D.weighted)
K.weighted50 = D2K(D.weighted50)
K.unweighted = D2K(D.unweighted)
K.BC = D2K(D.BC)
Ks <- list(BC = K.BC, U.UniFrac = K.unweighted, G.UniFrac.50 = K.weighted50, W.UniFrac = K.weighted)
```

Fit GEEMiRKAT:

```
set.seed(123)
out <- GEEMiRKAT(y, covs, id = id, Ks = Ks, model = "binomial")
out
```



## References
* Zhang J, et al. An adaptive microbiome regression-based kernel association test using generalized estimating equations for longitudinal microbiome studies. (under review)

* Chen J, et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. _Bioinformatics_ 2012;**28**:2106-2113.

* Escobar JS, et al. The gut microbiota of Colombians differs from that of Americans, Europeans and Asians. _BMC Microbiology_ 2014;**14**(1):311.

* Liu Y, et al. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. _American Journal of Human Genetics_ 2019;**104**:410-421.

* Koh H, et al. A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. _Frontiers in Genetics_ 2019;**10**:458.

* McMurdie PJ and Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;**8**(4):e61217

* Oksanen J, et al. vegan: Community Ecology Package. R package version 2.5-6. 2019. 

* Vázquez-Baeza Y., et al. Guiding longitudinal sampling in IBD cohorts. _Gut_ 2018;**67**:1743-1745.

* Wang L. GEE analysis of clustered binary data with diverging number of covariates. _Ann. Statist._ 2011;**39**:389–417.

* Wang L, et al. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. _Biometrics_ 2012;**68**(2):353-360.

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.


## Statement

Our code mainly refers to R packages, _GEEMiHC_, _OMiRKAT_ and _aGLMMMiRKAT_.
