# GEEMiRKAT

Type: Package

Title: Detecting sparse microbial association signals from longitudinal microbiome data based on generalized estimating equations

Version: 1.0

Author: Han Sun

Maintainer: Jianjun Zhang ; Han Sun sunh529@henau.edu.cn

Imports: phyloseq, GUniFrac, vegan, PGEE, MiRKAT, ACAT, devtools, GEEMiHC

Description: GEEMiHC is used for detecting sparse microbial association signals between microbiome and a host phenotype from longitudinal microbiome data.

License: GPL-2

Encoding: UTF-8

LazyData: true

URL: https://github.com/SunHan5/GEEMiRKAT



## Introduction

This R package, **GEEMiHC**, can be used for detecting sparse microbial association signals adaptively from longitudinal microbiome data. It can be applied to datasets with diverse types of outcomes to study the association between diverse types of host phenotype and microbiome, such BMI (Gaussian distribution), disease status (Binomial distribution) or number of tumors (Poisson distribution). Considering cross-sectional data as a special case of longitudinal data, it can be also applied to cross-sectional data, in which case the results will be consistent with MiHC.



## Installation 

phyloseq:
```
BiocManager::install("phyloseq")
```

GUniFrac:
```
install.packages("GUniFrac")
```

vegan:
```
install.packages("vegan")
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

GEEMiHC:
```
devtools::install_github("SunHan5/GEEMiHC", force=T)
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
_$Ind.pvs_: The p-values for the individual tests based on diverse distances.

_$Omn.pvs_: The p-values for the omnibus test.



## Example

Import requisite R packages: 

```
library(phyloseq)
library(GUniFrac)
library(vegan)
library(PGEE)
library(MiRKAT)
library(ACAT)
library(GEEMiHC)
library(GEEMiRKAT)
```


Import example microbiome data:

```
data(CD_longitudinal)
otu.tab <- CD_longitudinal@otu_table
tree <- CD_longitudinal@phy_tree
y <- sample_data(CD_longitudinal)$label
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(CD_longitudinal)$age)
covs[,2] <- as.factor(sample_data(CD_longitudinal)$smoker)
id <- sample_data(CD_longitudinal)$id
```

Fit GEEMiHC:

```
set.seed(123)
out <- GEEMiRKAT(y, covs, id = id, Ks = Ks, model = "binomial")
out
```



## References
* Sun H, et al. Detecting sparse microbial association signals adaptively from longitudinal microbiome data based on generalized estimating equations. _Briefings in Bioinformatics_ 2022;**23**:bbac149.

* Koh H, et al. A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. _Frontiers in Genetics_ 2019;**10**:458.

* Koh H and Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. _Microbiome_ 2020;**8**(1):63.

* McMurdie PJ and Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;**8**(4):e61217

* Vázquez-Baeza Y., et al. Guiding longitudinal sampling in IBD cohorts. _Gut_ 2018;**67**:1743-1745.

* Wang L. GEE analysis of clustered binary data with diverging number of covariates. _Ann. Statist._ 2011;**39**:389–417.

* Wang L, et al. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. _Biometrics_ 2012;**68**(2):353-360.

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.


## Statement

Our code mainly refers to R packages, _GEEMiHC_, _OMiRKAT_ and _aGLMMMiRKAT_.
