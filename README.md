
<!-- README.md is generated from README.Rmd. Please edit that file -->
dbr
===

[![Travis build status](https://travis-ci.org/chumbleycode/dbr.svg?branch=master)](https://travis-ci.org/chumbleycode/dbr) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/chumbleycode/dbr?branch=master&svg=true)](https://ci.appveyor.com/project/chumbleycode/dbr)

The goal of dbr is to ...

Installation
------------

You can install dbr with:

``` r
install.packages("devtools")
devtools::install_github("chumbleycode/dbr")
```

Usage
-----

The Bioconductor pkgs can be installed with:
============================================

``` r

install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8")
BiocManager::install("limma", version = "3.8")
BiocManager::install("GEOquery", version = "3.8")
```

the CRAN pkgs can be installed thus:
====================================

``` r
install.packages(tidyverse)
install.packages(broom)
install.packages(data.table)
install.packages(matrixStats)
```

This is a basic example which shows you how to solve a common problem:
======================================================================

``` r

########################################################
# LOAD PACKAGES
########################################################

library(tidyverse)
library(broom)
library(data.table)
library(matrixStats)
library(Biobase)
library(limma)
library(GEOquery)

# DIFFERENTIAL MOTIF BINDING
library(dbr)
########################################################
# DOWNLOAD OPEN SOURCE DATA THEN SPECIFY GENE-BY-GENE REGRESSION MODEL
########################################################

dat = getGEO("GSE77164")[[1]]
 
y <- dat %>% exprs
X <-
  dat %>%
  pData %>%
  select(age = `age:ch1`,
         soldier = `childsoldier:ch1`,
         edu = `educationlevel:ch1`)
X     <- model.matrix(~age + soldier + edu, data = X)
of_in <- "soldier1"
# "of_in" = scalar parameter of interest: currently must identify a *single*
# column of "design", i.e. db is not currently sufficiently general to handle
# general factor vectors.

########################################################
# ESTIMATE MODEL USING STANDARD limmma/edgeR PIPELINE
# BUT STORE AS A TIDY TOP TABLE (tT) WITH A ROW CALLED "gene"
########################################################

ttT <-
  lmFit(y, X) %>%
  eBayes %>%
  topTable(coef = of_in, n = Inf) %>%
  tidy_topTable

########################################################
# DO DB INFERENCE AND EXTRACT RESULTS.
# THE FUNCTION "extracter" GIVES A TABLE OF P-VALUES FOR DIFFERENT METHODS:
# p_uni  = SIMPLE UNIVARIATE REGRESSION OF B ON EACH TFBM
# p_cov  = MULTPLE REGRESSION OF B ON ALL TFBM's SIMULTANIOUSLY
# p_tot  = SIMPLE UNIVARIATE OF B ON THE TOTAL NUMBER SITES (NOT SPECIFIC TO ONE TFBM)
# p_par  = parametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
# p_npar = nonparametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
########################################################

# BECAUSE tT_sub IS NOT SPECIFIED THESE TWO EXAMPLES JUST CALCULATE
# UNIVARIATE (UNI) AND MULTIPLE (COV) REGRESSIONS OF ON TFBM MATRIX
# ONLY EXAMINE ONE TFBM: AR_
# IN THIS CASE ALL THREE REGRESSIONS ARE EQUIVALENT
ttT %>%
  infer_db(which_tfbms = "AR") %>%
  extract_db

# THE OTHER EXTREME (ALL TFBMS): BEWARE MULTIPLICITY
ttT %>%
  infer_db %>%
  extract_db

# NOTE: p_tot IS JUST ONE REGRESSION (OF DIFFERENTIAL EFFECT), THE RESULT IS
# REPLICATED FOR CONVENIENCE

########################################################
# TELIS AND REGRESSION:
# A MORE ELABORATE EXAMPLE WITH OPTIONAL PARAMETERS
########################################################
# HEURISTICAL FILTER ON UNCORRECTED P.VALUE, OR logFC, NOT AN INFERENCE
# OPTIONALLY SELECT TFBM MATRIX "which_matrix"
# OPTIONALLY SELECT INTERESTING TFBMS "which_tfbms"
# OPTIONAL SET OF METHODS (p_npar, p_par, which require tT_sub be specified)

# db
ttT %>%
  infer_db(ttT_sub        = filter(ttT, P.Value <= 0.05),
           which_matrix   = c("exonic1_utr1") ,
           which_tfbms    =  c("ALX3", "ALX4_TBX21", "AR") ,
           n_sim          = 10000) %>%
  extract_db(methods =  c("p_npar", "p_par"))

```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
