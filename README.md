
<!-- README.md is generated from README.Rmd. Please edit that file -->
dbr
===

<!-- [![Travis build status](https://travis-ci.org/chumbleycode/dbr.svg?branch=master)](https://travis-ci.org/chumbleycode/dbr) -->
<!-- [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/chumbleycode/dbr?branch=master&svg=true)](https://ci.appveyor.com/project/chumbleycode/dbr) -->
The goal of dbr is to ...

Installation
------------

You can install dbr with:

``` r
install.packages("devtools")
devtools::install_github("chumbleycode/dbr")
```

This is a basic example which shows you how to solve a common problem.
----------------------------------------------------------------------

First install the relevant packages

``` r
# Bioconductor pkgs can be installed with:
install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8")
BiocManager::install("limma", version = "3.8")
BiocManager::install("GEOquery", version = "3.8")

# CRAN pkgs can be installed thus:
install.packages(tidyverse)
install.packages(broom)
install.packages(data.table)
install.packages(matrixStats)
```

Then do differential expression and binding analysis. The latter asks whether variation in differential RNA expression over genes is related to variation in the binding site density of a particular transcription factor.

``` r

# LOAD PACKAGES
library(tidyverse)
library(limma)
library(dbr)

########################################################
# DOWNLOAD OPEN SOURCE DATA THEN SPECIFY GENE-BY-GENE REGRESSION MODEL
########################################################

dat = GEOquery::getGEO("GSE77164")[[1]]
 
y <- dat %>% Biobase::exprs()
X <-
  dat %>%
  Biobase::pData() %>%
  select(age = `age:ch1`,
         soldier = `childsoldier:ch1`,
         edu = `educationlevel:ch1`)
X     <- model.matrix(~ soldier + edu + age, data = X)

########################################################
# ESTIMATE MODEL USING STANDARD limmma/edgeR PIPELINE
# BUT STORE AS A TIDY TOP TABLE (tT) WITH A ROW CALLED "gene"
########################################################

# First, infer differential expression  over some exposure of interest, which must be a *single* column of "design" matrix X (i.e. db is not currently sufficiently general to handle general factor vectors).

of_in <- "soldier1"
ttT <-
  lmFit(y, X) %>%
  eBayes %>%
  tidy_topTable(of_in = of_in)

# Second, infer differential binding.
# By regression

# For a single Transcription factor binding motif "AR"
ttT %>%
  infer_db(which_tfbms = "AR") %>%
  extract_db
# For all known tfbms: beware multiplicity
ttT %>%
  infer_db %>%
  extract_db

# By Telis:
# Telis requires some heuristic to filter on uncorrected p.value, or logfc (not an inference). this requires that the variable ttt_sub, a scrict subset of the rows of ttt that we have already defined (e.g. the subset whose logfc exceeds some heuristic value.
# The additional options indicate that we have select our prefered tfbm matrix "which_matrix", tfbm hypothesis set "which_tfbms", and methods (p_npar, p_par, which require tt_sub be specified)
ttT %>%
  infer_db(ttT_sub        = filter(ttT, P.Value <= 0.05),
           which_matrix   = exonic1_utr1 ,
           which_tfbms    =  c("ALX3", "ALX4_TBX21", "AR") ,
           n_sim          = 10000) %>%
  extract_db(methods =  c("p_npar", "p_par"))


```
