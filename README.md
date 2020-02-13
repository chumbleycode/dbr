
<!-- README.md is generated from README.Rmd. Please edit that file -->
Differential Binding (dbr)
==========================

[![Travis build status](https://travis-ci.org/chumbleycode/dbr.svg?branch=master)](https://travis-ci.org/chumbleycode/dbr) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/chumbleycode/dbr?branch=master&svg=true)](https://ci.appveyor.com/project/chumbleycode/dbr)

The goal of dbr is to implicate a gene *r**e**g**u**l**a**t**o**r* - typically an upstream transcription factor - in the differential RNA expression observed between treatment groups. We then say that there is "differential binding" (DB) of the regulator over treatments. In practice, dbr asks whether the pattern of differential RNA expression over genes reflects (the per-gene count of DNA binding-site motifs for) some upstream gene regulator.

In addition to the raw gene-by-motif count matrices, the package currently provides some functions to augment the popular limma package. The package is in development and is likely to change. It currently reimplements the important TeLiS method of Cole et. al. (2005), which used a gene-set approach to infer DB. Our reimplementation incorporates the most up-to-date motif binding data, offers a non-parametric version of TeLiS (when the sampling frame of genes is small), and provides smooth compatibility with limma.

The package provides simple new functionality that aims to eschews the need to heuristically categorize genes, prior to DB analysis proper, as differentially expressed or not. The simplest - cheap and cheerful - approach is to simply regress gene-specific DE estimates on gene-specific binding-site counts over the entire relevant genome (genes for which it is possible to estimate DE over treatment or exposure groups). This approach will be validated and extended to multilevel modeling.

Installation
------------

You can install dbr with:

``` r
install.packages("devtools")
devtools::install_github("chumbleycode/dbr")
library(dbr)
```

Data: a simple example
----------------------

There are currently three TFBM matrices: utr1, exonic1, exonic\_utr1. Type "utr1" etc into the R console to see these. Get more info for each via ?utr1, ?exonic1, etc. Look for your DNA regulatory motifs of interest in the columns of these matrices. For example, recent literature has examined "a pre-specified set of TFs involved in inflammation (NF-kB and AP-1), IFN response (interferon-stimulated response elements; ISRE), SNS activity (CREB, which mediates SNS-induced b-adrenergic signaling), and glucocorticoid signaling (glucocorticoid receptor; GR)." In biomart nomenclature, "NF-kB" is is identified with NFKB1 or NFKB2. AP-1 is called JUN. ISRE is identified with the set of motifs including IRF2, IRF3, IRF4, 5, 7, 8, 9. CREB is identified with CREB3 or CREB3L1. GR is called NR3C1. This leaves us with 13 regulators plus one complex CEBPG::CREB3L1 (CEBPG\_CREB3L1), as follows. Examine the gene-by-motif count matrices in the R console with:

``` r
immune_tfbms = c("CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
                 "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
                 "NFKB2", "NR3C1")
utr1[, immune_tfbms] # the gene-by-motif matrix for immune motifs of interest
```

Analysis: a simple example.
---------------------------

We examine DB of some immune regulators amoung people with early-life stress (relative to unstressed) using data from Cole et al. (2016). Such analyses generally have two steps.

1.  Differential expression (DE): Estimate differential RNA expression across exposure groups. Here we use a linear model: the exposure must currently be a *single* column of "design" matrix of this linear model (dbr cannot currently handle treatments defined across multiple collumns, e.g. factors with many levels).
2.  Differential binding (DB): Infer dependence of the above, per-gene, estimates on the binding-site count of some regulator(s) of interest.

#### DE

``` r
# Load packages
library(tidyverse)
library(limma)

# Download open source data then specify gene-by-gene regression model
dat = GEOquery::getGEO("GSE77164")[[1]]

# Specify whole-genome regression of rna on design
y <- dat %>% Biobase::exprs()
X <- dat %>%
  Biobase::pData() %>%
  select(age = `age:ch1`,
         soldier = `childsoldier:ch1`,
         edu = `educationlevel:ch1`)
X <- model.matrix(~ soldier + edu + age, data = X) 

# Estimate DE using standard limmma/edger pipeline. 
ttT <-
  lmFit(y, X) %>%
  eBayes %>%
  tidy_topTable(of_in = "soldier1") # "soldier1" is one column of X
```

#### DB

Perhaps the simplest DB analysis is just a regression of gene-wise DE estimates on motif site count. This is an approximation to a full multilevel model.

##### Regression approach

``` r
# regress DE on one motif of interest
summary(lm(logFC ~ NR3C1, data = append_db(ttT))) 

# Or, use dbr to regress logFC on motif site count on all immune_motifs: beware multiple testing
ttT %>%
  infer_db(which_tfbms = immune_tfbms) %>%
  extract_db
```

Here p\_uni is the univariate p-value from a set of simple univariate regressions of logFC on each motif, and p\_cov is the corresponding p-value from a multivariate regression. The latter relates the *p**a**r**t**i**a**l* or conditional relation between a motif and DE, adjusting for the remaining motifs. Any NA's in the output for this column reflect colinearities in the design matrix (i.e. motifs are too highly related to be individually estimated).

##### An alternative approach (see Cole et al)

This approach requires that we first filter some genes to label as categorically DE, e.g. those with a high logFC. This filtering is not, in itself, a statistical inference. We give three examples of how to do this below.

``` r
# 1. genes showing > 20% difference in expression
# (Recalling that logFC is the estimated log2-fold-change of our effect) 
ttT_sub = filter(ttT, logFC >= log2(1.2))

# 2. top and bottom deciles (most extreme 20%)
ttT_sub = filter(ttT, ntile(logFC, 10) %in% c(1,10))
 
# 3. genes whose uncorrected p-values below 0.05 (not an inference):
ttT_sub = filter(ttT, P.Value <= 0.05)
```

Having chosen one of these, or defined your own, the filtered gene-subset enters as the first argument to infer\_db() below, like so:

``` r
ttT %>%
  infer_db(ttT_sub        = ttT_sub,
           which_tfbms    = immune_tfbms) %>%
  extract_db
```

Here "p\_par" is the 2 sided p-value for parametric TeLiS (par\_p\_over and par\_p\_under are the corresponding one-tailed values for over and motif-underrepresentation). If perm\_telis = TRUE, then p\_npar will give a (computationally costly) permutation p-value for TeLiS.

##### Chi-squared test

Other ways to relate DE labels to transcription factor motif, e.g. CREB3 motif. Having categorized genes as "DE" or not, we can examine the relation between this label and the motif count as follows, for example.

``` r
# Chi-squared
append_db(ttT,ttT_sub = ttT_sub) %>% 
  select(gene_subset, CREB3) %>% 
  table %>% 
  chisq.test
```

##### Some other ways to filter genes for TeLiS:

"In addition to a priori gene sets, we also assessed empirical transcriptome correlates of each demographic factor (Figure 2A). Each demographic parameter was associated with hundreds of genes showing &gt; 20% difference in expression across the observed range of variation (although the statistical significance of these individual transcript associations varied, with some dimensions such sex, race, and BMI showing robust effects at a genome-wide false discovery rate (FDR) of 5%, whereas others failed to yield significant effects after correction for multiple testing; Supporting Information Figure S1)."
