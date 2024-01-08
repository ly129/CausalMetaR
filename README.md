
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMetafoR

<!-- badges: start -->
<!-- badges: end -->

The `CMetafoR` package provides robust and efficient methods for
estimating causal effects in a target population using a multi-source
dataset. The multi-source data can be a collection of trials,
observational studies, or a combination of both, which have the same
data structure (outcome, treatment, and covariates). The target
population can be based on an internal dataset or an external dataset where
only covariate information is available. The causal estimands available
are average treatment effects and subgroup treatment effects.

## Installation

You can install the development version of `CMetafoR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ly129/CMetafoR")
```
