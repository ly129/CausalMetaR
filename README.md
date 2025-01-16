
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CausalMetaR

<!-- badges: start -->
<!-- badges: end -->

The `CausalMetaR` package provides robust and efficient methods for
estimating causal effects in a target population using a multi-source
dataset. The multi-source data can be a collection of trials,
observational studies, or a combination of both, which have the same
data structure (outcome, treatment, and covariates). The target
population can be based on an internal dataset or an external dataset
where only covariate information is available. The causal estimands
available are average treatment effects (ATEs) and subgroup treatment
effects (STEs).

## Installation

You can install the released version of `CausalMetaR` from
[CRAN](https://CRAN.R-project.org/package=CausalMetaR) with:

``` r
install.packages("CausalMetaR")
```

You can install the development version of `CausalMetaR` from
[GitHub](https://github.com/ly129/CausalMetaR) with:

``` r
# install.packages("devtools")
devtools::install_github("ly129/CausalMetaR")
```

## Usage

Refer to [Wang et al. (2024)](https://doi.org/10.48550/arXiv.2402.04341)
for a detailed guide on using the package. Below, we include the code
used in the examples in the paper. The examples are based on the
multi-source dataset `dat_multisource` included in the package.

Loading the package:

``` r
library(CausalMetaR)
```

Specifying the working models:

``` r
outcome_model_args <- list(family = gaussian(),
                           SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"))
treatment_model_args <- list(family = binomial(),
                             SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"))
external_model_args = list(family = binomial(),
                           SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"))
```

Setting a random number seed for reproducibility:

``` r
set.seed(1234)
```

The examples below estimate ATEs and STEs in external and internal
target populations. Each example may take a few minutes to run on
standard laptop.

### Example 1: Estimating the ATE in the external target population

``` r
result_ae <- ATE_external(
  Y = dat_multisource$Y,
  S = dat_multisource$S,
  A = dat_multisource$A,
  X = dat_multisource[, 1:10],
  X_external = dat_external[, 1:10],
  outcome_model_args = outcome_model_args,
  treatment_model_args = treatment_model_args,
  external_model_args = external_model_args,
  cross_fitting = TRUE,
  replications = 5)
#> Loading required package: nnls
#> Loading required namespace: nnet
result_ae
#> AVERAGE TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION
#> 
#> Treatment effect (mean difference) estimates:
#> ---------------------------------------------
#>  Estimate     SE Lower 95% CI Upper 95% CI
#>    6.6294 0.1535       5.8616       7.3972
```

### Example 2: Estimating ATEs in the internal target populations

``` r
result_ai <- ATE_internal(
  Y = dat_multisource$Y,
  S = dat_multisource$S,
  A = dat_multisource$A,
  X = dat_multisource[, 1:10],
  outcome_model_args = outcome_model_args,
  treatment_model_args = treatment_model_args,
  cross_fitting = TRUE,
  replications = 5)
result_ai
#> AVERAGE TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS
#> 
#> Treatment effect (mean difference) estimates:
#> ---------------------------------------------
#>  Source Estimate     SE Lower 95% CI Upper 95% CI
#>       A   6.5874 0.1903       6.2145       6.9603
#>       B   7.7556 0.2577       7.2506       8.2606
#>       C   7.2916 0.3594       6.5872       7.9960
```

### Example 3: Estimating STEs in the external target population

``` r
result_se <- STE_external(
  Y = dat_multisource$Y, 
  S = dat_multisource$S, 
  A = dat_multisource$A, 
  X = dat_multisource[, 2:10], 
  EM = dat_multisource$EM,
  X_external = dat_external[, 2:10], 
  EM_external = dat_external$EM, 
  outcome_model_args = outcome_model_args, 
  treatment_model_args = treatment_model_args, 
  external_model_args = external_model_args,
  cross_fitting = TRUE,
  replications = 5)
result_se
#> SUBGROUP TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION
#> 
#> Treatment effect (mean difference) estimates:
#> ---------------------------------------------
#>  Subgroup Estimate     SE Lower 95% CI Upper 95% CI Lower 95% SCB Upper 95% SCB
#>         a   7.0787 0.3563       5.9088       8.2485        5.5453        8.6121
#>         b   5.5207 0.2321       4.5764       6.4650        4.2830        6.7585
#>         c   7.5709 0.1805       6.7382       8.4037        6.4794        8.6625
#>         d   6.5748 0.2253       5.6446       7.5051        5.3556        7.7941
#>         e   5.3741 0.3382       4.2343       6.5139        3.8802        6.8681
```

### Example 4: Estimating STEs in the internal target populations

``` r
result_si <- STE_internal(
  Y = dat_multisource$Y, 
  S = dat_multisource$S, 
  A = dat_multisource$A, 
  X = dat_multisource[, 2:10], 
  EM = dat_multisource$EM,
  outcome_model_args = outcome_model_args, 
  treatment_model_args = treatment_model_args,
  cross_fitting = TRUE,
  replications = 5)
result_si
#> SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS
#> 
#> Treatment effect (mean difference) estimates:
#> ---------------------------------------------
#>  Source Subgroup Estimate     SE Lower 95% CI Upper 95% CI Lower 95% SCB
#>       A        a   6.9197 0.5001       5.9395       7.8999        5.6345
#>                b   5.4340 0.3681       4.7126       6.1555        4.4880
#>                c   7.5452 0.3097       6.9383       8.1522        6.7493
#>                d   6.5053 0.3630       5.7939       7.2168        5.5724
#>                e   5.4595 0.5215       4.4373       6.4816        4.1192
#>       B        a   8.2134 0.7554       6.7328       9.6939        6.2720
#>                b   6.8396 0.5381       5.7850       7.8942        5.4567
#>                c   8.7995 0.4232       7.9699       9.6290        7.7118
#>                d   7.7098 0.4529       6.8223       8.5974        6.5460
#>                e   6.4405 0.6975       5.0734       7.8076        4.6479
#>       C        a   8.0544 1.2049       5.6929      10.4159        4.9579
#>                b   6.2151 0.6711       4.8998       7.5304        4.4904
#>                c   8.2620 0.6426       7.0026       9.5215        6.6106
#>                d   7.1427 0.6538       5.8612       8.4242        5.4624
#>                e   6.0910 0.8718       4.3823       7.7997        3.8504
#>  Upper 95% SCB
#>         8.2050
#>         6.3801
#>         8.3411
#>         7.4382
#>         6.7998
#>        10.1547
#>         8.2225
#>         9.8872
#>         8.8737
#>         8.2332
#>        11.1509
#>         7.9398
#>         9.9135
#>         8.8231
#>         8.3316
```

## Reference

To cite this package, use:

> Wang G, McGrath S, Lian Y. CausalMetaR: An R package for performing
> causally interpretable meta-analyses. *Research Synthesis Methods*. In
> press.
