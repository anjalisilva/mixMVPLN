
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixMVPLN

Mixtures of Matrix Variate Poisson-log Normal Model for Clustering
Three-way Count Data

<!-- badges: start -->
<!-- https://www.codefactor.io/repository/github/anjalisilva/mixmvpln/issues -->
<!-- [![CodeFactor](https://www.codefactor.io/repository/github/anjalisilva/mixmvpln/badge)](https://www.codefactor.io/repository/github/anjalisilva/mixmvpln)-->

[![GitHub
issues](https://img.shields.io/github/issues/anjalisilva/mixMVPLN)](https://github.com/anjalisilva/mixMVPLN/issues)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
![GitHub language
count](https://img.shields.io/github/languages/count/anjalisilva/mixMVPLN)
![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/anjalisilva/mixMVPLN/master)

<!-- https://shields.io/category/license -->
<!-- badges: end -->

## Description

`mixMVPLN` is an R package for performing model-based clustering of
three-way count data using mixtures of matrix variate Poisson-log normal
(MVPLN) distributions ([Silva et al.,
2018](https://arxiv.org/abs/1807.08380)). Three different frameworks are
available for parameter estimation of the mixtures of MVPLN models: 1)
method based on Markov chain Monte Carlo expectation-maximization
algorithm (MCMC-EM), 2) method based on variational Gaussian
approximations (VGAs), and 3) a hybrid approach that combines both the
variational approximation-based approach and MCMC-EM-based approach.
Information criteria (AIC, BIC, AIC3 and ICL) are offered for model
selection. Also included is a function for simulating data from MVPLN
model.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("anjalisilva/mixMVPLN", build_vignettes = TRUE)
library("mixMVPLN")
```

## Overview

To list all functions available in the package:

``` r
ls("package:mixMVPLN")
```

`mixMVPLN` contains 5 functions. For the purpose of generating
simulation data via mixtures of MVPLN: *mvplnDataGenerator*. For
carrying out clustering of count data using mixtures of MVPLN via 1)
method based on MCMC-EM with parallelization: *mvplnMCMCclus*; 2) method
based on VGAs: *mvplnVGAclus*; and 3) the hybrid approach that combines
both the VGAs and MCMC-EM-based approach: *mvplnHybriDclus*. For
visualization of clustering results, there is the: *mvplnVisualize*.

An overview of the package is illustrated below:

<div style="text-align:center">

<img src="inst/extdata/Overview_mixMVPLN.png" width="800" height="450"/>

<div style="text-align:left">

<div style="text-align:left">

<div style="text-align:left">


## Details

Matrix variate distributions offer a natural way for modeling matrices.
Extensions of matrix variate distributions in the context of mixture
models have given rise to mixtures of matrix variate distributions.

The multivariate Poisson-log normal (MPLN) distribution was proposed in
1989 ([Aitchison and Ho,
1989](https://www.jstor.org/stable/2336624?seq=1)). A multivariate
Poisson-log normal mixture model for clustering of count data was
proposed by [Silva et al.,
2019](https://pubmed.ncbi.nlm.nih.gov/31311497/). Here this work is
extended and a mixture of matrix variate Poisson-log normal (MVPLN)
distribution for clustering three-way count data is proposed by [Silva
et al., 2018](https://arxiv.org/abs/1807.08380). A mixture of MVPLN
distribution is a multivariate log normal mixture of independent Poisson
distributions. The MVPLN distribution can account for both the
correlations between variables (p) and the correlations between
occasions (r), as two different covariance matrices are used for the two
modes.

The MCMC-EM algorithm via Stan is used for parameter estimation. Coarse
grain parallelization is employed, such that when a range of
components/clusters (g = 1,…,G) are considered, each component/cluster
is run on a different processor. This can be performed because each
component/cluster size is independent from another. All
components/clusters in the range to be tested have been parallelized to
run on a seperate core using the *parallel* R package. The number of
cores used for clustering is can be user-specified or calculated using
*parallel::detectCores() - 1*. To check the convergence of MCMC chains,
the potential scale reduction factor and the effective number of samples
are used. The Heidelberger and Welch’s convergence diagnostic
(Heidelberger and Welch, 1983) is used to check the convergence of the
MCMC-EM algorithm.

The AIC, BIC, AIC3 and ICL are used for model selection. Starting values
(argument: initMethod) and the number of iterations for each chain
(argument: nInitIterations) play an important role to the successful
operation of this algorithm.

## Tutorials

For tutorials, refer to the vignette:

``` r
browseVignettes("mixMVPLN")
```

or see [A tour of
mixMVPLN](https://github.com/anjalisilva/mixMVPLN/blob/master/vignettes/Introduction_mixMVPLN_MCMEM.md).

## Citation for Package

``` r
citation("mixMVPLN")
```

Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2018).
Finite mixtures of matrix-variate Poisson-log normal distributions for
three-way count data. arXiv preprint arXiv:1807.08380.

``` r
A BibTeX entry for LaTeX users is

  @misc{,
    title = {Finite mixtures of matrix-variate Poisson-log normal distributions for three-way count data},
    author = {A. Silva and S. J. Rothstein and P. D. McNicholas and S. Subedi},
    note = {arXiv preprint arXiv:1206.4768},
    year = {2018},
    url = {https://arxiv.org/abs/1807.08380},
  }
```

## Package References

[Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2018).
Finite mixtures of matrix-variate Poisson-log normal distributions for
three-way count data. arXiv preprint
arXiv:1807.08380.](https://arxiv.org/abs/1807.08380)

## Other References

[Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal
distribution. *Biometrika.*](https://www.jstor.org/stable/2336624?seq=1)

[Akaike, H. (1973). Information theory and an extension of the maximum
likelihood principle. In *Second International Symposium on Information
Theory*, New York, NY, USA, pp. 267–281. Springer
Verlag.](https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15)

[Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture
model for clustering with the integrated classification likelihood.
*IEEE Transactions on Pattern Analysis and Machine Intelligence*
22.](https://hal.inria.fr/inria-00073163/document)

[Bozdogan, H. (1994). Mixture-model cluster analysis using model
selection criteria and a new informational measure of complexity. In
*Proceedings of the First US/Japan Conference on the Frontiers of
Statistical Modeling: An Informational Approach: Volume 2 Multivariate
Statistical Modeling*, pp. 69–113. Dordrecht: Springer
Netherlands.](https://link.springer.com/chapter/10.1007/978-94-011-0800-3_3)

[Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method
for differential expression analysis of RNA-seq data. *Genome Biology*
11,
R25.](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)

[Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics* 6.](https://www.jstor.org/stable/2958889?seq=1)

[Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A
multivariate Poisson-log normal mixture model for clustering
transcriptome sequencing data. *BMC
Bioinformatics.*](https://pubmed.ncbi.nlm.nih.gov/31311497/)

## Maintainer

- Anjali Silva (<anjali@alumni.uoguelph.ca>).

## Contributions

`mixMVPLN` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/mixMVPLN/issues).

## Acknowledgments

- This work was initially started at University of Guelph, Ontario,
  Canada and was funded by Ontario Graduate Fellowship (Silva), Queen
  Elizabeth II Graduate Scholarship (Silva), Arthur Richmond Memorial
  Scholarship (Silva), and Discovery Grant from Natural Sciences and
  Engineering Research Council of Canada (Dang).

- Later work of Silva was conducted at the Princess Margaret Cancer
  Centre - University Health Network, Ontario, Canada and was supported
  by CIHR Postdoctoral Fellowship and resources received for
  Postgraduate Affiliation Program of Vector Institute. Later work of
  Dang was conducted at Carleton University, Ontario Canada and was
  supported by Canada Research Chair Program.

- We acknowledge Steven J. Rothstein (UGuelph), Paul D. McNicholas
  (McMasterU), Xiaoke Qin (CarletonU) and Marcelo Ponce (UToronto) for
  their input.
