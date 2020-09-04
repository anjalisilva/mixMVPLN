# `mixMVPLN`

## Description
`mixMVPLN` is an R package for performing model-based clustering of three-way count data using mixtures of matrix variate Poisson-log normal (MVPLN) distributions proposed by [Silva et al., 2018](https://arxiv.org/abs/1807.08380).

Main function, __*mvplnClustering*__, can be used to carry out model-based clustering using mixtures of MVPLN model. Markov chain Monte Carlo expectation-maximization algorithm (MCMC-EM) is used for parameter estimation. Information criteria (AIC, BIC, AIC3 and ICL) are offered for model selection. Function __*mvplnDataGenerator*__ is available to generate simlulation data. For more information, see details section below.

## Installation

To install the latest version of the package:

``` r
require("devtools")
install_github("anjalisilva/mixMVPLN", build_vignettes = TRUE)
library("mixMVPLN")
```

## Overview

`mixMVPLN` contains 2 functions. For the purpose of generating simulation data via mixtures of MVPLN: *mvplnDataGenerator*. For carrying out clustering of count data using mixtures of MVPLN with parallelization: *mvplnClustering*. 

To list all functions available in the package:

``` r
ls("package:mixMVPLN")
```

## Details

Matrix variate distributions offer a natural way for modeling matrices. Extensions of matrix variate distributions in the context of mixture models have given rise to mixtures of matrix variate distributions.

The multivariate Poisson-log normal (MPLN) distribution was proposed in 1989 ([Aitchison and Ho, 1989](https://www.jstor.org/stable/2336624?seq=1)). A multivariate Poisson-log normal mixture model for clustering of count data was proposed by [Silva et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31311497/). Here this work is extended and a mixture of matrix variate Poisson-log normal (MVPLN) distribution for clustering three-way count data is proposed by [Silva et al., 2018](https://arxiv.org/abs/1807.08380). A mixture of MVPLN distribution is a multivariate log normal mixture of independent Poisson distributions. The MVPLN distribution can account for both the correlations between variables (p) and the correlations between occasions (r), as two different covariance matrices are used for the two modes. 

The MCMC-EM algorithm via Stan is used for parameter estimation. Coarse grain parallelization is employed, such that when a range of components/clusters (g = 1,...,G) are considered, each component/cluster is run on a different processor. This can be performed because each component/cluster size is independent from another. All components/clusters in the range to be tested have been parallelized to run on a seperate core using the *parallel* R package. The number of cores used for clustering is can be user-specified or calculated using *parallel::detectCores() - 1*. To check the convergence of MCMC chains, the potential scale reduction factor and the effective number of samples are used. The Heidelberger and Welch’s convergence diagnostic (Heidelberger and Welch, 1983) is used to check the convergence of the MCMC-EM algorithm. 

The AIC, BIC, AIC3 and ICL are used for model selection. Starting values (argument: initMethod) and the number of iterations for each chain (argument: nInitIterations) play an important role to the successful operation of this algorithm.

  
## Tutorials  

For tutorials, refer to the vignette:

``` r
browseVignettes("mixMVPLN")
```

or see [A tour of mixMVPLN](https://github.com/anjalisilva/mixMVPLN/blob/master/vignettes/Introduction_mixMVPLN.md).
  


## Citation for Package
``` r
citation("mixMVPLN")
```
Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2018). Finite mixtures of matrix-variate Poisson-log normal distributions for three-way count data. arXiv preprint arXiv:1807.08380.
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


## References for Package

* [Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2018). Finite mixtures of matrix-variate Poisson-log normal distributions for three-way count data. arXiv preprint arXiv:1807.08380.](https://arxiv.org/abs/1807.08380)

* [Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data. *BMC Bioinformatics.*](https://pubmed.ncbi.nlm.nih.gov/31311497/)

* [Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution. *Biometrika.*](https://www.jstor.org/stable/2336624?seq=1)

* For others, refer to help page of inidividual functions via `?function` or `help(function)`.


## Maintainer

* Anjali Silva (anjali.silva@uhnresearch.ca). 


## Contributions

`mixMVPLN` welcomes issues, enhancement requests, and other contributions. To submit an issue, use the [GitHub issues](https://github.com/anjalisilva/mixMVPLN/issues).


## Acknowledgments

* This work was funded by Natural Sciences and Engineering Research Council of Canada, Queen Elizabeth II Graduate Scholarship, and Arthur Richmond Memorial Scholarship. This work was conducted at the University of Guelph, Ontario, Canada. 
