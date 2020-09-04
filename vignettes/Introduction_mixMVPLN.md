
# A quick tour of mixMVPLN

## Introduction

[**mixMVPLN**](https://arxiv.org/abs/1807.08380) is an R package for model-based clustering based on mixtures of matrix variate Poisson-log normal (MVPLN) distributions. It is applicable for clustering of three-way count data. [**mixMVPLN**](https://arxiv.org/abs/1807.08380) provides functions for parameter estimation via the Markov chain Monte Carlo expectation-maximization (MCMC-EM) algorithm. Information criteria (AIC, BIC, AIC3 and ICL) are offered for model selection. Also included is a function for simulating data from this model. Function *mvplnClustering* within [**mixMVPLN**](https://arxiv.org/abs/1807.08380) makes use of the **parallel** R package to run each component/cluster (G) in parallel, as each G is independent from another. 

This document gives a quick tour of [**mixMVPLN**](https://arxiv.org/abs/1807.08380) (version 0.1.0) functionalities. It was written in R Markdown, using the **knitr** package for production. See `help(package = "mixMVPLN")` for further details and references provided by `citation("mixMVPLN")`. To download [**mixMVPLN**](https://arxiv.org/abs/1807.08380), use the following commands:

``` r
require("devtools")
install_github("anjalisilva/mixMVPLN", build_vignettes = TRUE)
library("mixMVPLN")
```

<br>


## Data Simulation

The function *mvplnDataGenerator* permits to simulate data from a mixture of MVPLN distributions. See *?mvplnDataGenerator* for more information, an example, and references. To simulate a dataset from a mixture of MVPLN with 50 units, 3 total repsonses, and 2 occasions, with two mixture components, each with a mixing proportion of 0.79 and 0.21, respectively, let us use *mvplnDataGenerator*. Here **clusterGeneration** R package is used to generate positive definite covariance matrices for illustration purposes. 

``` r
set.seed(1234) # for reproducibility, setting seed
true_G <- 2 # number of total components/clusters
true_r <- 2 # number of total occasions
true_p <- 3 # number of total responses
true_n <- 50 # number of total units


# Mu is a r x p matrix
true_M1 <- matrix(rep(6, (true_r * true_p)),
                  ncol = true_p,
                  nrow = true_r, byrow = TRUE)

true_M2 <- matrix(rep(1, (true_r * true_p)),
                  ncol = true_p,
                  nrow = true_r,
                  byrow = TRUE)

true_M_all <- rbind(true_M1, true_M2)

# Phi is a r x r matrix
# Loading needed packages for generating data
if (!require(clusterGeneration)) install.packages("clusterGeneration")
# Covariance matrix containing variances and covariances between r occasions
true_Phi1 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
                                                  dim = true_r,
                                                  rangeVar = c(1, 1.7))$Sigma
true_Phi1[1, 1] <- 1 # For identifiability issues

true_Phi2 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
                                                  dim = true_r,
                                                  rangeVar = c(0.7, 0.7))$Sigma
true_Phi2[1, 1] <- 1 # For identifiability issues
true_Phi_all <- rbind(true_Phi1, true_Phi2)

# Omega is a p x p matrix
# Covariance matrix containing variances and covariances between p responses
true_Omega1 <- genPositiveDefMat("unifcorrmat", dim = true_p,
                                 rangeVar = c(1, 1.7))$Sigma

true_Omega2 <- genPositiveDefMat("unifcorrmat", dim = true_p,
                                 rangeVar = c(0.7, 0.7))$Sigma

true_Omega_all <- rbind(true_Omega1,true_Omega2)

# Simulating data 
simulatedMVData <- mvplnDataGenerator(nOccasions = true_r,
                                 nResponses = true_p,
                                 nUnits = true_n,
                                 mixingProportions = c(0.79, 0.21),
                                 matrixMean = true_M_all,
                                 phi = true_Phi_all,
                                 omega = true_Omega_all)
```
<br>

The generated dataset can be checked:

``` r
length(simulatedMVData$dataset) # 50 units
class(simulatedMVData$dataset) # list with length of 50
typeof(simulatedMVData$dataset) # list
summary(simulatedMVData$dataset) # summary of data
dim(simulatedMVData$dataset[[1]]) # dimension of first unit is 2 x 3
                                  # 2 occasions and 3 responses
```

<br>

<div style="text-align:left">

## Clustering

<div style="text-align:left">
Once the count data is available, clustering can be performed using the *mvplnClustering* function. See *?mvplnClustering* for more information, an example, and references. Here, clustering will be performed using the above generated dataset. 

Coarse grain parallelization is employed in *mvplnClustering*, such that when a range of components/clusters (g = 1,...,G) are considered, each component/cluster size is run on a different processor. This can be performed because each component/cluster size is independent from another. All components/clusters in the range to be tested have been parallelized to run on a seperate core using the *parallel* R package. The number of cores used for clustering can be specified by user. Otherwise, it is internally determined using *parallel::detectCores() - 1*.

Below, clustering of *simulatedMVData$dataset* is performed for g = 1:2 with 300 iterations and *kmeans* initialization with 2 initialization runs. 

``` r
clusteringResults <- mvplnClustering(dataset = simulatedMVData$dataset,
                                     membership = simulatedMVData$truemembership,
                                     gmin = 1,
                                     gmax = 2,
                                     nChains = 3,
                                     nIterations = 300,
                                     initMethod = "kmeans",
                                     nInitIterations = 2,
                                     normalize = "Yes")
```

The model selected by BIC for this dataset can be viewed as follows.

``` r
clusteringResults$BIC.all$BICmodelselected

# Cross tabulation of BIC selected model labels with true lables
table(mplnResults$BIC.all$BICmodelselected_labels, simulatedMVData$truemembership)
```

<br>

<div style="text-align:left">

## References

[Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2018). Finite mixtures of matrix-variate Poisson-log normal distributions for three-way count data. arXiv preprint arXiv:1807.08380.](https://arxiv.org/abs/1807.08380)



