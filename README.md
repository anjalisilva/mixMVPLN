# mixMVPLN

For clustering three-way count data using mixtures of matrix variate Poisson-log normal distributions.

### Announcements

Version 1.0 was released. 

## Getting started

### Description

Carries out model-based clustering using mixtures of matrix variate Poisson-log normal (MVPLN) model. Markov chain Monte Carlo expectation-maximization algorithm (MCMC-EM) is used for parameter estimation. Information criteria (AIC, BIC, AIC3 and ICL) are offered for model selection. 

### Usage

```R
MVPLNClustering(dataset=testing_dataset, membership=NA, Gmin, Gmax, n_chains=3, n_iterations=NA, init_method="kmeans", n_init_iterations=5)

```
### Arguments

```
dataset             A list of length n, where each entry contain a matrix of counts sized r x p. Assign test dataset using the variable name 'testing_dataset'.
Gmin                Lowest number of components/clusters to be tested. Should be <= Gmax. 
Gmax                Largest number of components/clusters to be tested. Should be >= Gmin. 
n_chains            A positive integer specifying the number of Markov chains. Recommended value is >= 3.  
n_iterations        A positive integer specifying the number of iterations for each chain (including warmup). The warmup is equal to 1/2 of n_iterations. Recommended value is >= number of observations (n).
membership          A vector with length equal to the number of observations (n), indicating the true memberships of each observation. If not available, use NA. 
init_method         Type of initialization strategy to be used. Methods include "kmeans", "random", "medoids", "clara", and "fanny". Default is "kmeans". 
n_init_iterations   Number of runs to be used for the init_method, with a default value of 5. If no initialization, set to 0. 
```

## Details

Output of MVPLNClustering() is an S3 object of class MVPLN. 

Matrix variate distributions offer a natural way for modeling matrices. Extensions of matrix variate distributions in the context of mixture models have given rise to mixtures of matrix variate distributions. The mixture of matrix variate Poisson-log normal (MVPLN) distribution is a multivariate log normal mixture of independent Poisson distributions. The MVPLN distribution can account for both the correlations between variables (p) and the correlations between occasions (r), as two different covariance matrices are used for the two modes. 

The mixtures of MVPLN distribution is used for clustering count data from RNA sequencing using the approach of [Silva et al., 2017](https://arxiv.org/abs/1711.11190v1). The MCMC-EM algorithm via Stan is used for parameter estimation. Coarse grain parallelization is employed, such that when a range of components/clusters (G) are considered, each G is run on a different processor. To check the convergence of MCMC chains, the potential scale reduction factor and the effective number of samples are used. The Heidelberger and Welch’s convergence diagnostic (Heidelberger and Welch, 1983) is used to check the convergence of the MCMC-EM algorithm. 

The AIC, BIC, AIC3 and ICL are used for model selection. Starting values (init_method) and the number of iterations for each chain (n_iterations) play an important role to the successful operation of this algorithm.


## Value

```
dataset                 The input matrix of size n x d.
dimensionality          The number of variables (d). 
normalization_factors   A vector of length d containing the normalization factors for the differences in library sizes.
Gmin                    Lowest number of components/clusters tested.
Gmax                    Largest number of components/clusters tested.
initalization_method    Type of initialization strategy used.
allresults              List of objects of class MPLN giving the results for all models Gmin,..., Gmax.
loglikelihood           Log-likelihood values calculated for each of the fitted models Gmin,..., Gmax.
n_parameters            Number of total parameters for each of the fitted models Gmin,..., Gmax.
truelabels              A vector of length n giving the true labels, if provided. 
ICL.all                 All ICL results, including ICL value, model selected and cluster labels. 
BIC.all                 All BIC results, including BIC value, model selected and cluster labels. 
AIC.all                 All AIC results, including AIC value, model selected and cluster labels. 
AIC3.all                All AIC3 results, including AIC3 value, model selected and cluster labels. 
SlopeHeuristics         All slope heuristics results. 
Djumpmodelselected      The model selected via Djump. 
DDSEmodelselected       The model selected via DDSE.
totaltime               Total time. 
```

## Examples

```R
# Generating data
set.seed(1)

true_G <- 2
true_r <- 2 
true_p <- 3 
true_n <- 100 

# Mu is a r x p matrix
true_M1 <- matrix(rep(6,(true_r*true_p)), ncol=true_p, nrow=true_r, byrow=TRUE)
true_M2 <- matrix(rep(1,(true_r*true_p)), ncol=true_p, nrow=true_r, byrow=TRUE)
true_M_all <- rbind(true_M1, true_M2)

# Phi is a r x r matrix
# Covariance matrix containing variances and covariances btw r varieties
true_Phi1 <- genPositiveDefMat("unifcorrmat",dim=true_r, rangeVar=c(1,1.7))$Sigma
true_Phi1[1,1] <- 1 # for identifiability issues
true_Phi2 = genPositiveDefMat("unifcorrmat",dim=true_r, rangeVar=c(0.7,0.7))$Sigma
true_Phi2[1,1] = 1 # for identifiability issues
true_Phi_all = rbind(true_Phi1,true_Phi2)

# Omega is a p x p matrix 
# Covariance matrix containing variances and covariances btw p growth stages
true_Omega1 <- genPositiveDefMat("unifcorrmat",dim=true_p, rangeVar=c(1,1.7))$Sigma
true_Omega2 <- genPositiveDefMat("unifcorrmat",dim=true_p, rangeVar=c(0.7,0.7))$Sigma
true_Omega_all <- rbind(true_Omega1,true_Omega2)

simulated_counts <- Datagenerator(i=1, r=true_r, p=true_p, n=true_n, pi_g=c(0.79,0.21), mu=true_M_all, phi=true_Phi_all, omega=true_Omega_all)

# Clustering data for G = 1:3

testing_dataset <- simulated_counts # Assign test dataset using the variable name 'testing_dataset'
clustering_results <- MVPLNClustering(dataset=testing_dataset$dataset, membership=simulated_counts$truemembership, Gmin=1, Gmax=3, n_chains=3, n_iterations=300, init_method="kmeans", n_init_iterations=5)

```

## Authors

* Anjali Silva 
* [Sanjeena Subedi](https://sanjeenadang.wordpress.com/)

## License

This project is licensed under the MIT License.

## Acknowledgments

* Dr. Marcelo Ponce, SciNet HPC Consortium, University of Toronto, ON, Canada for all the computational support. 
* This work was funded by Natural Sciences & Engineering Research Council of Canada and Ontario Graduate Fellowship.

## References

[Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2017). A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data. arXiv preprint arXiv:1711.11190.](https://arxiv.org/abs/1711.11190v1)

## Maintainer

* Anjali Silva (anjali@alumni.uoguelph.ca)


