% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvplnVGAclus.R
\name{mvplnVGAclus}
\alias{mvplnVGAclus}
\title{Clustering Using mixtures of MVPLN via variational Gaussian approximations}
\usage{
mvplnVGAclus(
  dataset,
  membership = "none",
  gmin,
  gmax,
  initMethod = "kmeans",
  nInitIterations = 0,
  normalize = "Yes"
)
}
\arguments{
\item{dataset}{A list of length nUnits, containing Y_j matrices.
A matrix Y_j has size r x p, and the dataset will have 'j' such
matrices with j = 1,...,n. If a Y_j has all zeros, such Y_j will
be removed prior to cluster analysis.}

\item{membership}{A numeric vector of size length(dataset) containing
the cluster membership of each Y_j matrix. If not available, leave
as "none".}

\item{gmin}{A positive integer specifying the minimum number of
components to be considered in the clustering run.}

\item{gmax}{A positive integer, >gmin, specifying the maximum number of
components to be considered in the clustering run.}

\item{initMethod}{A method for initialization. Current options are
"kmeans", "random", "medoids", "clara", or "fanny". If nInitIterations
is set to >= 1, then this method will be used for initialization.
Default is "kmeans".}

\item{nInitIterations}{A positive integer or zero, specifying the number
of initialization runs prior to clustering. The run with the highest
loglikelihood is used to initialize values, if nInitIterations > 1.
Default value is set to 0, in which case no initialization will
be performed.}

\item{normalize}{A character string "Yes" or "No" specifying
if normalization should be performed. Currently, normalization
factors are calculated using TMM method of edgeR package.
Default is "Yes", for which normalization will be performed.}
}
\value{
Returns an S3 object of class mvplnParallel with results.
\itemize{
  \item dataset - The input dataset on which clustering is performed.
  \item nUnits - Number of units in the input dataset.
  \item nVariables - Number of variables in the input dataset.
  \item nOccassions - Number of occassions in the input dataset.
  \item normFactors - A vector of normalization factors used for
     input dataset.
  \item gmin - Minimum number of components considered in the
     clustering run.
  \item gmax - Maximum number of components considered in the
     clustering run.
  \item initalizationMethod - Method used for initialization.
  \item allResults - A list with all results.
  \item loglikelihood - A vector with value of final log-likelihoods
     for each cluster size.
  \item nParameters - A vector with number of parameters for each
     cluster size.
  \item trueLabels - The vector of true labels, if provided by user.
  \item ICL_all - A list with all ICL model selection results.
  \item BIC_all - A list with all BIC model selection results.
  \item AIC_all - A list with all AIC model selection results.
  \item AIC3_all - A list with all AIC3 model selection results.
  \item totalTime - Total time used for clustering and model selection.
}
}
\description{
Performs clustering using mixtures of matrix variate Poisson-log normal
(MVPLN) via variational Gaussian approximations. Model selection can
be done using AIC, AIC3, BIC and ICL.
}
\details{
Starting values (argument: initMethod) and the number of
    iterations for each MCMC chain (argument: nIterations) play an
    important role in the successful operation of this algorithm.
    Occasionally an error may result due to singular matrix. In that
    case rerun the method and may set a different seed to ensure
    a different set of initialization values.
}
\examples{
\dontrun{
# Example 1
set.seed(1234)
trueG <- 2 # number of total G
truer <- 2 # number of total occasions
truep <- 3 # number of total responses
trueN <- 1000 # number of total units

# Mu is a r x p matrix
trueM1 <- matrix(rep(6, (truer * truep)),
                 ncol = truep,
                 nrow = truer, byrow = TRUE)

trueM2 <- matrix(rep(1, (truer * truep)),
                 ncol = truep,
                 nrow = truer,
                 byrow = TRUE)

trueMall <- rbind(trueM1, trueM2)

# Phi is a r x r matrix
# Loading needed packages for generating data
# if (!require(clusterGeneration)) install.packages("clusterGeneration")
# library("clusterGeneration")

# Covariance matrix containing variances and covariances between r occasions
# truePhi1 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#                                                   dim = truer,
#                                                   rangeVar = c(1, 1.7))$Sigma
truePhi1 <- matrix(c(1.075551, -0.488301, -0.488301, 1.362777), nrow = 2)
truePhi1[1, 1] <- 1 # For identifiability issues

# truePhi2 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#                                                   dim = truer,
#                                                   rangeVar = c(0.7, 0.7))$Sigma
truePhi2 <- matrix(c(0.7000000, 0.6585887, 0.6585887, 0.7000000), nrow = 2)
truePhi2[1, 1] <- 1 # For identifiability issues
truePhiall <- rbind(truePhi1, truePhi2)

# Omega is a p x p matrix
# Covariance matrix containing variances and covariances between p responses
# trueOmega1 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#                                    rangeVar = c(1, 1.7))$Sigma
trueOmega1 <- matrix(c(1.0526554, 1.0841910, -0.7976842,
                       1.0841910,  1.1518811, -0.8068102,
                       -0.7976842, -0.8068102,  1.4090578),
                       nrow = 3)
# trueOmega2 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#                                    rangeVar = c(0.7, 0.7))$Sigma
trueOmega2 <- matrix(c(0.7000000, 0.5513744, 0.4441598,
                       0.5513744, 0.7000000, 0.4726577,
                       0.4441598, 0.4726577, 0.7000000),
                       nrow = 3)
trueOmegaAll <- rbind(trueOmega1, trueOmega2)

# Generated simulated data
sampleData <- mixMVPLN::mvplnDataGenerator(nOccasions = truer,
                                           nResponses = truep,
                                           nUnits = trueN,
                                           mixingProportions = c(0.79, 0.21),
                                           matrixMean = trueMall,
                                           phi = truePhiall,
                                           omega = trueOmegaAll)

# Clustering simulated matrix variate count data
clusteringResults <- mixMVPLN::mvplnVGAclus(
                      dataset = sampleData$dataset,
                      membership = sampleData$truemembership,
                      gmin = 1,
                      gmax = 3,
                      initMethod = "kmeans",
                      nInitIterations = 3,
                      normalize = "Yes")

# Example 2
trueG <- 1 # number of total G
truer <- 2 # number of total occasions
truep <- 3 # number of total responses
trueN <- 1000 # number of total units
truePiG <- 1L # mixing proportion for G = 1

# Mu is a r x p matrix
trueM1 <- matrix(c(6, 5.5, 6, 6, 5.5, 6),
                 ncol = truep,
                 nrow = truer,
                 byrow = TRUE)
trueMall <- rbind(trueM1)
# Phi is a r x r matrix
# truePhi1 <- clusterGeneration::genPositiveDefMat(
#                               "unifcorrmat",
#                                dim = truer,
#                               rangeVar = c(0.7, 1.7))$Sigma
truePhi1 <- matrix(c(1.3092747, 0.3219674,
                     0.3219674, 1.3233794), nrow = 2)
truePhi1[1, 1] <- 1 # for identifiability issues
truePhiall <- rbind(truePhi1)

# Omega is a p x p matrix
# trueOmega1 <- genPositiveDefMat(
#                    "unifcorrmat",
#                     dim = truep,
#                     rangeVar = c(1, 1.7))$Sigma
trueOmega1 <- matrix(c(1.1625581, 0.9157741, 0.8203499,
                       0.9157741, 1.2216287, 0.7108193,
                       0.8203499, 0.7108193, 1.2118854), nrow = 3)
trueOmegaAll <- rbind(trueOmega1)

# Generated simulated data
sampleData2 <- mixMVPLN::mvplnDataGenerator(
                         nOccasions = truer,
                         nResponses = truep,
                         nUnits = trueN,
                         mixingProportions = truePiG,
                         matrixMean = trueMall,
                         phi = truePhiall,
                         omega = trueOmegaAll)

# Clustering simulated matrix variate count data
clusteringResults <- mixMVPLN::mvplnVGAclus(
                      dataset = sampleData2$dataset,
                      membership = sampleData2$truemembership,
                      gmin = 1,
                      gmax = 2,
                      initMethod = "clara",
                      nInitIterations = 2,
                      normalize = "Yes")
# Example 3
# mixture of six independent Poisson distributions with G = 2
set.seed(23456)
totalSet <- 1 # Number of datasets
data <- list()
for (i in 1:totalSet) {
 N <- 1000 # biological samples e.g. genes
 d <- 6 # dimensionality e.g. conditions*replicates = total samples
 G <- 2

 piG <- c(0.45, 0.55) # mixing proportions

 tmu1 <- c(1000, 1500, 1500, 1000, 1500, 1000)
 tmu2 <- c(1000, 1000, 1000, 1500, 1000, 1200)
 z <- t(rmultinom(n = N, size = 1, prob = piG))
 para <- list()
 para[[1]] <- NULL # for pi
 para[[2]] <- list() # for mu
 para[[3]] <- list() # for Z
 names(para) <- c("pi", "mu", "Z")

 para[[1]] <- piG
 para[[2]] <- list(tmu1, tmu2)
 para[[3]] <- z

 Y <- matrix(0, ncol = d, nrow = N)

 for (g in 1:G) {
  obs <- which(para[[3]][, g] == 1)
  for (j in 1:d) {
   Y[obs, j] <- rpois(length(obs), para[[2]][[g]][j])
  }
 }

 Y_mat <- list()
 for (obs in 1:N) {
  Y_mat[[obs]] <- matrix(Y[obs, ], nrow = 2, byrow = TRUE)
 }

 data[[i]] <- list()
 data[[i]][[1]] <- para
 data[[i]][[2]] <- Y
 data[[i]][[3]] <- Y_mat
}
outputExample3 <- mixMVPLN::mvplnVGAclus(dataset = data[[1]][[3]],
                               membership = "none",
                               gmin = 1,
                               gmax = 3,
                               initMethod = "kmeans",
                               nInitIterations = 1,
                               normalize = "Yes")

# Example 4
# mixture of six independent negative binomial distributions with G = 2
set.seed(23456)
N1 <- 2000 # biological samples e.g. genes
d1 <- 6 # dimensionality e.g. conditions*replicates = total samples
piG1 <- c(0.79, 0.21) # mixing proportions for G=2
trueMu1 <- rep(c(1000, 500), 3)
trueMu2 <- sample(c(1000, 500), 6, replace=TRUE)
trueMu <- list(trueMu1, trueMu2)

z <- t(rmultinom(N1, size = 1, piG1))

nG <- vector("list", length = length(piG1))
dataNB <- matrix(NA, ncol = d1, nrow = N1)

for (i in 1:length(piG1)) { # This code is not generalized, for G = 2 only
  nG[[i]] <- which(z[, i] == 1)
  for (j in 1:d1) {
    dataNB[nG[[i]], j] <- rnbinom(n = length(nG[[i]]),
                                  mu = trueMu[[i]][j],
                                  size = 100)
  }
}

dataNBmat <- list()
for (obs in 1:N1) {
  dataNBmat[[obs]] <- matrix(dataNB[obs,],
                             nrow = 2,
                             byrow = TRUE)
}

outputExample4 <- mixMVPLN::mvplnVGAclus(dataset= dataNBmat,
                                  membership = "none",
                                  gmin = 1,
                                  gmax = 2,
                                  initMethod = "kmeans",
                                  nInitIterations = 1,
                                  normalize = "Yes")
}

}
\references{
Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
\emph{Biometrika} 76.

Akaike, H. (1973). Information theory and an extension of the maximum likelihood
principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
USA, pp. 267–281. Springer Verlag.

Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.

Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
clustering with the integrated classification likelihood. \emph{IEEE Transactions
on Pattern Analysis and Machine Intelligence} 22.

Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
Conference on the Frontiers of Statistical Modeling: An Informational Approach:
Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.

Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.

Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
6.

Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
\href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}

Silva, A. et al. (2018). Finite Mixtures of Matrix Variate Poisson-Log Normal Distributions
for Three-Way Count Data. \href{https://arxiv.org/abs/1807.08380}{arXiv preprint arXiv:1807.08380}.
}
\author{
{Anjali Silva, \email{anjali@alumni.uoguelph.ca}, Sanjeena Dang,
         \email{sanjeenadang@cunet.carleton.ca}. }
}
