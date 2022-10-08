#' Generating Data Using Mixtures of MVPLN
#'
#' This function simulates data from a mixture of MVPLN model. Each
#' dataset will have 'n' random matrices or units, each matrix with
#' dimension r x p, where 'r' is the number of occasions and 'p' is
#' the number of responses/variables.
#'
#' @param nOccasions A positive integer indicating the number of
#'    occassions. A matrix Y_j has size r x p, and the dataset will
#'    have 'j' such matrices with j = 1,...,n. Here, Y_j matrix is
#'    said to contain k ∈ {1,...,p} responses/variables over
#'    i ∈ {1,...,r} occasions.
#' @param nResponses A positive integer indicating the number of
#'    responses/variables. A matrix Y_j has size r x p, and the
#'    dataset will have 'j' such matrices with j = 1,...,n. Here,
#'    Y_j matrix is said to contain k ∈ {1,...,p} responses/variables
#'    over i ∈ {1,...,r} occasions.
#' @param nUnits A positive integer indicating the number of units. A
#'    matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n.
#' @param mixingProportions A numeric vector that length equal to the
#'    number of total components, indicating the proportion of each
#'    component. Vector content should sum to 1.
#' @param matrixMean A matrix of size r x p for each component/cluster,
#'    giving the matrix of means (M). All matrices should be combined
#'    via rbind. See example.
#' @param phi A matrix of size r x r, which is the covariance matrix
#'    containing the variances and covariances between 'r' occasions,
#'    for each component/cluster. All matrices should be combined via
#'    rbind. See example.
#' @param omega A matrix of size p x p, which is the covariance matrix
#'     containing the variance and covariances of 'p' responses/variables,
#'     for each component/cluster. All matrices should be combined via
#'     rbind. See example.
#' @return Returns an S3 object of class mvplnDataGenerator with results.
#' \itemize{
#'   \item dataset - Simulated dataset with 'n' matrices, each matrix
#'         with dimension r x p, where 'r' is the number of occasions
#'         and 'p' is the number of responses/variables.
#'   \item truemembership - A numeric vector indicating the membership of
#'         each observation.
#'   \item units - A positive integer indicating the number of units used
#'         for simulating the data.
#'   \item occassions - A positive integer indicating the number of
#'         occassions used for simulating the data.
#'   \item variables - A positive integer indicating the number of
#'         responses/variables used for simulating the data.
#'   \item mixingProportions - A numeric vector indicating the mixing
#'      proportion of each component.
#'   \item means - Matrix of mean used for simulating the data.
#'   \item phi - Covariance matrix containing the variances and
#'         covariances between 'r' occasions used for simulating the
#'         data.
#'   \item psi - Covariance matrix containing the variance and
#'         covariances of 'p' responses/variables used for simulating
#'         the data.
#'}
#'
#' @examples
#' # Example 1
#' # Generating simulated matrix variate count data
#' set.seed(1234)
#' trueG <- 2 # number of total G
#' truer <- 2 # number of total occasions
#' truep <- 3 # number of total responses
#' trueN <- 100 # number of total units
#' truePiG <- c(0.79, 0.21) # mixing proportions
#'
#' # Mu is a r x p matrix
#' trueM1 <- matrix(rep(6, (truer * truep)),
#'                  ncol = truep,
#'                  nrow = truer, byrow = TRUE)
#'
#' trueM2 <- matrix(rep(1, (truer * truep)),
#'                  ncol = truep,
#'                  nrow = truer,
#'                  byrow = TRUE)
#'
#' trueMall <- rbind(trueM1, trueM2)
#'
#' # Phi is a r x r matrix
#' # Loading needed packages for generating data
#' # if (!require(clusterGeneration)) install.packages("clusterGeneration")
#' # library("clusterGeneration")
#'
#' # Covariance matrix containing variances and covariances between r occasions
#' # truePhi1 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#' #                                                   dim = truer,
#' #                                                   rangeVar = c(1, 1.7))$Sigma
#' truePhi1 <- matrix(c(1.075551, -0.488301,
#'                    -0.488301, 1.362777), nrow = 2)
#' truePhi1[1, 1] <- 1 # For identifiability issues
#'
#' # truePhi2 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#' #                                                   dim = truer,
#' #                                                   rangeVar = c(0.7, 0.7))$Sigma
#' truePhi2 <- matrix(c(0.7000000, 0.6585887,
#'                      0.6585887, 0.7000000), nrow = 2)
#' truePhi2[1, 1] <- 1 # For identifiability issues
#' truePhiall <- rbind(truePhi1, truePhi2)
#'
#' # Omega is a p x p matrix
#' # Covariance matrix containing variances and covariances between p responses
#' # trueOmega1 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#' #                                    rangeVar = c(1, 1.7))$Sigma
#' trueOmega1 <- matrix(c(1.0526554, 1.0841910, -0.7976842,
#'                        1.0841910,  1.1518811, -0.8068102,
#'                        -0.7976842, -0.8068102,  1.4090578),
#'                        nrow = 3)
#' # trueOmega2 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#' #                                    rangeVar = c(0.7, 0.7))$Sigma
#' trueOmega2 <- matrix(c(0.7000000, 0.5513744, 0.4441598,
#'                        0.5513744, 0.7000000, 0.4726577,
#'                        0.4441598, 0.4726577, 0.7000000),
#'                        nrow = 3)
#' trueOmegaAll <- rbind(trueOmega1, trueOmega2)
#'
#' # Generated simulated data
#' sampleData <- mixMVPLN::mvplnDataGenerator(nOccasions = truer,
#'                                            nResponses = truep,
#'                                            nUnits = trueN,
#'                                            mixingProportions = truePiG,
#'                                            matrixMean = trueMall,
#'                                            phi = truePhiall,
#'                                            omega = trueOmegaAll)
#'
#' # Example 2
#' trueG <- 1 # number of total G
#' truer <- 2 # number of total occasions
#' truep <- 3 # number of total responses
#' trueN <- 1000 # number of total units
#' truePiG <- 1L # mixing proportion for G = 1
#'
#' # Mu is a r x p matrix
#' trueM1 <- matrix(c(6, 5.5, 6, 6, 5.5, 6),
#'                  ncol = truep,
#'                  nrow = truer,
#'                  byrow = TRUE)
#' trueMall <- rbind(trueM1)

#' # Phi is a r x r matrix
#' set.seed(1)
#' # truePhi1 <- clusterGeneration::genPositiveDefMat(
#' #                               "unifcorrmat",
#' #                                dim = truer,
#' #                               rangeVar = c(0.7, 1.7))$Sigma
#' truePhi1 <- matrix(c(1.3092747, 0.3219674,
#'                      0.3219674, 1.3233794), nrow = 2)
#' truePhi1[1, 1] <- 1 # for identifiability issues
#' truePhiall <- rbind(truePhi1)
#'
#' # Omega is a p x p matrix
#' set.seed(1)
#' # trueOmega1 <- genPositiveDefMat(
#' #                    "unifcorrmat",
#' #                     dim = truep,
#' #                     rangeVar = c(1, 1.7))$Sigma
#' trueOmega1 <- matrix(c(1.1625581, 0.9157741, 0.8203499,
#'                        0.9157741, 1.2216287, 0.7108193,
#'                        0.8203499, 0.7108193, 1.2118854), nrow = 3)
#' trueOmegaAll <- rbind(trueOmega1)
#'
#' # Generated simulated data
#' set.seed(1)
#' sampleData2 <- mixMVPLN::mvplnDataGenerator(
#'                          nOccasions = truer,
#'                          nResponses = truep,
#'                          nUnits = trueN,
#'                          mixingProportions = truePiG,
#'                          matrixMean = trueMall,
#'                          phi = truePhiall,
#'                          omega = trueOmegaAll)
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}, Sanjeena Dang,
#'          \email{sanjeenadang@cunet.carleton.ca}. }
#'
#' @references
#' Silva, A. et al. (2018). Finite Mixtures of Matrix Variate Poisson-Log Normal Distributions
#' for Three-Way Count Data. \href{https://arxiv.org/abs/1807.08380}{arXiv preprint arXiv:1807.08380}.
#'
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76. \href{https://www.jstor.org/stable/2336624?seq=1}{Link}.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}.
#'
#' @export
#' @import stats
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
mvplnDataGenerator <- function(nOccasions,
                               nResponses,
                               nUnits,
                               mixingProportions,
                               matrixMean,
                               phi,
                               omega) {
  # Checking user input
  if(is.numeric(nOccasions) != TRUE) {
    stop("nOccasions argument should be an integer of class numeric.")
  }

  if(nOccasions%%1 != 0) {
    stop("nOccasions argument should be an integer of class numeric.")
  }

  if(length(nOccasions) > 1) {
    stop("nOccasions argument should be an integer of class numeric.")
  }

  if(is.numeric(nResponses) != TRUE) {
     stop("nResponses argument should be an integer of class numeric.")
  }

  if(nResponses%%1 != 0) {
    stop("nResponses argument should be an integer of class numeric.")
  }

  if(length(nResponses) > 1) {
    stop("nResponses argument should be an integer of class numeric.")
  }

  if(is.numeric(nUnits) != TRUE) {
    stop("nUnits argument should be an integer of class numeric.")
  }

  if(nUnits%%1 != 0) {
    stop("nUnits argument should be an integer of class numeric.")
  }

  if(length(nUnits) > 1) {
    stop("nUnits argument should be an integer of class numeric.")
  }

  if (sum(mixingProportions) != 1) {
    stop("mixingProportions argument should be a vector that sum to 1.")
  }

  if(is.matrix(matrixMean) != TRUE) {
    stop("matrixMean argument should be of class matrix.")
  }

  if(ncol(matrixMean) != nResponses) {
    stop("matrixMean argument should a matrix of size r x p for each component/cluster.")
  }

  if(is.matrix(phi) != TRUE) {
    stop("phi argument should be of class matrix.")
  }

  if(ncol(phi) != nOccasions) {
    stop("phi argument should be a matrix of size r x r for each component/cluster.")
  }

  if(is.matrix(omega) != TRUE) {
    stop("omega argument should be of class matrix.")
  }

  if(ncol(omega) != nResponses) {
    stop("omega argument should be a matrix of size p x p for each component/cluster.")
  }


  # Begin calculations - generate z
  z <- t(stats::rmultinom(n = nUnits, size = 1, prob = mixingProportions))

  # For saving data
  Y <- Y2 <- normFactors <- list()
  theta <- nG <- vector("list", length = length(mixingProportions))

  # Generating theta
  for (g in 1:length(mixingProportions)) {
    nG[[g]] <- which(z[, g] == 1)
    for (u in 1:length(nG[[g]]) ) {
      # generating V
      theta[[nG[[g]][u]]] <- matrix(mvtnorm::rmvnorm(n = 1,
                                                      sigma = kronecker(phi[((g - 1) * nOccasions + 1):(g * nOccasions), ],
                                                                        omega[((g - 1) * nResponses + 1):(g * nResponses), ]) ),
                                     nrow = nOccasions,
                                     ncol = nResponses,
                                     byrow = TRUE)

      # Adding M1
      theta[[nG[[g]][u]]] <-  theta[[nG[[g]][u]]] +
                              matrixMean[((g - 1) * nOccasions + 1):(g * nOccasions), ]
    }
  }

  # Doing exp(theta) to get counts
  for (j in 1:nUnits) {
    Y[[j]] <- Y2[[j]] <- matrix(ncol = nResponses, nrow = nOccasions)
    for (i in 1:nOccasions) {
      for (k in 1:nResponses) {
        Y[[j]][i, k] <- stats::rpois(n = 1, lambda = exp(theta[[j]][i, k]))
      }
    }
  }

  # Unlist to generate 2D dataset for norm factor generation
  UnlistDataset <- matrix(NA,
                          ncol = nOccasions * nResponses,
                          nrow = nUnits)
  for (u in 1:nUnits) {
    for (i in 1:nOccasions) {
      UnlistDataset[u, ((i - 1) * nResponses + 1):(i *  nResponses)] = Y[[u]][i, ]
    }
  }


  # Generate norm factors
  normFactors <- matrix(log(as.vector(edgeR::calcNormFactors(as.matrix(UnlistDataset),
                        method = "TMM"))),
                        nrow = nOccasions,
                        ncol = nResponses,
                        byrow = TRUE) # added one because if the entire
                                      # column is zero, this gives an error

  # Add norm factors and generate final counts
  for (j in 1:nUnits) {
    for (i in 1:nOccasions) {
      for (k in 1:nResponses) {
        Y2[[j]][i, k] <- stats::rpois(n = 1,
                                      lambda = exp(theta[[j]][i, k] +
                                      normFactors[i, k]))
      }
    }
  }

  results <- list(dataset = Y2,
                  truemembership = mclust::map(z),
                  units = nUnits,
                  occassions = nOccasions,
                  variables = nResponses,
                  mixingProportions = mixingProportions,
                  means = matrixMean,
                  phi = phi,
                  omega = omega,
                  data2D = UnlistDataset)
  class(results) <- "mvplnDataGenerator"
  return(results)
}
# [END]
