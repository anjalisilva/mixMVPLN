#' Generating Data Using Mixtures of MVPLN
#'
#' This function simulates data from a mixture of MVPLN model. Each
#' dataset will have 'n' random matrices or units, each matrix with
#' dimension r x p, where 'r' is the number of occasions and 'p' is
#' the number of responses/variables.
#'
#' @param nOccasions A positive integer indicating the number of occassions.
#'    A matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n. Here, Y_j matrix is said to contain k ∈ {1,...,p}
#'    responses/variables over i ∈ {1,...,r} occasions.
#' @param nResponses A positive integer indicating the number of responses/variables.
#'    A matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n. Here, Y_j matrix is said to contain k ∈ {1,...,p}
#'    responses/variables over i ∈ {1,...,r} occasions.
#' @param nUnits A positive integer indicating the number of units. A
#'    matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n.
#' @param mixingProportions A numeric vector that length equal to the number of total
#'    components, indicating the proportion of each component. Vector content should
#'    sum to 1.
#' @param matrixMean A matrix of size r x p for each component/cluster, giving the
#'    matrix of means (M). All matrices should be combined via rbind. See example.
#' @param phi A matrix of size r x r, which is the covariance matrix containing the
#'    variances and covariances between 'r' occasions, for each component/cluster.
#'    All matrices should be combined via rbind. See example.
#' @param omega A matrix of size p x p, which is the covariance matrix containing the
#'    variance and covariances of 'p' responses/variables, for each component/cluster.
#'    All matrices should be combined via rbind. See example.
#'
#' @return Returns an S3 object of class mvplnDataGenerator with results.
#' \itemize{
#'   \item dataset - Simulated dataset with 'n' matrices, each matrix with
#'         dimension r x p, where 'r' is the number of occasions and 'p' is
#'         the number of responses/variables.
#'   \item truemembership - A numeric vector indicating the membership of
#'         each observation.
#'   \item units - A positive integer indicating the number of units used for
#'         simulating the data.
#'   \item occassions - A positive integer indicating the number of occassions
#'         used for simulating the data.
#'   \item variables - A positive integer indicating the number of
#'         responses/variables used for simulating the data.
#'   \item mixingProportions - A numeric vector indicating the mixing
#'      proportion of each component.
#'   \item means - Matrix of mean used for simulating the data.
#'   \item phi - Covariance matrix containing the variances and covariances
#'         between 'r' occasions used for simulating the data.
#'   \item psi - Covariance matrix containing the variance and covariances
#'         of 'p' responses/variables used for simulating the data.
#'}
#'
#' @examples
#' # Generating simulated data
#' set.seed(1234)
#' true_G <- 2 # number of total G
#' true_r <- 2 # number of total occasions
#' true_p <- 3 # number of total responses
#' true_n <- 100 # number of total units
#'
#' # M is a r x p matrix
#' true_M1 <- matrix(rep(6, (true_r*true_p)),
#'                   ncol = true_p, nrow = true_r, byrow = TRUE)
#' true_M2 <- matrix(rep(1, (true_r*true_p)),
#'                   ncol = true_p, nrow = true_r, byrow = TRUE)
#' true_M_all <- rbind(true_M1, true_M2)
#'
#' # Phi is a r x r matrix
#' # Covariance matrix containing variances and covariances between r occasions
#'
#' true_Phi1 <- matrix(c(1.435610, -1.105615, -1.105615,  1.426492), nrow = 2)
#' true_Phi1[1, 1] <- 1 # for identifiability issues
#'
#' true_Phi2 <- matrix(c(0.7000000, 0.1727312, 0.1727312, 0.7000000), nrow = 2)
#' true_Phi2[1, 1] <- 1 # for identifiability issues
#' true_Phi_all <- rbind(true_Phi1, true_Phi2)

#'
#' # Omega is a p x p matrix
#' # Covariance matrix containing variance and covariances of p responses/variables
#'
#'   true_Omega1 <- matrix(c(1.4855139,  0.9049113, -0.9063446,
#'                           0.9049113,  1.3814824, -1.2298301,
#'                           -0.9063446, -1.2298301,  1.1979135), nrow = 3)
#'
#'  true_Omega2 <- matrix(c(0.7000000, 0.5379098, 0.4837924,
#'                          0.5379098, 0.7000000, 0.4089374,
#'                          0.4837924, 0.4089374, 0.7000000), nrow = 3)
#'  true_Omega_all <- rbind(true_Omega1, true_Omega2)
#'
#'
#' sampleData <- mixMVPLN::mvplnDataGenerator(nOccasions = true_r,
#'                                            nResponses = true_p,
#'                                            nUnits = true_n,
#'                                            mixingProportions = c(0.79, 0.21),
#'                                            matrixMean = true_M_all,
#'                                            phi = true_Phi_all,
#'                                            omega = true_Omega_all)
#'
#' @author Anjali Silva, \email{anjali.silva@uhnresearch.ca}
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
  if((is.numeric(nOccasions) != TRUE) || (nOccasions%%1 != 0) || (length(nOccasions) >1)) {
    stop("nOccasions argument should be an integer of class numeric.")
  }

  if((is.numeric(nResponses) != TRUE) || (nResponses%%1 != 0) || (length(nResponses) >1)) {
    stop("nResponses argument should be an integer of class numeric.")
  }

  if((is.numeric(nUnits) != TRUE) || (nUnits%%1 != 0) || (length(nUnits) >1)) {
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
  theta <- n_g <- vector("list", length = length(mixingProportions))

  # Generating theta
  for (g in 1:length(mixingProportions)) {
    n_g[[g]] <- which(z[, g] == 1)
    for (u in 1:length(n_g[[g]]) ) {
      # generating V
      theta[[n_g[[g]][u]]] <- matrix(mvtnorm::rmvnorm(n = 1,
                                                      sigma = kronecker(phi[((g - 1) * nOccasions + 1):(g * nOccasions), ],
                                                                        omega[((g - 1) * nResponses + 1):(g * nResponses), ]) ),
                                     nrow = nOccasions,
                                     ncol = nResponses)

      # Adding M1
      theta[[n_g[[g]][u]]] <-  theta[[n_g[[g]][u]]] + matrixMean[((g - 1) * nOccasions + 1):(g * nOccasions), ]
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
  UnlistDataset <- matrix(NA, ncol = nOccasions * nResponses, nrow = nUnits)
  for (u in 1:nUnits) {
    for (i in 1:nOccasions) {
      UnlistDataset[u, ((i - 1) * nResponses + 1):(i *  nResponses)] = Y[[u]][i, ]
    }
  }


  # Generate norm factors
  normFactors <- matrix(log(as.vector(edgeR::calcNormFactors(as.matrix(UnlistDataset), method = "TMM"))),
                        nrow = nOccasions,
                        ncol = nResponses,
                        byrow = TRUE) # added one because if the entire column is zero, this gives an error

  # Add norm factors and generate final counts
  for (j in 1:nUnits) {
    for (i in 1:nOccasions) {
      for (k in 1:nResponses) {
        Y2[[j]][i, k] <- stats::rpois(n = 1, lambda = exp(theta[[j]][i, k] + normFactors[i, k]))
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
                  psi = omega)
  class(results) <- "mvplnDataGenerator"
  return(results)
}

