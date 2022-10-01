#' Clustering Using mixtures of MVPLN via MCMC-EM
#'
#' Performs clustering using mixtures of matrix variate Poisson-log normal
#' (MVPLN) via MCMC-EM with parallelization. Model selection can be done
#' using AIC, AIC3, BIC and ICL.
#'
#' @param dataset A list of length nUnits, containing Y_j matrices.
#'    A matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n. If a Y_j has all zeros, such Y_j will
#'    be removed prior to cluster analysis.
#' @param membership A numeric vector of size length(dataset) containing the
#'    cluster membership of each Y_j matrix. If not available, leave as
#'    "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each MCMC chain (including warmup). The value should be greater than
#'    40. The upper limit will depend on size of dataset.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be considered. Default is 0.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#' @param numNodes A positive integer indicating the number of nodes to be
#'     used from the local machine to run the clustering algorithm. Else
#'     leave as NA, so default will be detected as
#'     parallel::makeCluster(parallel::detectCores() - 1).
#'
#' @return Returns an S3 object of class mvplnParallel with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item nUnits - Number of units in the input dataset.
#'   \item nVariables - Number of variables in the input dataset.
#'   \item nOccassions - Number of occassions in the input dataset.
#'   \item normFactors - A vector of normalization factors used for
#'      input dataset.
#'   \item gmin - Minimum number of components considered in the clustering
#'      run.
#'   \item gmax - Maximum number of components considered in the clustering
#'      run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item loglikelihood - A vector with value of final log-likelihoods for
#'      each cluster size.
#'   \item nParameters - A vector with number of parameters for each
#'      cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICL_all - A list with all ICL model selection results.
#'   \item BIC_all - A list with all BIC model selection results.
#'   \item AIC_all - A list with all AIC model selection results.
#'   \item AIC3_all - A list with all AIC3 model selection results.
#'   \item totalTime - Total time used for clustering and model selection.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1
#' # Generating simulated matrix variate count data
#' set.seed(1234)
#' trueG <- 2 # number of total G
#' truer <- 2 # number of total occasions
#' truep <- 3 # number of total responses
#' trueN <- 100 # number of total units
#' truePiG <- c(0.6, 0.4) # mixing proportions for G=2
#'
#' # Mu is a r x p matrix
#' trueM1 <- matrix(rep(6.2, (truer * truep)),
#'                  ncol = truep,
#'                  nrow = truer,
#'                  byrow = TRUE)
#'
#' trueM2 <- matrix(rep(1.5, (truer * truep)),
#'                  ncol = truep,
#'                  nrow = truer,
#'                  byrow = TRUE)
#'
#' trueMall <- rbind(trueM1, trueM2)
#'
#' # Phi is a r x r matrix
#' # Loading needed packages for generating data
#' # install.packages("clusterGeneration")
#' # library("clusterGeneration")
#'
#' set.seed(1)
#' truePhi1 <- matrix(rep(0, truer * truer), truer, truer)
#' diag(truePhi1) <- diag(clusterGeneration::genPositiveDefMat(
#'                        "unifcorrmat",
#'                         dim = truer,
#'                         rangeVar = c(1, 1.7))$Sigma)
#' truePhi1[1, 1] <- 1 # for identifiability issuess
#' truePhi2 <- matrix(rep(0,truer * truer), truer, truer)
#' diag(truePhi2) <- diag(clusterGeneration::genPositiveDefMat(
#'                        "unifcorrmat",
#'                         dim = truer,
#'                         rangeVar = c(0.7, 0.7))$Sigma)
#' truePhi2[1, 1] <- 1 # for identifiability issues
#' truePhiall <- rbind(truePhi1, truePhi2)
#'
#'
#' # Omega is a p x p matrix
#' set.seed(1)
#' trueOmega1 <- matrix(rep(0, truep * truep), truep, truep)
#' diag(trueOmega1) <- diag(clusterGeneration::genPositiveDefMat(
#'                          "unifcorrmat",
#'                           dim = truep,
#'                           rangeVar = c(1, 1.7))$Sigma)
#' trueOmega2 <- matrix(rep(0, truep * truep), truep, truep)
#' diag(trueOmega2) <- diag(clusterGeneration::genPositiveDefMat(
#'                          "unifcorrmat",
#'                           dim = truep,
#'                          mrangeVar = c(0.6, 0.9))$Sigma)
#' trueOmegaAll <- rbind(trueOmega1, trueOmega2)
#'
#' # Generated simulated data
#' sampleData <- mixMVPLN::mvplnDataGenerator(
#'                          nOccasions = truer,
#'                          nResponses = truep,
#'                          nUnits = trueN,
#'                          mixingProportions = truePiG,
#'                          matrixMean = trueMall,
#'                          phi = truePhiall,
#'                          omega = trueOmegaAll)
#'
#'
#' # Clustering simulated matrix variate count data
#' clusteringResults <- mixMVPLN::mvplnMCMCclus(dataset = sampleData$dataset,
#'                                       membership = sampleData$truemembership,
#'                                       gmin = 1,
#'                                       gmax = 3,
#'                                       nChains = 3,
#'                                       nIterations = 300,
#'                                       initMethod = "kmeans",
#'                                       nInitIterations = 1,
#'                                       normalize = "Yes")
#'
#'
#' # Example 2
#' trueG <- 1 # 1 cluster
#' truer <- 2 # variety
#' truep <- 3 # growth stages
#' truen <- 100 # genes
#' truePiG <- c(1) # mixing proportions for G = 1
#'
#' # Mu is a r x p matrix
#' trueM1 <- matrix(c(6, 5.5, 6, 6, 5.5, 6),
#'                  ncol = truep,
#'                  nrow = truer,
#'                  byrow = TRUE)
#' trueMall <- rbind(trueM1)

#' # Phi is a r x r matrix
#' set.seed(1)
#' truePhi1 <- clusterGeneration::genPositiveDefMat(
#'                               "unifcorrmat",
#'                               dim = truer,
#'                               rangeVar = c(0.7, 1.7))$Sigma
#' truePhi1[1, 1] <- 1 # for identifiability issues
#' truePhiall <- rbind(truePhi1)
#'
#' # Omega is a p x p matrix
#' set.seed(1)
#' trueOmega1 <- genPositiveDefMat(
#'                    "unifcorrmat",
#'                     dim = truep,
#'                     rangeVar = c(1, 1.7))$Sigma
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
#' # Clustering simulated matrix variate count data
#' clusteringResults <- mixMVPLN::mvplnMCMCclus(dataset = sampleData2$dataset,
#'                                       membership = sampleData2$truemembership,
#'                                       gmin = 1,
#'                                       gmax = 3,
#'                                       nChains = 3,
#'                                       nIterations = 300,
#'                                       initMethod = "kmeans",
#'                                       nInitIterations = 1,
#'                                       normalize = "Yes")
#' }
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}, Sanjeena Dang,
#'          \email{sanjeenadang@cunet.carleton.ca}. }
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Silva, A. et al. (2018). Finite Mixtures of Matrix Variate Poisson-Log Normal Distributions
#' for Three-Way Count Data. \href{https://arxiv.org/abs/1807.08380}{arXiv preprint arXiv:1807.08380}.
#'
#' @export
#' @import coda
#' @import cluster
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstan stan_model
#' @import parallel
#' @importFrom utils tail
#'
mvplnMCMCclus <- function(dataset,
                            membership = "none",
                            gmin,
                            gmax,
                            nChains = 3,
                            nIterations = NA,
                            initMethod = "kmeans",
                            nInitIterations = 2,
                            normalize = "Yes",
                            numNodes = NA) {

  ptm <- base::proc.time()

  # Performing checks
  if (typeof(unlist(dataset)) != "double" & typeof(unlist(dataset)) != "integer") {
    stop("dataset should be a list of count matrices.");}

  if (any((unlist(dataset) %% 1 == 0) == FALSE)) {
    stop("dataset should be a list of count matrices.")
  }

  if (is.list(dataset) != TRUE) {
    stop("dataset needs to be a list of matrices.")
  }

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > length(dataset)) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (is.numeric(nChains) != TRUE) {
    stop("nChains should be a positive integer of class numeric specifying
      the number of Markov chains.")
  }

  if(nChains < 3) {
    cat("nChains is less than 3. Note: recommended minimum number of
      chains is 3.")
  }

  if(is.numeric(nIterations) != TRUE) {
    stop("nIterations should be a positive integer of class numeric,
      specifying the number of iterations for each Markov chain
      (including warmup).")
  }

  if(is.numeric(nIterations) == TRUE && nIterations < 40) {
    stop("nIterations argument should be greater than 40.")
  }

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }


  # Checking if missing membership values
  # First check for the case in which G = 1, otherwise check
  #   if missing cluster memberships
  if(all(membership != "none") && length(unique(membership)) != 1) {
    if(all(membership != "none") &&
       all((diff(sort(unique(membership))) == 1) != TRUE) ) {
      stop("Cluster memberships in the membership vector
        are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
    }
  }


  if(all(membership != "none") && length(membership) != length(dataset)) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  if (is.character(initMethod) == TRUE) {
    if(initMethod != "kmeans" & initMethod != "random" & initMethod != "medoids" & initMethod != "medoids" & initMethod != "clara" & initMethod != "fanny") {
      stop("initMethod should of class character, specifying
      either: kmeans, random, medoids, clara, or fanny.")
    }
  } else if (is.character(initMethod) != TRUE) {
    stop("initMethod should of class character, specifying
      either: kmeans, random, medoids, clara, or fanny.")
  }

  if (is.numeric(nInitIterations) != TRUE) {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (is.character(normalize) != TRUE) {
    stop("normalize should be a string of class character specifying
      if normalization should be performed.")
  }

  if (normalize != "Yes" && normalize != "No") {
    stop("normalize should be a string indicating Yes or No, specifying
       if normalization should be performed.")
  }

  if((is.logical(numNodes) != TRUE) && (is.numeric(numNodes) != TRUE)) {
    stop("numNodes should be a positive integer indicating the number
         of nodes to be used from the local machine to run the
         clustering algorithm. Else leave as NA.")
  }


  # Check numNodes and calculate the number of cores and initiate cluster
  if(is.na(numNodes) == TRUE) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
  } else if (is.numeric(numNodes) == TRUE) {
    cl <- parallel::makeCluster(numNodes)
  } else {
    stop("numNodes should be a positive integer indicating the number of
        nodes to be used from the local machine. Else leave as NA, so the
        numNodes will be determined by
        parallel::makeCluster(parallel::detectCores() - 1).")
  }

  n <- length(dataset)
  p <- ncol(dataset[[1]])
  r <- nrow(dataset[[1]])
  d <- p * r


  # Changing dataset into a n x rp dataset
  TwoDdataset <- matrix(NA, ncol = d, nrow = n)
  sample_matrix <- matrix(c(1:d), nrow = r, byrow = TRUE)
  for (u in 1:n) {
    for (e in 1:p) {
      for (s in 1:r) {
        TwoDdataset[u, sample_matrix[s, e]] <- dataset[[u]][s, e]
      }
    }
  }

  # Check if entire row is zero
  if(any(rowSums(TwoDdataset) == 0)) {
    TwoDdataset <- TwoDdataset[- which(rowSums(TwoDdataset) == 0), ]
    n <- nrow(TwoDdataset)
    membership <- membership[- which(rowSums(TwoDdataset) == 0)]
  }

  if(all(is.na(membership) == TRUE)) {
    membership <- "Not provided" }


  # Calculating normalization factors
  if(normalize == "Yes") {
    normFactors <- log(as.vector(edgeR::calcNormFactors(as.matrix(TwoDdataset),
                                                        method = "TMM")))
  } else if(normalize == "No") {
    normFactors <- rep(0, d)
  } else{
    stop("Argument normalize should be 'Yes' or 'No' ")
  }


  # Construct a Stan model
  # Stan model was developed with help of Sanjeena Dang
  # <sdang@math.binghamton.edu>
  stancode <- 'data{int<lower=1> d; // Dimension of theta
  int<lower=0> N; //Sample size
  int y[N,d]; //Array of Y
  vector[d] mu;
  matrix[d,d] Sigma;
  vector[d] normfactors;
  }
  parameters{
  matrix[N,d] theta;
  }
  model{ //Change values for the priors as appropriate
  for (n in 1:N){
  theta[n,]~multi_normal(mu,Sigma);
  }
  for (n in 1:N){
  for (k in 1:d){
  real z;
  z=exp(normfactors[k]+theta[n,k]);
  y[n,k]~poisson(z);
  }
  }
  }'

  mod <- rstan::stan_model(model_code = stancode, verbose = FALSE)


  # Constructing parallel code
  mvplnParallelCode <- function(g) {
    ## ** Never use set.seed(), use clusterSetRNGStream() instead,
    # to set the cluster seed if you want reproducible results
    # clusterSetRNGStream(cl = cl, iseed = g)
    mvplnParallelRun <- callingClustering(n = n,
                                          r = r,
                                          p = p,
                                          d = d,
                                          dataset = TwoDdataset,
                                          gmin = g,
                                          gmax = g,
                                          nChains = nChains,
                                          nIterations = nIterations,
                                          initMethod = initMethod,
                                          nInitIterations = nInitIterations,
                                          normFactors = normFactors,
                                          model = mod)
    return(mvplnParallelRun)
  }



  # Doing clusterExport
  parallel::clusterExport(cl = cl,
                          varlist = c("callingClustering",
                                      "calcLikelihood",
                                      "calcZvalue",
                                      "dataset",
                                      "initMethod",
                                      "initializationRun",
                                      "mvplnParallelCode",
                                      "mvplnCluster",
                                      "mod",
                                      "nInitIterations",
                                      "normFactors",
                                      "nChains",
                                      "nIterations",
                                      "parameterEstimation",
                                      "p",
                                      "r",
                                      "stanRun",
                                      "TwoDdataset"),
                          envir = environment())

  # Doing clusterEvalQ
  parallel::clusterEvalQ(cl, library(parallel))
  parallel::clusterEvalQ(cl, library(rstan))
  parallel::clusterEvalQ(cl, library(Rcpp))
  parallel::clusterEvalQ(cl, library(mclust))
  parallel::clusterEvalQ(cl, library(mvtnorm))
  parallel::clusterEvalQ(cl, library(edgeR))
  parallel::clusterEvalQ(cl, library(clusterGeneration))
  parallel::clusterEvalQ(cl, library(coda))




  parallelRun <- list()
  cat("\nRunning parallel code...")
  parallelRun <- parallel::clusterMap(cl = cl,
                                     fun = mvplnParallelCode,
                                     g = gmin:gmax)
  cat("\nDone parallel code.")
  #parallel::stopCluster(cl)


  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- k <- ll <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }



    # save the final log-likelihood
    ll[g] <- unlist(tail(parallelRun[[g]]$allresults$loglikelihood, n = 1))

    k[g] <- calcParameters(g = clustersize,
                           r = r,
                           p = p)

    # starting model selection
    if (g == max(1:(gmax - gmin + 1))) {
      bic <- BICFunction(ll = ll,
                         k = k,
                         n = n,
                         run = parallelRun,
                         gmin = gmin,
                         gmax = gmax)

      icl <- ICLFunction(bIc = bic,
                         gmin = gmin,
                         gmax = gmax,
                         run = parallelRun)

      aic <- AICFunction(ll = ll,
                         k = k,
                         run = parallelRun,
                         gmin = gmin,
                         gmax = gmax )

      aic3 <- AIC3Function(ll = ll,
                           k = k,
                           run = parallelRun,
                           gmin = gmin,
                           gmax = gmax)
    }
  }



  final <- proc.time() - ptm
  RESULTS <- list(dataset = dataset,
                  nUnits = n,
                  nVariables = p,
                  nOccassions = r,
                  normFactors = normFactors,
                  gmin = gmin,
                  gmax = gmax,
                  initalizationMethod = initMethod,
                  allResults = parallelRun,
                  loglikelihood = ll,
                  nParameters = k,
                  trueLabels = membership,
                  ICL.all = icl,
                  BIC.all = bic,
                  AIC.all = aic,
                  AIC3.all = aic3,
                  totalTime = final)

  class(RESULTS) <- "mvplnParallel"
  return(RESULTS)
}

# Log likelihood calculation
calcLikelihood <- function(dataset,
                           z,
                           G,
                           PI,
                           normFactors,
                           M,
                           Sigma,
                           thetaStan) {
  n <- nrow(dataset)
  like <- matrix(NA, nrow = n, ncol = G)
  d <- ncol(dataset)
  like <- sapply(c(1:G), function(g) sapply(c(1:n), function(i) z[i, g] * (log(PI[g]) +
                                                                             t(dataset[i, ]) %*% (thetaStan[[g]][i, ] + normFactors) -
                                                                             sum(exp(thetaStan[[g]][i, ] + normFactors)) - sum(lfactorial(dataset[i, ])) -
                                                                             d / 2 * log(2 * pi) - 1 / 2 * log(det(Sigma[((g - 1) * d + 1):(g * d), ])) -
                                                                             0.5 * t(thetaStan[[g]][i, ] -
                                                                                       M[g, ]) %*% solve(Sigma[((g - 1) * d + 1):(g * d), ])
                                                                           %*% (thetaStan[[g]][i, ] - M[g, ])) ))
  loglike <- sum(rowSums(like))
  return(loglike)
  # Developed by Anjali Silva
}

# Parameter estimation
parameterEstimation <- function(r,
                                p,
                                G,
                                z,
                                fit,
                                nIterations,
                                nChains) {
  d <- r * p
  n <- nrow(z)

  M <- matrix(NA, ncol = p, nrow = r * G) #initialize M = rxp matrix
  Sigma <- do.call("rbind", rep(list(diag(r * p)), G)) # sigma
  Phi <- do.call("rbind", rep(list(diag(r)), G))
  Omega <- do.call("rbind", rep(list(diag(p)), G))

  # generating stan values
  EThetaOmega <- E_theta_phi <- E_theta_phi2 <- list()

  # update parameters
  theta_mat <- lapply(as.list(c(1:G)),
                      function(g) lapply(as.list(c(1:n)),
                                         function(o) as.matrix(fit[[g]])[, c(o,sapply(c(1:n),
                                                                                      function(i) c(1:(d - 1)) * n + i)[, o]) ]) )

  thetaStan <- lapply(as.list(c(1:G)),
                       function(g) t(sapply(c(1:n),
                                            function(o) colMeans(as.matrix(fit[[g]])[, c(o, sapply(c(1:n),
                                                                                                   function(i) c(1:(d-1))*n+i)[,o]) ]))) )


  # updating mu
  M <- lapply(as.list(c(1:G)), function(g) matrix(data = colSums(z[, g] * thetaStan[[g]]) / sum(z[, g]),
                                                  ncol=p,
                                                  nrow=r,
                                                  byrow = T) )
  M <- do.call(rbind,M)

  # updating phi
  E_theta_phi <- lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)),
                                                            function(i) lapply(as.list(c(1:nrow(theta_mat[[g]][[i]]))),
                                                                               function(e) ((matrix(theta_mat[[g]][[i]][e,], r, p, byrow = T)) - M[((g - 1) * r + 1):(g * r), ]) %*% solve(Omega[((g - 1) * p + 1):(g*p), ]) %*% t((matrix(theta_mat[[g]][[i]][e, ], r, p, byrow = T)) - M[((g - 1) * r + 1):(g * r), ]) ) ) )
  E_theta_phi2 <- lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)),
                                                             function(i) z[i, g] * Reduce("+", E_theta_phi[[g]][[i]]) / ((0.5 * nIterations) * nChains) ) )
  Phi <- lapply(as.list(c(1:G)), function(g) (Reduce("+", E_theta_phi2[[g]]) / (p * sum(z[, g]))))
  Phi <- do.call(rbind, Phi)

  # updating omega
  EThetaOmega <- lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)),
                                                              function(i) lapply(as.list(1:nrow(theta_mat[[g]][[i]])),
                                                                                 function(e) t((matrix(theta_mat[[g]][[i]][e,], r, p, byrow = T)) - M[((g - 1) * r + 1):(g * r), ]) %*% solve(Phi[((g - 1) * r + 1):(g * r), ]) %*% ((matrix(theta_mat[[g]][[i]][e, ], r, p, byrow = T)) - M[((g - 1) * r + 1):(g * r), ]) ) ) )
  EThetaOmega2 <- lapply(as.list(c(1:G)) ,function(g) lapply(as.list(c(1:n)),
                                                               function(i)  z[i, g] * Reduce("+", EThetaOmega[[g]][[i]]) / ((0.5 * nIterations) * nChains) ) )
  Omega <- lapply(as.list(c(1:G)), function(g) Reduce("+", EThetaOmega2[[g]]) / (r * sum(z[, g]))  )
  Omega <- do.call(rbind, Omega)

  Sigma <- lapply(as.list(c(1:G)), function(g) kronecker(Phi[((g - 1) * r + 1):(g * r), ], Omega[((g - 1) * p + 1):(g * p), ]) )
  Sigma <- do.call(rbind,Sigma)

  results <- list(Mu = M,
                  Sigma = Sigma,
                  Omega = Omega,
                  Phi = Phi,
                  theta = thetaStan)
  class(results) <- "RStan"
  return(results)
  # Developed by Anjali Silva
}


# Parameter calculation
calcParameters <- function(g,
                           r,
                           p) {

  muPara <- r * p * g
  sigmaPara <- g / 2 *( (r * (r + 1)) + (p * (p + 1)))
  piPara <- g - 1 # If g-1 parameters, you can do 1-these to get the last one

  # total parameters are
  paraTotal <- muPara + sigmaPara + piPara

  return(paraTotal)
  # Developed by Anjali Silva
}


# z value calculation
calcZvalue <- function(PI,
                       dataset,
                       M,
                       G,
                       Sigma,
                       thetaStan,
                       normFactors) {

  d <- ncol(dataset)
  n <- nrow(dataset)
  forz <- sapply(c(1:G), function(g) sapply(c(1:n), function(i) PI[g] * exp(t(dataset[i,]) %*% (thetaStan[[g]][i, ] + normFactors)
                                                                            - sum(exp(thetaStan[[g]][i, ] + normFactors))
                                                                            - sum(lfactorial(dataset[i, ])) -
                                                                              d / 2 * log(2 * pi) - 1 / 2 * log(det(Sigma[((g - 1) * d + 1):(g * d), ])) -
                                                                              0.5 * t(thetaStan[[g]][i, ] - M[g, ]) %*% solve(Sigma[((g - 1) * d + 1):(g * d), ]) %*%
                                                                              (thetaStan[[g]][i, ] - M[g, ])) ))

  if (G == 1) {
    errorpossible <- Reduce(intersect, list(which(forz == 0), which(rowSums(forz) == 0)))
    zvalue <- forz / rowSums(forz)
    zvalue[errorpossible, ] <- 1
  }else {zvalue <- forz / rowSums(forz)}

  return(zvalue)
  # Developed by Anjali Silva
}

# Calling clustering
callingClustering <- function(n, r, p, d,
                              dataset,
                              gmin,
                              gmax,
                              nChains,
                              nIterations = NA,
                              initMethod = NA,
                              nInitIterations = NA,
                              normFactors,
                              model) {
  ptm_inner <- proc.time()

  for (gmodel in 1:(gmax - gmin + 1)) {
    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }

    if(nInitIterations != 0) {
      initializeruns <- initializationRun(r = r,
                                          p = p,
                                          gmodel = clustersize,
                                          mod = model,
                                          dataset = dataset,
                                          initMethod = initMethod,
                                          nInitIterations = nInitIterations,
                                          nChains = nChains,
                                          nIterations = nIterations,
                                          normFactors = normFactors)
      allruns <- mvplnCluster(r = r,
                              p = p,
                              z = NA,
                              dataset = dataset,
                              G = clustersize,
                              mod = model,
                              normFactors = normFactors,
                              nChains = nChains,
                              nIterations = NA,
                              initialization = initializeruns)
    } else if(nInitIterations == 0) {
      allruns <- mvplnCluster( r = r,
                               p = p,
                               z = mclust::unmap(stats::kmeans(x = log(dataset + 1 / 3),
                                                               centers = clustersize)$cluster),
                               dataset = dataset,
                               G = clustersize,
                               mod = model,
                               normFactors = normFactors,
                               nChains = nChains,
                               nIterations = nIterations,
                               initialization = NA)
    }
  }

  finalInner <- proc.time() - ptm_inner
  RESULTS <- list(gmin = gmin,
                  gmax = gmax,
                  initalization_method = initMethod,
                  allresults = allruns,
                  totaltime = finalInner)


  class(RESULTS) <- "mvplnCallingClustering"
  return(RESULTS)
  # Developed by Anjali Silva
}

# Initialization
initializationRun <- function(r,
                              p,
                              dataset,
                              gmodel,
                              mod,
                              normFactors,
                              nChains,
                              nIterations,
                              initMethod,
                              nInitIterations) {
  z <- init_runs <- list()
  logLInit <- vector()
  n <- nrow(dataset)
  d <- r*p

  for(iterations in seq_along(1:nInitIterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)

    if (initMethod == "kmeans" | is.na(initMethod)) {
      z[[iterations]] <- mclust::unmap(stats::kmeans(x = log(dataset + 1 / 3),
                                                     centers = gmodel)$cluster)
    } else if (initMethod == "random") {
      if(gmodel == 1) { # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times = n),
                                     ncol = gmodel,
                                     nrow = n)
      } else { # generating z if g>1
        zConv = 0
        while(! zConv) { # ensure that dimension of z is same as G (i.e.
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(stats::rmultinom(n = n,
                                                size = 1,
                                                prob = rep(1 / gmodel, gmodel)))
          if(length(which(colSums(z[[iterations]]) > 0)) == gmodel) {
            zConv <- 1
          }
        }
      }
    } else if (initMethod == "medoids") {
      z[[iterations]] <- mclust::unmap(cluster::pam(x = log(dataset + 1 / 3),
                                                    k = gmodel)$cluster)
    } else if (initMethod == "clara") {
      z[[iterations]] <- mclust::unmap(cluster::clara(x = log(dataset + 1 / 3),
                                                      k = gmodel)$cluster)
    } else if (initMethod == "fanny") {
      z[[iterations]] <- mclust::unmap(cluster::fanny(x = log(dataset + 1 / 3),
                                                      k = gmodel)$cluster)
    }

    init_runs[[iterations]] <- mvplnCluster(r = r,
                                            p = p,
                                            z = z[[iterations]],
                                            dataset = dataset,
                                            G = gmodel,
                                            mod = mod,
                                            normFactors = normFactors,
                                            nChains = nChains,
                                            nIterations = nIterations,
                                            initialization = "init")

    logLInit[iterations] <-
      unlist(utils::tail((init_runs[[iterations]]$loglikelihood), n = 1))
  }

  initialization <- init_runs[[which(logLInit == max(logLInit, na.rm = TRUE))[1]]]

  return(initialization)
  # Developed by Anjali Silva
}



# clustering function
mvplnCluster <- function(r, p, z,
                         dataset,
                         G,
                         mod,
                         normFactors,
                         nChains,
                         nIterations,
                         initialization) {
  n <- nrow(dataset)

  PhiAllOuter <- OmegaAllOuter <- MAllOuter <- SigmaAllOuter <- list()
  medianMuOuter <- medianSigmaOuter <- medianPhiOuter <- medianOmegaOuter <- list()
  convOuter <- 0
  itOuter <- 2
  obs <- PI <- logL <- normMuOuter <- normSigmaOuter <- vector()

  if (all(is.na(initialization)) == TRUE  || all(initialization == "init")) {
    # initialize Phi; rxr times G
    PhiAllOuter[[1]] <- Phi <- do.call("rbind", rep(list(diag(r)), G))

    # initialize Omega; pxp times G
    OmegaAllOuter[[1]] <- Omega <- do.call("rbind", rep(list(diag(p)), G))

    M <- matrix(NA, ncol = p, nrow = r*G) # initialize M = rxp matrix
    Sigma <- do.call("rbind", rep(list(diag(r * p)), G)) # sigma (rp by rp)
    for (g in seq_along(1:G)) {
      M[((g - 1) * r + 1):(g * r), ] <- matrix(log(mean(dataset[c(which(z[, g] == 1)), ])),
                                               ncol = p,
                                               nrow = r)
      Sigma[((g - 1) * (r * p) + 1):(g * (r * p)), ] <- (cov(log(dataset[c(which(z[, g] == 1)), ] + (1 / 3)))) #diag(r*p)
    }
    MAllOuter[[1]] <- M
    SigmaAllOuter[[1]] <- Sigma
  } else{ # running after initialization has been done
    MAllOuter[[1]] <- M <- initialization$finalmu
    PhiAllOuter[[1]] <- Phi <- initialization$finalphi
    OmegaAllOuter[[1]] <- Omega <- initialization$finalomega
    SigmaAllOuter[[1]] <- Sigma <- initialization$finalsigma
    z <- initialization$probaPost
    nIterations <- initialization$FinalRstanIterations
  }


  # start the loop
  while(! convOuter) {
    obs <- apply(z, 2, sum) # number of observations in each group
    PI <- sapply(obs, function(x) x / n)  # obtain probability of each group

    stanresults <- stanRun(r = r,
                           p = p,
                           dataset = dataset,
                           G = G,
                           nChains = nChains,
                           nIterations = nIterations,
                           mod = mod,
                           Mu = MAllOuter[[itOuter-1]],
                           Sigma = SigmaAllOuter[[itOuter-1]],
                           normFactors = normFactors)

    # update parameters
    paras <- parameterEstimation(r = r,
                            p = p,
                            G = G,
                            z = z,
                            fit = stanresults$fitrstan,
                            nIterations = stanresults$nIterations,
                            nChains = nChains)


    MAllOuter[[itOuter]] <- paras$Mu
    SigmaAllOuter[[itOuter]] <- paras$Sigma
    PhiAllOuter[[itOuter]] <- paras$Phi
    OmegaAllOuter[[itOuter]] <- paras$Omega

    thetaStan <- paras$theta

    vectorizedM <- t(sapply(c(1:G), function(g)
      ( MAllOuter[[itOuter]][((g - 1) * r + 1):(g * r), ]) ))

    logL[itOuter] <- calcLikelihood(dataset = dataset,
                                     z = z,
                                     G = G,
                                     PI = PI,
                                     normFactors = normFactors,
                                     M = vectorizedM,
                                     Sigma = SigmaAllOuter[[itOuter]],
                                     thetaStan = thetaStan)
    # plot(logL[-1], xlab="iteration", ylab=paste("Initialization logL value for",G) )


    threshold_outer <- 12
    if(itOuter > (threshold_outer + 1)) {

      if (all(coda::heidel.diag(logL[- 1], eps = 0.1, pvalue = 0.05)[, c(1, 4)] == 1) || itOuter > 50) {
        programclust <- vector()
        programclust <- map(z)

        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[! K]) > 0) { # J[!K] tells which are empty clusters
          z <- z[, -J[! K]]
          programclust <- map(z)
        }

        convOuter <- 1
      }
    }


    if(convOuter != 1) { # only update until convergence, not after
      # updating z value
      z <- calcZvalue(PI = PI,
                      dataset = dataset,
                      M = vectorizedM,
                      G = G,
                      Sigma = SigmaAllOuter[[itOuter]],
                      thetaStan = thetaStan,
                      normFactors = normFactors)
      itOuter <- itOuter + 1
      nIterations <- nIterations + 10
    }
  } # end of loop


  results <- list(finalmu = MAllOuter[[itOuter]] +
                    do.call("rbind", rep(list(matrix(normFactors, byrow = TRUE, ncol = p)), G)),
                  finalsigma = SigmaAllOuter[[itOuter]],
                  finalphi = PhiAllOuter[[itOuter]],
                  finalomega = OmegaAllOuter[[itOuter]],
                  allmu = lapply(MAllOuter, function(x) (x + do.call("rbind", rep(list(matrix(normFactors, byrow = TRUE, ncol = p)), G)))),
                  allsigma = SigmaAllOuter,
                  allphi = PhiAllOuter,
                  allomega = OmegaAllOuter,
                  clusterlabels = programclust,
                  iterations = itOuter,
                  FinalRstanIterations = nIterations,
                  proportion = PI,
                  loglikelihood = logL,
                  probaPost = z)

  class(results) <- "mvplnCluster"
  return(results)
  # Developed by Anjali Silva
}

# Stan sampling
stanRun <- function(r,
                    p,
                    dataset,
                    G,
                    nChains,
                    nIterations,
                    mod,
                    Mu,
                    Sigma,
                    normFactors) {
  n <- nrow(dataset)

  fitrstan <- list()
  for (g in seq_along(1:G)) {
    data1 <- list(d = ncol(dataset),
                  N = nrow(dataset),
                  y = dataset,
                  mu = as.vector(Mu[((g - 1) * r + 1):(g * r), ]),
                  Sigma = Sigma[((g - 1) * (r * p) + 1):(g *(r * p)), ],
                  normfactors = as.vector(normFactors))
    stanproceed <- 0
    try <- 1

    while (! stanproceed){

      fitrstan[[g]] <- rstan::sampling(object = mod,
                                       data = data1,
                                       iter = nIterations,
                                       chains = nChains,
                                       verbose = FALSE,
                                       refresh = -1)
      if (all(rstan::summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) == TRUE &&
          all(rstan::summary(fitrstan[[g]])$summary[, "n_eff"] > 100) == TRUE) {
        stanproceed <- 1
      } else if(all(rstan::summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) != TRUE ||
                all(rstan::summary(fitrstan[[g]])$summary[,"n_eff"] > 100) != TRUE) {
        if(try == 5) { # stop after 5 tries
          stanproceed <- 1
        }
        nIterations <- nIterations + 100
        try <- try + 1
      }
    }
  }

  results <- list(fitrstan = fitrstan,
                  nIterations = nIterations)
  class(results) <- "RStan"
  return(results)
  # Developed by Anjali Silva
}

# [END]
