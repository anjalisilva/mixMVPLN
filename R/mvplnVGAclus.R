#' Clustering Using mixtures of MVPLN via variational Gaussian approximations
#'
#' Performs clustering using mixtures of matrix variate Poisson-log normal
#' (MVPLN) via variational Gaussian approximations. Model selection can
#' be done using AIC, AIC3, BIC and ICL.
#'
#' @param dataset A list of length nUnits, containing Y_j matrices.
#'    A matrix Y_j has size r x p, and the dataset will have 'j' such
#'    matrices with j = 1,...,n. If a Y_j has all zeros, such Y_j will
#'    be removed prior to cluster analysis.
#' @param membership A numeric vector of size length(dataset) containing
#'    the cluster membership of each Y_j matrix. If not available, leave
#'    as "none".
#' @param gmin A positive integer specifying the minimum number of
#'    components to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each MCMC chain (including warmup). The value should be greater
#'    than 40. The upper limit will depend on size of dataset.
#' @param initMethod A method for initialization. Current options are
#'    "none", "kmeans", "random", "medoids", "clara", or "fanny".
#'    Default is "none".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be considered. Default is 0.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization
#'     factors are calculated using TMM method of edgeR package.
#'     Default is "Yes".
#'
#' @return Returns an S3 object of class mvplnParallel with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item nUnits - Number of units in the input dataset.
#'   \item nVariables - Number of variables in the input dataset.
#'   \item nOccassions - Number of occassions in the input dataset.
#'   \item normFactors - A vector of normalization factors used for
#'      input dataset.
#'   \item gmin - Minimum number of components considered in the
#'      clustering run.
#'   \item gmax - Maximum number of components considered in the
#'      clustering run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item loglikelihood - A vector with value of final log-likelihoods
#'      for each cluster size.
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
#' trueN <- 1000 # number of total units
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
#' # Clustering simulated matrix variate count data
#' clusteringResults <- mixMVPLN::mvplnVGAclus(
#'                       dataset = sampleData$dataset,
#'                       membership = sampleData$truemembership,
#'                       gmin = 1,
#'                       gmax = 3,
#'                       initMethod = "kmeans",
#'                       nInitIterations = 2,
#'                       normalize = "Yes")
#'
#' # Example 2
#' trueG <- 1 # 1 cluster
#' truer <- 2 # variety
#' truep <- 3 # growth stages
#' trueN <- 1000 # genes
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
#' clusteringResults <- mixMVPLN::mvplnVGAclus(
#'                       dataset = sampleData2$dataset,
#'                       membership = sampleData2$truemembership,
#'                       gmin = 1,
#'                       gmax = 2,
#'                       initMethod = "kmeans",
#'                       nInitIterations = 2,
#'                       normalize = "Yes")
#'
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
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils tail
#'
mvplnVGAclus <- function(dataset,
                          membership = "none",
                          gmin,
                          gmax,
                          initMethod = "kmeans",
                          nInitIterations = 0,
                          normalize = "Yes") {

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
    if(initMethod != "none" & initMethod != "kmeans" & initMethod != "random" & initMethod != "medoids" & initMethod != "medoids" & initMethod != "clara" & initMethod != "fanny") {
      stop("initMethod should of class character, specifying
      either: none, kmeans, random, medoids, clara, or fanny.")
    }
  } else if (is.character(initMethod) != TRUE) {
    stop("initMethod should of class character, specifying
      either: none, kmeans, random, medoids, clara, or fanny.")
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


  n <- length(dataset)
  p <- ncol(dataset[[1]])
  r <- nrow(dataset[[1]])
  d <- p * r


  # Changing dataset into a n x rp dataset
  TwoDdataset <- matrix(NA, ncol = d, nrow = n)
  sampleMatrix <- matrix(c(1:d), nrow = r, byrow = TRUE)
  for (u in 1:n) {
    for (e in 1:p) {
      for (s in 1:r) {
        TwoDdataset[u, sampleMatrix[s, e]] <- dataset[[u]][s, e]
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
    normFactors <- as.vector(edgeR::calcNormFactors(as.matrix(TwoDdataset),
                                                    method = "TMM"))

  } else if(normalize == "No") {
    normFactors <- rep(0, d)
  } else {
    stop("Argument normalize should be 'Yes' or 'No' ")
  }

  parallelFA <- function(G, dataset,
                         TwoDdataset,
                         r, p, d, N,
                         normFactors,
                         nInitIterations,
                         initMethod) {

    # arranging normalization factors
    libMat <- matrix(normFactors, N, d, byrow = T)
    libMatList <- list()
    for (i in 1:N) {
      libMatList[[i]] <- t(matrix(libMat[i, ], nrow = p))
    }

    if (initMethod == "none") {

      # Initialization
      mu <- omega <- phi <- list() # mu is M;
      delta <- kappa <- sigma <- isigma <- iphi <- iomega <- list()
      # delta  is variational parameter Delta
      # kappa is variational parameter kappa
      # sigma is Psi in the math - kronecker of omega and Phi
      # isigma is inverse of Psi
      # iphi is inverse of Phi
      # iomega is inverse of omega

      m <- S <- list()
      # m is vectorized xi
      # S is Psi

      # Other intermediate items initialized
      Sk <- array(0, c(d, d, G) )
      start <- GX <- dGX <- zS <- list()

      iKappa <- iDelta <- startList <- list()
      # iKappa is inverse of kappa
      # iDelta is inverse of Delta

      kMeansResults <- kmeans(log(TwoDdataset+1),
                              centers = G,
                              nstart = 50,
                              iter.max = 20)$cluster
      zValue <- mclust::unmap(kMeansResults) ### Starting value for Z
      piG <- colSums(zValue) / N

      for (g in 1:G) {
        obs <- which(zValue[, g] == 1)
        mu[[g]] <- colMeans(log(TwoDdataset[obs, ] + 1 / 6))
        sigma[[g]] <- var(log(TwoDdataset[obs, ] + 1 / 6))
        isigma[[g]] <- solve(sigma[[g]])
        phi[[g]] <- diag(r) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        omega[[g]] <- diag(p) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        iphi[[g]] <- solve(phi[[g]])
        iomega[[g]] <- solve(omega[[g]])
      }

      for (g in 1:G) {
        start[[g]] <- log(TwoDdataset + 1/6) ###Starting value for M
        m[[g]] <- log(TwoDdataset + 1/6)
        S[[g]] <- list()
        delta[[g]] <- list()
        kappa[[g]] <- list()
        startList[[g]] <- list()
        for (i in 1:N) {
          startList[[g]][[i]] <- log(dataset[[i]])
          delta[[g]][[i]] <- diag(r) * 0.001
          kappa[[g]][[i]] <- diag(p) * 0.001
          S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
        }
      }
    } else {

    # Initialize based on specified initMethod method
    outputInitialization <- list()
    checklogL <- vector()
    for (initIt in 1:nInitIterations) {
       set.seed(initIt)
       outputInitialization[[initIt]] <- initializationRun(
                                                G,
                                                dataset,
                                                TwoDdataset,
                                                r, p, d, N,
                                                normFactors,
                                                nInitIterations,
                                                initMethod)
       checklogL[initIt] <- outputInitialization[[initIt]]$finalLogLik
    }

    # select init run with highest logL
    maxRun <- which.max(checklogL)

    mu <- outputInitialization[[maxRun]]$mu
    omega <- outputInitialization[[maxRun]]$omega
    phi <- outputInitialization[[maxRun]]$phi
    delta <- outputInitialization[[maxRun]]$delta
    kappa <- outputInitialization[[maxRun]]$kappa
    sigma <- outputInitialization[[maxRun]]$sigma
    isigma <- outputInitialization[[maxRun]]$isigma
    iphi <- outputInitialization[[maxRun]]$iphi
    iomega <- outputInitialization[[maxRun]]$iomega
    m <- outputInitialization[[maxRun]]$m
    S <- outputInitialization[[maxRun]]$S
    Sk <- outputInitialization[[maxRun]]$Sk
    start <- outputInitialization[[maxRun]]$start
    GX <- outputInitialization[[maxRun]]$GX
    dGX <- outputInitialization[[maxRun]]$dGX
    zS <- outputInitialization[[maxRun]]$zS
    iKappa <- outputInitialization[[maxRun]]$iKappa
    iDelta <- outputInitialization[[maxRun]]$iDelta
    startList <- outputInitialization[[maxRun]]$startList
    zValue <- outputInitialization[[maxRun]]$zValue
    piG <- outputInitialization[[maxRun]]$piG
  }

    # start clustering after initialization
    it <- 1
    aloglik <- loglik <- NULL
    checks <- aloglik[c(1:5)] <- 0
    itMax <- 200

    while (checks == 0) {

      for (g in 1:G) {
        GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
        iDelta[[g]] <- iKappa[[g]] <- list()
        deltaO <- delta[[g]]
        kappaO <- kappa[[g]]

        for (i in 1:N) {
          iDelta[[g]][[i]] <- diag(c(t(diag(kappa[[g]][[i]])) %*%
                               t(exp(log(libMatList[[i]]) +
                               startList[[g]][[i]] +
                               0.5 * diag(delta[[g]][[i]]) %*%
                               t(diag(kappa[[g]][[i]])))))) +
                               iphi[[g]] * sum(diag(iomega[[g]] %*%
                               kappa[[g]][[i]]))
          delta[[g]][[i]] <- p * solve(iDelta[[g]][[i]])


          iKappa[[g]][[i]] <- diag(c(t(diag(delta[[g]][[i]])) %*%
                               (exp(log(libMatList[[i]]) +
                               startList[[g]][[i]] +
                               t(0.5 * diag(kappa[[g]][[i]]) %*%
                               t(diag(delta[[g]][[i]]))))))) +
                               iomega[[g]] * sum(diag(iphi[[g]] %*%
                               delta[[g]][[i]]))

          kappa[[g]][[i]] <- solve(iKappa[[g]][[i]])


          S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
          zS[[g]][[i]] <- zValue[i, g] * S[[g]][[i]]
          GX[[g]][[i]] <- TwoDdataset[i, ] -
                          exp(start[[g]][i, ] +
                          log(libMat[i,]) +
                          0.5 * diag(S[[g]][[i]])) -
                          (isigma[[g]]) %*% (start[[g]][i, ] - mu[[g]])
          m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
          startList[[g]][[i]] <- t(matrix(m[[g]][i, ], nrow = p))
        }
        start[[g]] <- m[[g]]

        mu[[g]] <- colSums(zValue[, g] * m[[g]]) / sum(zValue[, g]) # this is xi

        # Updating Sample covariance
        muMat <- t(matrix(mu[[g]], nrow = p))
        phiM <- list()
        for (i in 1:N) {
          phiM[[i]] <- zValue[i,g]*(startList[[g]][[i]] - muMat) %*%
                        iomega[[g]] %*% t(startList[[g]][[i]] - muMat) +
                        zValue[i, g] * delta[[g]][[i]] * sum(diag(iomega[[g]] %*%
                        kappa[[g]][[i]]))
        }
        phi[[g]] <- Reduce("+", phiM) / sum(zValue[, g] * p)
        iphi[[g]] <- solve(phi[[g]])

        omegaM <- list()
        for (i in 1:N) {
          omegaM[[i]] <- zValue[i, g] * t(startList[[g]][[i]] - muMat) %*%
                          iphi[[g]] %*% (startList[[g]][[i]] - muMat) +
                          zValue[i, g] * kappa[[g]][[i]] * sum(diag(iphi[[g]] %*%
                          delta[[g]][[i]]))
        }
        omega[[g]] <- Reduce("+", omegaM) / sum(zValue[, g] * r)
        iomega[[g]] <- solve(omega[[g]])
        sigma[[g]] <- phi[[g]] %x% omega[[g]]
        isigma[[g]] <- iphi[[g]] %x% iomega[[g]]
      }

      piG <- colSums(zValue) / N
      # Internal functions
      funFive <- function(x, y = isigma[[g]]) {
        sum(diag(x %*% y))
      }

      FMatrix <- matrix(NA, ncol = G, nrow = N)

      for (g in 1:G) {
        two <- rowSums(exp(m[[g]] + log(libMat) +
               0.5 * matrix(unlist(lapply(S[[g]], diag)),
               ncol = d, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], funFive))
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        FMatrix[, g] <- piG[g] * exp(rowSums(m[[g]] * TwoDdataset) -
                  two - rowSums(lfactorial(TwoDdataset)) +
                  rowSums(log(libMat) * TwoDdataset) - 0.5 *
                  mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]],
                  inverted = TRUE) - five + six + 0.5 *
                  log(det(isigma[[g]])) - d/2)
      }

      loglik[it] <- sum(log(rowSums(FMatrix)))

      # Plotting
      # if (it > 3) {
      #   plot(loglik, type = "l")
      # }

      zValue <- FMatrix / rowSums(FMatrix)
      if (it <= 5) {
        zValue[zValue == "NaN"] <- 0
      }

      if (it > 5) {
        # Aitkaine's stopping criterion
        # Print for check; commented out
        # cat("\n At it:", it ,"loglik[it - 1]:",
        # loglik[it - 1], "- loglik[it - 2]:", loglik[it - 2],
        # "is", loglik[it - 1] - loglik[it - 2], "\n")
        if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
          a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
          addTo <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
          aloglik[it] <- loglik[it - 1] + addTo
          if (abs(aloglik[it] - aloglik[it - 1]) < 0.05) {
            checks <- 1
          } else {
            checks <- checks
          }
        }
      }

      it <- it + 1
      if (it == itMax) {
        checks <- 1
      }
      finalPhi <- finalOmega <- list()
      for (g in 1:G) {
        finalPhi[[g]] <- phi[[g]] / diag(phi[[g]])[1]
        finalOmega[[g]] <- omega[[g]] * diag(phi[[g]])[1]
      }
    }

    programclust <- mclust::map(zValue)

    FinalGResults <- list(mu = mu,
                          sigma = sigma,
                          phi = finalPhi,
                          omega = finalOmega,
                          probaPost = zValue,
                          loglikelihood = loglik,
                          proportion = piG,
                          clusterlabels = programclust,
                          iterations = it)

    class(FinalGResults) <- "FinalGResults"

    return(FinalGResults)
  }

  parallelFAOutput <- list()
  for(g in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }

    parallelFAOutput[[g]] <- parallelFA(
                 G = clustersize,
                 dataset = dataset,
                 TwoDdataset = TwoDdataset,
                 r = r,
                 p = p,
                 d = d,
                 N = n,
                 normFactors = normFactors,
                 nInitIterations = nInitIterations,
                 initMethod = initMethod)
    }


    # cluster data result extracting
    BIC <- ICL <- AIC <- AIC3 <- k <- ll <- vector()
    for(g in seq_along(1:(gmax - gmin + 1))) {
      #print(g)
      if(length(1:(gmax - gmin + 1)) == gmax) {
        clustersize <- g
      } else if(length(1:(gmax - gmin + 1)) < gmax) {
        clustersize <- seq(gmin, gmax, 1)[g]
      }

      # save the final log-likelihood
      ll[g] <- unlist(tail(parallelFAOutput[[g]]$loglik, n = 1))

      k[g] <- calcParameters(g = clustersize,
                             r = r,
                             p = p)

      # starting model selection
      if (g == max(1:(gmax - gmin + 1))) {
        bic <- BICFunction(ll = ll,
                           k = k,
                           n = n,
                           run = parallelFAOutput,
                           gmin = gmin,
                           gmax = gmax,
                           parallel = FALSE)

        icl <- ICLFunction(bIc = bic,
                           gmin = gmin,
                           gmax = gmax,
                           run = parallelFAOutput,
                           parallel = FALSE)

        aic <- AICFunction(ll = ll,
                           k = k,
                           run = parallelFAOutput,
                           gmin = gmin,
                           gmax = gmax,
                           parallel = FALSE)

        aic3 <- AIC3Function(ll = ll,
                             k = k,
                             run = parallelFAOutput,
                             gmin = gmin,
                             gmax = gmax,
                             parallel = FALSE)
      }
    }



    final <- base::proc.time() - ptm

    RESULTS <- list(dataset = dataset,
                    nUnits = n,
                    nVariables = p,
                    nOccassions = r,
                    normFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = parallelFAOutput,
                    loglikelihood = ll,
                    nParameters = k,
                    trueLabels = membership,
                    ICLAll = icl,
                    BICAll = bic,
                    AICAll = aic,
                    AIC3All = aic3,
                    totalTime = final)

    class(RESULTS) <- "mvplnParallel"
    return(RESULTS)

}

initializationRun <- function(G,
                              dataset,
                              TwoDdataset,
                              r, p, d, N,
                              normFactors,
                              nInitIterations,
                              initMethod) {

    # arranging normalization factors
    libMat <- matrix(normFactors, N, d, byrow = T)
    libMatList <- list()
    for (i in 1:N) {
      libMatList[[i]] <- t(matrix(libMat[i, ], nrow = p))
    }

    # Initialization
      mu <- omega <- phi <- list() # mu is M;
      delta <- kappa <- sigma <- isigma <- iphi <- iomega <- list()
      # delta  is variational parameter Delta
      # kappa is variational parameter kappa
      # sigma is Psi in the math - kronecker of omega and Phi
      # isigma is inverse of Psi
      # iphi is inverse of Phi
      # iomega is inverse of omega

      m <- S <- list()
      # m is vectorized xi
      # S is Psi

      # Other intermediate items initialized
      Sk <- array(0, c(d, d, G) )
      start <- GX <- dGX <- zS <- list()

      iKappa <- iDelta <- startList <- list()
      # iKappa is inverse of kappa
      # iDelta is inverse of Delta


      if (initMethod == "kmeans") {
        # cat("\n initMethod == kmeans \n")
        zValue <- mclust::unmap(stats::kmeans(x = log(TwoDdataset + 1 / 3),
                                              centers = G)$cluster)

      } else if (initMethod == "random") {
        # cat("\n initMethod == random \n")
        if(G == 1) { # generating z if g = 1
          zValue <- as.matrix(rep.int(1, times = n),
                                       ncol = G,
                                       nrow = n)
        } else { # generating z if g>1
          zConv = 0
          while(! zConv) { # ensure that dimension of z is same as G (i.e.
            # if one column contains all 0s, then generate z again)
            zValue <- t(stats::rmultinom(n = n,
                                         size = 1,
                                         prob = rep(1 / G, G)))
            if(length(which(colSums(zValue) > 0)) == G) {
              zConv <- 1
            }
          }
        }
      } else if (initMethod == "medoids") {
        # cat("\n initMethod == medoids \n")
        zValue <- mclust::unmap(cluster::pam(x = log(TwoDdataset + 1 / 3),
                                                      k = G)$cluster)
      } else if (initMethod == "clara") {
        # cat("\n initMethod == clara \n")
        zValue <- mclust::unmap(cluster::clara(x = log(TwoDdataset + 1 / 3),
                                                        k = G)$cluster)
      } else if (initMethod == "fanny") {
        # cat("\n initMethod == fanny \n")
        zValue <- mclust::unmap(cluster::fanny(x = log(TwoDdataset + 1 / 3),
                                                        k = G)$cluster)
      }

      piG <- colSums(zValue) / N

      for (g in 1:G) {
        obs <- which(zValue[, g] == 1)
        mu[[g]] <- colMeans(log(TwoDdataset[obs, ] + 1 / 6))
        sigma[[g]] <- var(log(TwoDdataset[obs, ] + 1 / 6))
        isigma[[g]] <- solve(sigma[[g]])
        phi[[g]] <- diag(r) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        omega[[g]] <- diag(p) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        iphi[[g]] <- solve(phi[[g]])
        iomega[[g]] <- solve(omega[[g]])
      }

      for (g in 1:G) {
        start[[g]] <- log(TwoDdataset + 1/6) ###Starting value for M
        m[[g]] <- log(TwoDdataset + 1/6)
        S[[g]] <- list()
        delta[[g]] <- list()
        kappa[[g]] <- list()
        startList[[g]] <- list()
        for (i in 1:N) {
          startList[[g]][[i]] <- log(dataset[[i]])
          delta[[g]][[i]] <- diag(r) * 0.001
          kappa[[g]][[i]] <- diag(p) * 0.001
          S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
        }
      }


    it <- 1
    aloglik <- loglik <- NULL
    checks <- aloglik[c(1:5)] <- 0
    itMax <- 200

    while (checks == 0) {

      for (g in 1:G) {
        GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
        iDelta[[g]] <- iKappa[[g]] <- list()
        deltaO <- delta[[g]]
        kappaO <- kappa[[g]]

        for (i in 1:N) {
          iDelta[[g]][[i]] <- diag(c(t(diag(kappa[[g]][[i]])) %*%
                                       t(exp(log(libMatList[[i]]) +
                                               startList[[g]][[i]] +
                                               0.5 * diag(delta[[g]][[i]]) %*%
                                               t(diag(kappa[[g]][[i]])))))) +
            iphi[[g]] * sum(diag(iomega[[g]] %*%
                                   kappa[[g]][[i]]))
          delta[[g]][[i]] <- p * solve(iDelta[[g]][[i]])


          iKappa[[g]][[i]] <- diag(c(t(diag(delta[[g]][[i]])) %*%
                                       (exp(log(libMatList[[i]]) +
                                              startList[[g]][[i]] +
                                              t(0.5 * diag(kappa[[g]][[i]]) %*%
                                                  t(diag(delta[[g]][[i]]))))))) +
            iomega[[g]] * sum(diag(iphi[[g]] %*%
                                     delta[[g]][[i]]))

          kappa[[g]][[i]] <- solve(iKappa[[g]][[i]])


          S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
          zS[[g]][[i]] <- zValue[i, g] * S[[g]][[i]]
          GX[[g]][[i]] <- TwoDdataset[i, ] -
            exp(start[[g]][i, ] +
                  log(libMat[i,]) +
                  0.5 * diag(S[[g]][[i]])) -
            (isigma[[g]]) %*% (start[[g]][i, ] - mu[[g]])
          m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
          startList[[g]][[i]] <- t(matrix(m[[g]][i, ], nrow = p))
        }
        start[[g]] <- m[[g]]

        mu[[g]] <- colSums(zValue[, g] * m[[g]]) / sum(zValue[, g]) # this is xi

        # Updating Sample covariance
        muMat <- t(matrix(mu[[g]], nrow = p))
        phiM <- list()
        for (i in 1:N) {
          phiM[[i]] <- zValue[i,g]*(startList[[g]][[i]] - muMat) %*%
            iomega[[g]] %*% t(startList[[g]][[i]] - muMat) +
            zValue[i, g] * delta[[g]][[i]] * sum(diag(iomega[[g]] %*%
                                                        kappa[[g]][[i]]))
        }
        phi[[g]] <- Reduce("+", phiM) / sum(zValue[, g] * p)
        iphi[[g]] <- solve(phi[[g]])

        omegaM <- list()
        for (i in 1:N) {
          omegaM[[i]] <- zValue[i, g] * t(startList[[g]][[i]] - muMat) %*%
            iphi[[g]] %*% (startList[[g]][[i]] - muMat) +
            zValue[i, g] * kappa[[g]][[i]] * sum(diag(iphi[[g]] %*%
                                                        delta[[g]][[i]]))
        }
        omega[[g]] <- Reduce("+", omegaM) / sum(zValue[, g] * r)
        iomega[[g]] <- solve(omega[[g]])
        sigma[[g]] <- phi[[g]] %x% omega[[g]]
        isigma[[g]] <- iphi[[g]] %x% iomega[[g]]
      }

      piG <- colSums(zValue) / N
      # Internal functions
      funFive <- function(x, y = isigma[[g]]) {
        sum(diag(x %*% y))
      }

      FMatrix <- matrix(NA, ncol = G, nrow = N)

      for (g in 1:G) {
        two <- rowSums(exp(m[[g]] + log(libMat) +
                             0.5 * matrix(unlist(lapply(S[[g]], diag)),
                                          ncol = d, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], funFive))
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        FMatrix[, g] <- piG[g] * exp(rowSums(m[[g]] * TwoDdataset) -
                                       two - rowSums(lfactorial(TwoDdataset)) +
                                       rowSums(log(libMat) * TwoDdataset) - 0.5 *
                                       mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]],
                                                   inverted = TRUE) - five + six + 0.5 *
                                       log(det(isigma[[g]])) - d/2)
      }

      loglik[it] <- sum(log(rowSums(FMatrix)))

      zValue <- FMatrix / rowSums(FMatrix)
      if (it <= 5) {
        zValue[zValue == "NaN"] <- 0
      }

      if (it > 5) {
        # Aitkaine's stopping criterion
        if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
          a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
          addTo <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
          aloglik[it] <- loglik[it - 1] + addTo
          if (abs(aloglik[it] - aloglik[it - 1]) < 0.05) {
            checks <- 1
          } else {
            checks <- checks
          }
        }
      }

      it <- it + 1
      if (it == itMax) {
        checks <- 1
      }
      finalPhi <- finalOmega <- list()
      for (g in 1:G) {
        finalPhi[[g]] <- phi[[g]] / diag(phi[[g]])[1]
        finalOmega[[g]] <- omega[[g]] * diag(phi[[g]])[1]
      }
    }

    programclust <- mclust::map(zValue)

    initRunOutput <- list(mu = mu,
                          sigma = sigma,
                          omega = omega,
                          phi = phi,
                          omega = omega,
                          delta = delta,
                          kappa = kappa,
                          isigma = isigma,
                          iphi = iphi,
                          iomega = iomega,
                          m = m,
                          S = S,
                          Sk = Sk,
                          start = start,
                          GX = GX,
                          dGX = dGX,
                          zS = zS,
                          iKappa = iKappa,
                          iDelta = iDelta,
                          startList = startList,
                          zValue = zValue,
                          loglik = loglik,
                          finalLogLik = loglik[it-1],
                          piG = piG,
                          clusterlabels = programclust,
                          iterations = it)

    class(initRunOutput) <- "initializationRun"

    return(initRunOutput)
  }


# [END]

