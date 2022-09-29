#' Clustering Using mixtures of MVPLN via variational Gaussian approximations
#'
#' Performs clustering using mixtures of matrix variate Poisson-log normal
#' (MVPLN) via variational Gaussian approximations. Model selection can be done
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
#' # Generating simulated matrix variate count data
#' set.seed(1234)
#' trueG <- 2 # number of total G
#' truer <- 2 # number of total occasions
#' truep <- 3 # number of total responses
#' trueN <- 100 # number of total units
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
#' if (!require(clusterGeneration)) install.packages("clusterGeneration")
#'
#' # Covariance matrix containing variances and covariances between r occasions
#' truePhi1 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#'                                                   dim = truer,
#'                                                   rangeVar = c(1, 1.7))$Sigma
#' truePhi1[1, 1] <- 1 # For identifiability issues
#'
#' truePhi2 <- clusterGeneration::genPositiveDefMat("unifcorrmat",
#'                                                   dim = truer,
#'                                                   rangeVar = c(0.7, 0.7))$Sigma
#' truePhi2[1, 1] <- 1 # For identifiability issues
#' truePhiall <- rbind(truePhi1, truePhi2)
#'
#' # Omega is a p x p matrix
#' # Covariance matrix containing variances and covariances between p responses
#' trueOmega1 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#'                                    rangeVar = c(1, 1.7))$Sigma
#'
#' trueOmega2 <- clusterGeneration::genPositiveDefMat("unifcorrmat", dim = truep,
#'                                    rangeVar = c(0.7, 0.7))$Sigma
#' trueOmegaAll <- rbind(trueOmega1, trueOmega2)
#'
#' # Generated simulated data
#' sampleData <- mixMVPLN::mvplnDataGenerator(nOccasions = truer,
#'                                            nResponses = truep,
#'                                            nUnits = trueN,
#'                                            mixingProportions = c(0.79, 0.21),
#'                                            matrixMean = trueMall,
#'                                            phi = truePhiall,
#'                                            omega = trueOmegaAll)
#'
#' # Clustering simulated matrix variate count data
#' clusteringResults <- mixMVPLN::mvplnVGAclus(dataset = sampleData$dataset,
#'                                       membership = sampleData$truemembership,
#'                                       gmin = 1,
#'                                       gmax = 3,
#'                                       initMethod = "kmeans",
#'                                       nInitIterations = 1,
#'                                       normalize = "Yes")
#' }
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}, Sanjeena Dang,
#'          \email{sanjeenadang@cunet.carleton.ca}. }
#'
#' @references
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


  if(all(membership != "none") && is.numeric(membership) != TRUE) {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") &&
     all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
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




  parallelFA <- function(G, TwoDdataset,
                         r, p, d, N,
                         normalize) {

    ###### Parameter Updates ####

    # Calculating normalization factors
    if(normalize == "Yes") {
      normFactors <- as.vector(edgeR::calcNormFactors(as.matrix(TwoDdataset),
                                                      method = "TMM"))

    } else if(normalize == "No") {
      normFactors <- rep(0, d)
    } else {
      stop("Argument normalize should be 'Yes' or 'No' ")
    }


    lib_mat <- matrix(normFactors,N,d,byrow=T)
    lib_mat_list <- list()
    for (i in 1:N){
      lib_mat_list[[i]]<-t(matrix(lib_mat[i,],nrow=p))
    }

    Y_list <- list()
    for (i in 1:N){
      Y_list[[i]]<-t(matrix(TwoDdataset[i,],nrow=p))
    }






    #### Initialization ###
    mu <- list() ###This is M
    #psi <- list() #### This is kronecker of omega and Phi
    omega <- list() ### This is omega
    phi <- list() ### This is Phi
    delta <- list() #### This is variational parameter Delta
    kappa <- list() #### This is variational parameter kappa
    sigma <- list()  #### This is Psi in the math - kronecker of omega and Phi
    isigma <- list() ###Inverse of Psi
    iphi <- list() ### This is inverse of Phi
    iomega <- list() ### This is inverse of omega
    m <- list()### Vectorized xi
    S <- list() ###It is Psi
    # P <- list()
    # Q <- list()

    ###Other intermediate items initialized
    start <- list()
    Sk <- array(0, c(d,d,G) )
    GX <- list()
    dGX <- list()
    z_S <- list()

    i_kappa <- list() ### inverse of kappa
    i_delta <- list() ### inverse of Delta
    start_list <- list()

    k_means <- kmeans(log(TwoDdataset+1),
                      centers = G,
                      nstart = 50, iter.max = 20)$cluster
    z <- mclust::unmap(k_means) ### Starting value for Z
    pi_g <- colSums(z) / N

    for (g in 1:G) {
      obs <- which(z[, g] ==1 )
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
      start_list[[g]] <- list()
      for (i in 1:N) {
        start_list[[g]][[i]] <- log(Y_list[[i]])
        delta[[g]][[i]] <- diag(r)*0.001
        kappa[[g]][[i]] <- diag(p)*0.001
        S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
      }
    }



    checks <- 0
    it <- 1
    aloglik <- loglik <- NULL
    aloglik[c(1:5)] <- 0
    it_max <- 200


    while (checks==0) {


      for (g in 1:G) {
        GX[[g]] <- list()
        dGX[[g]] <- list()
        z_S[[g]] <- list()
        i_delta[[g]] <- list()
        i_kappa[[g]] <- list()
        delta_o <- delta[[g]]
        kappa_o <- kappa[[g]]

        for (i in 1:N) {
          i_delta[[g]][[i]] <- diag(c(t(diag(kappa[[g]][[i]])) %*%
                               t(exp(log(lib_mat_list[[i]]) +
                               start_list[[g]][[i]] +
                               0.5 * diag(delta[[g]][[i]]) %*%
                               t(diag(kappa[[g]][[i]])))))) +
                               iphi[[g]] * sum(diag(iomega[[g]] %*%
                               kappa[[g]][[i]]))
          delta[[g]][[i]] <- p * solve(i_delta[[g]][[i]])


          i_kappa[[g]][[i]] <- diag(c(t(diag(delta[[g]][[i]])) %*%
                               (exp(log(lib_mat_list[[i]]) +
                               start_list[[g]][[i]] +
                               t(0.5 * diag(kappa[[g]][[i]]) %*%
                               t(diag(delta[[g]][[i]]))))))) +
                               iomega[[g]] * sum(diag(iphi[[g]] %*%
                               delta[[g]][[i]]))


          kappa[[g]][[i]] <- solve(i_kappa[[g]][[i]])


          S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
          z_S[[g]][[i]] <- z[i, g] * S[[g]][[i]]
          GX[[g]][[i]] <- TwoDdataset[i, ] -
                          exp(start[[g]][i, ] +
                          log(lib_mat[i,]) +
                          0.5 * diag(S[[g]][[i]])) -
                          (isigma[[g]]) %*% (start[[g]][i, ] - mu[[g]])
          m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
          start_list[[g]][[i]] <- t(matrix(m[[g]][i, ], nrow = p))
        }
        start[[g]] <- m[[g]]

        mu[[g]] <- colSums(z[, g] * m[[g]]) / sum(z[, g]) # this is xi



        # Updating Sample covariance
        mu_mat <- t(matrix(mu[[g]], nrow = p))
        phi_m <- list()
        for (i in 1:N) {
          phi_m[[i]] <- z[i,g]*(start_list[[g]][[i]] - mu_mat) %*%
                        iomega[[g]] %*% t(start_list[[g]][[i]] - mu_mat) +
                        z[i, g] * delta[[g]][[i]] * sum(diag(iomega[[g]] %*%
                        kappa[[g]][[i]]))
        }
        phi[[g]] <- Reduce("+", phi_m) / sum(z[, g]*p)
        iphi[[g]] <- solve(phi[[g]])

        omega_m <- list()
        for (i in 1:N) {
          omega_m[[i]] <- z[i, g] * t(start_list[[g]][[i]] - mu_mat) %*%
                          iphi[[g]] %*% (start_list[[g]][[i]] - mu_mat) +
                          z[i, g] * kappa[[g]][[i]] * sum(diag(iphi[[g]] %*%
                          delta[[g]][[i]]))
        }
        omega[[g]] <- Reduce("+", omega_m) / sum(z[, g] * r)
        iomega[[g]] <- solve(omega[[g]])
        sigma[[g]] <- phi[[g]] %x% omega[[g]]
        isigma[[g]] <- iphi[[g]] %x% iomega[[g]]
      }





      pi_g <- colSums(z) / N
      ### Some useful functions
      fun_five <- function(x, y = isigma[[g]]) {
        sum(diag(x %*% y))
      }

      FMatrix  <- matrix(NA, ncol = G, nrow = N)

      for (g in 1:G) {
        two <- rowSums(exp(m[[g]] + log(lib_mat) +
               0.5 * matrix(unlist(lapply(S[[g]], diag)),
               ncol = d, byrow = TRUE)))
        five <- 0.5 * unlist(lapply(S[[g]], fun_five))
        six <- 0.5 * log(unlist(lapply(S[[g]], det)))
        FMatrix [, g] <- pi_g[g] * exp(rowSums(m[[g]] * TwoDdataset) -
                  two - rowSums(lfactorial(TwoDdataset)) +
                  rowSums(log(lib_mat) * TwoDdataset) - 0.5 *
                  mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]],
                  inverted = TRUE) - five + six + 0.5 *
                  log(det(isigma[[g]])) - d/2)
      }

      loglik[it] <- sum(log(rowSums(FMatrix )))

      z <- FMatrix  / rowSums(FMatrix )
      if (it <= 5) {
        z[z == "NaN"] <- 0
      }

      if (it > 5) {
        #Aitkaine's stopping criterion
        if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
          a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
          add_to <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
          aloglik[it] <- loglik[it - 1] + add_to
          if (abs(aloglik[it] - aloglik[it - 1]) < 0.05) checks <- 1 else checks <- checks
        }
      }
      #print(it)
      it <- it + 1
      if (it == it_max) checks <- 1
      final_phi <- list()
      final_omega <- list()
      for (g in 1:G) {
        final_phi[[g]] <- phi[[g]] / diag(phi[[g]])[1]
        final_omega[[g]] <- omega[[g]] * diag(phi[[g]])[1]
      }
    }


    return(list(pi_g = pi_g,
                mu = mu,
                sigma = sigma,
                z = z,
                loglik = loglik,
                kmeans = k_means,
                phi = final_phi,
                omega = final_omega))
  }

  parallelFA(G = 2,
              TwoDdataset = bozzo2016,
              r=r, # variety
              p= p, # growth stages
              d = d,
              N = n,
              normalize = normalize)



    for(g in seq_along(1:(gmax - gmin + 1))) {
      parallel_FA(G = g,
                  TwoDdataset = TwoDdataset,
                  r = r,
                  p = p,
                  n = n,
                  d = d,
                  normFactors = normFactors)
    }


    # cluster data result extracting
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
                    mu = mu,
                    sigma = sigma,
                    z = z,
                    kmeans = kMeans,
                    phi = finalPhi,
                    omega = finalOmega,
                    loglikelihood = loglik,
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
