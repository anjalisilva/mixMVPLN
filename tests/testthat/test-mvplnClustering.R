context("Checking for parallel clustering performance")
library(mixMVPLN)

test_that("Checking clustering results", {

  set.seed(1234)
  true_G <- 2 # number of total G
  true_r <- 2 # number of total occasions
  true_p <- 3 # number of total responses
  true_n <- 50 # number of total units

  # M is a r x p matrix
  true_M1 <- matrix(rep(6, (true_r * true_p)),
                    ncol = true_p, nrow = true_r,
                    byrow = TRUE)
  true_M2 <- matrix(rep(1, (true_r * true_p)),
                    ncol = true_p, nrow = true_r,
                    byrow = TRUE)
  true_M_all <- rbind(true_M1, true_M2)

  # Phi is a r x r matrix
  # Covariance matrix containing variances and covariances between r occasions
  # true_Phi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = true_r,
  #                                                   rangeVar = c(1, 1.7))$Sigma
  true_Phi1 <- matrix(c(1.435610, -1.105615, -1.105615,  1.426492), nrow = 2)
  true_Phi1[1, 1] <- 1 # for identifiability issues

  # true_Phi2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = true_r,
  #                                                   rangeVar = c(0.7, 0.7))$Sigma
  true_Phi2 <- matrix(c(0.7000000, 0.1727312, 0.1727312, 0.7000000), nrow = 2)
  true_Phi2[1, 1] <- 1 # for identifiability issues
  true_Phi_all <- rbind(true_Phi1, true_Phi2)

  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  # true_Omega1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = true_p,
  #                                                     rangeVar = c(1, 1.7))$Sigma
  true_Omega1 <- matrix(c(1.4855139,  0.9049113, -0.9063446,
                          0.9049113,  1.3814824, -1.2298301,
                          -0.9063446, -1.2298301,  1.1979135), nrow = 3)

  # true_Omega2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = true_p,
  #                                                     rangeVar = c(0.7, 0.7))$Sigma
  true_Omega2 <- matrix(c(0.7000000, 0.5379098, 0.4837924,
                          0.5379098, 0.7000000, 0.4089374,
                          0.4837924, 0.4089374, 0.7000000), nrow = 3)
  true_Omega_all <- rbind(true_Omega1, true_Omega2)

  simulatedMVdata <- mvplnDataGenerator(nOccasions = true_r,
                                        nResponses = true_p,
                                        nUnits = true_n,
                                        mixingProportions = c(0.79, 0.21),
                                        matrixMean = true_M_all,
                                        phi = true_Phi_all,
                                        omega = true_Omega_all)

  # Clustering simulated matrix variate count data
  clusteringResults <- mvplnClustering(dataset = simulatedMVdata$dataset,
                                       membership = simulatedMVdata$truemembership,
                                       gmin = 1,
                                       gmax = 1,
                                       nChains = 3,
                                       nIterations = 250,
                                       initMethod = "kmeans",
                                       nInitIterations = 0,
                                       normalize = "Yes",
                                       numNodes = 2)

  # Setting numNodes = 2 based on the following entry, otherwise error.
  # "NB: you can’t use unexported functions and you shouldn’t open new graphics
  # devices or use more than two cores. Individual examples shouldn’t
  # take more than 5s."
  # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa

  testthat::expect_that(length(clusteringResults), equals(17))
  testthat::expect_that(clusteringResults, is_a("mvplnParallel"))
  testthat::expect_that(clusteringResults$initalizationMethod, equals("kmeans"))
  testthat::expect_that(clusteringResults$true_labels, equals(simulatedMVdata$trueMembership))
  testthat::expect_that(trunc(clusteringResults$ICL.all$ICLmodelselected), equals(1))
  testthat::expect_that(trunc(clusteringResults$AIC.all$AICmodelselected), equals(1))
  testthat::expect_that(trunc(clusteringResults$BIC.all$BICmodelselected), equals(1))
})

context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
  set.seed(1234)
  true_G <- 2 # number of total G
  true_r <- 2 # number of total occasions
  true_p <- 3 # number of total responses
  true_n <- 70 # number of total units

  # M is a r x p matrix
  true_M1 <- matrix(rep(6, (true_r * true_p)),
                    ncol = true_p, nrow = true_r,
                    byrow = TRUE)
  true_M2 <- matrix(rep(1, (true_r * true_p)),
                    ncol = true_p, nrow = true_r,
                    byrow = TRUE)
  true_M_all <- rbind(true_M1, true_M2)

  # Phi is a r x r matrix
  library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  # true_Phi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = true_r,
  #                                                   rangeVar = c(1, 1.7))$Sigma
  true_Phi1 <- matrix(c(1.435610, -1.105615, -1.105615,  1.426492), nrow = 2)
  true_Phi1[1, 1] <- 1 # for identifiability issues

  # true_Phi2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = true_r,
  #                                                   rangeVar = c(0.7, 0.7))$Sigma
  true_Phi2 <- matrix(c(0.7000000, 0.1727312, 0.1727312, 0.7000000), nrow = 2)
  true_Phi2[1, 1] <- 1 # for identifiability issues
  true_Phi_all <- rbind(true_Phi1, true_Phi2)

  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  # true_Omega1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = true_p,
  #                                                     rangeVar = c(1, 1.7))$Sigma
  true_Omega1 <- matrix(c(1.4855139,  0.9049113, -0.9063446,
                          0.9049113,  1.3814824, -1.2298301,
                          -0.9063446, -1.2298301,  1.1979135), nrow = 3)

  # true_Omega2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = true_p,
  #                                                     rangeVar = c(0.7, 0.7))$Sigma
  true_Omega2 <- matrix(c(0.7000000, 0.5379098, 0.4837924,
                          0.5379098, 0.7000000, 0.4089374,
                          0.4837924, 0.4089374, 0.7000000), nrow = 3)
  true_Omega_all <- rbind(true_Omega1, true_Omega2)

  simulatedMVdata <- mvplnDataGenerator(nOccasions = true_r,
                                        nResponses = true_p,
                                        nUnits = true_n,
                                        mixingProportions = c(0.79, 0.21),
                                        matrixMean = true_M_all,
                                        phi = true_Phi_all,
                                        omega = true_Omega_all)

  # dataset provided as character
  testthat::expect_error(mvplnClustering(
    dataset = "simulatedMVdata$dataset",
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # dataset provided as logical
  testthat::expect_error(mvplnClustering(
    dataset = NA,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect size for true membership
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership[-1],
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect class for true membership
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input type for gmin
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = "1",
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input for gmin and gmax
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 5,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))


  # Incorrect input type for nChains
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = "3",
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input type for nIterations
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = "300",
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input type for initMethod
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "other",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))


  # Incorrect input type for initMethod
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = NA,
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input type for nInitIterations
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = "1",
    normalize = "Yes",
    numNodes = 2))

  # Incorrect input type for normalize
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "other",
    numNodes = 2))

  # Incorrect input type for numNodes
  testthat::expect_error(mvplnClustering(
    dataset = simulatedMVdata$dataset,
    membership = simulatedMVdata$truemembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 300,
    initMethod = "kmeans",
    nInitIterations = 1,
    normalize = "Yes",
    numNodes = "2"))


})

