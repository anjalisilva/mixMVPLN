context("Checking for parallel clustering performance")
library(mixMVPLN)

test_that("Checking clustering results", {

  set.seed(1234)
  true_G <- 2 # number of total G
  true_r <- 2 # number of total occasions
  true_p <- 3 # number of total responses
  true_n <- 20 # number of total units

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
  true_Phi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                    dim = true_r,
                                                    rangeVar = c(1, 1.7))$Sigma
  true_Phi1[1, 1] <- 1 # for identifiability issues

  true_Phi2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                    dim = true_r,
                                                    rangeVar = c(0.7, 0.7))$Sigma
  true_Phi2[1, 1] <- 1 # for identifiability issues
  true_Phi_all <- rbind(true_Phi1, true_Phi2)

  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  true_Omega1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                      dim = true_p,
                                                      rangeVar = c(1, 1.7))$Sigma
  true_Omega2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                      dim = true_p,
                                                      rangeVar = c(0.7, 0.7))$Sigma
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
                                       gmax = 2,
                                       nChains = 3,
                                       nIterations = 300,
                                       initMethod = "kmeans",
                                       nInitIterations = 1,
                                       normalize = "Yes",
                                       numNodes = 2)

  # Setting numNodes = 2 based on the following entry, otherwise error.
  # "NB: you can’t use unexported functions and you shouldn’t open new graphics
  # devices or use more than two cores. Individual examples shouldn’t
  # take more than 5s."
  # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa

  expect_that(length(clusteringResults), equals(16))
  expect_that(clusteringResults, is_a("mplnParallel"))
  expect_that(clusteringResults$initalization_method, equals("kmeans"))
  numPara <- c(27)
  expect_that(clusteringResults$numb_of_parameters, equals(numPara))
  expect_that(clusteringResults$true_labels, equals(simulatedCounts$trueMembership))
  expect_that(trunc(clusteringResults$ICL_all$ICLmodelselected), equals(2))
  expect_that(trunc(clusteringResults$AIC_all$AICmodelselected), equals(2))
  expect_that(trunc(clusteringResults$BIC_all$BICmodelselected), equals(2))
})

context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- mplnDataGenerator(nObservations = 500,
                                        dimensionality = 6,
                                        mixingProportions = c(0.79, 0.21),
                                        mu = rbind(trueMu1, trueMu2),
                                        sigma = rbind(trueSigma1, trueSigma2),
                                        produceImage = "No")

  # dataset provided as character
  expect_error(mplnParallel(dataset = "dataset",
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnParallel(dataset = NA,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership[-1],
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = "trueMembership",
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = "1",
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 5,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = "3",
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = "1000",
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "other",
                            nInitIterations = 3,
                            normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = NA,
                            nInitIterations = 3,
                            normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = "3",
                            normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnParallel(dataset = simulatedCounts$dataset,
                            membership = simulatedCounts$trueMembership,
                            gmin = 1,
                            gmax = 2,
                            nChains = 3,
                            nIterations = 1000,
                            initMethod = "kmeans",
                            nInitIterations = 3,
                            normalize = "Other"))

})




