context("Checking Variational Gaussian Approximations Approach")
library(mixMVPLN)

test_that("Checking clustering results", {

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
  clusteringResultsVGA <- mixMVPLN::mvplnVGAclus(
    dataset = sampleData2$dataset,
    membership = sampleData2$truemembership,
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes")



  # Setting numNodes = 2 based on the following entry, otherwise error.
  # "NB: you can’t use unexported functions and you shouldn’t open new graphics
  # devices or use more than two cores. Individual examples shouldn’t
  # take more than 5s."
  # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa

  expect_type(clusteringResultsVGA, "list")
  expect_length(clusteringResultsVGA, 17)
  expect_s3_class(clusteringResultsVGA, "mvplnVGA")
  expect_identical(clusteringResultsVGA$nUnits, 1000L)
  expect_identical(clusteringResultsVGA$nVariables, 3L)
  expect_identical(clusteringResultsVGA$nOccassions, 2L)
  expect_named(clusteringResultsVGA, c("dataset", "nUnits",
                                  "nVariables", "nOccassions",
                                  "normFactors", "gmin", "gmax",
                                  "initalizationMethod", "allResults",
                                  "loglikelihood", "nParameters",
                                  "trueLabels", "ICLAll",
                                  "BICAll", "AICAll",
                                  "AIC3All", "totalTime"))
  expect_output(str(clusteringResultsVGA), "List of 17")
  expect_vector(clusteringResultsVGA$trueLabels, ptype = double(), size = 1000)
  expect_identical(clusteringResultsVGA$BICAll$BICmodelselected, 1)
})

context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
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


  # dataset provided as character
  testthat::expect_error(mvplnVGAclus(
      dataset = "simulatedMVdata$dataset",
      membership = sampleData2$truemembership,
      gmin = 1,
      gmax = 2,
      initMethod = "clara",
      nInitIterations = 2,
      normalize = "Yes"))

  # dataset provided as logical
  testthat::expect_error(mvplnVGAclus(
    dataset = NA,
    membership = sampleData2$truemembership,
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes"))

  # Incorrect size for true membership
  testthat::expect_error(mvplnVGAclus(
    dataset = "simulatedMVdata$dataset",
    membership = sampleData2$truemembership[-1],
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes"))

  # Incorrect class for true membership
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes"))


  # Incorrect input type for gmin
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = "1",
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes"))


  # Incorrect input for gmin and gmax
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 5,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = 2,
    normalize = "Yes"))

  # Incorrect input type for nChains
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = "2",
    normalize = "Yes"))

  # Incorrect input type for nInitIterations
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    initMethod = "clara",
    nInitIterations = "2",
    normalize = "Yes"))


  # Incorrect input type for initMethod
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    initMethod = "other",
    nInitIterations = 2,
    normalize = "Yes"))

  # Incorrect input type for normalize
  testthat::expect_error(mvplnVGAclus(
    dataset = simulatedMVdata$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    initMethod = "other",
    nInitIterations = 2,
    normalize = "other"))
})
# [END]
