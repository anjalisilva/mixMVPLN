context("Checking for data simulation")
library(mixMVPLN)

test_that("Data generation is as expected", {

  # Generating simulated data
  set.seed(1234)
  trueG <- 2 # number of total G
  truer <- 2 # number of total occasions
  truep <- 3 # number of total responses
  trueN <- 1000 # number of total units

  # M is a r x p matrix
  trueM1 <- matrix(rep(6, (truer * truep)),
                    ncol = truep, nrow = truer,
                    byrow = TRUE)
  trueM2 <- matrix(rep(1, (truer * truep)),
                    ncol = truep, nrow = truer,
                    byrow = TRUE)
  trueMall <- rbind(trueM1, trueM2)

  # Phi is a r x r matrix
  # library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  # truePhi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = truer,
  #                                                   rangeVar = c(1, 1.7))$Sigma
  truePhi1 <- matrix(c(1.435610, -1.105615,
                       -1.105615,  1.426492), nrow = 2)
  truePhi1[1, 1] <- 1 # for identifiability issues

  # truePhi2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = truer,
  #                                                   rangeVar = c(0.7, 0.7))$Sigma
  truePhi2 <- matrix(c(0.7000000, 0.1727312,
                       0.1727312, 0.7000000), nrow = 2)
  truePhi2[1, 1] <- 1 # for identifiability issues
  truePhiAll <- rbind(truePhi1, truePhi2)

  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  # trueOmega1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = truep,
  #                                                    rangeVar = c(1, 1.7))$Sigma
  trueOmega1 <- matrix(c(1.4855139,  0.9049113, -0.9063446,
                         0.9049113,  1.3814824, -1.2298301,
                        -0.9063446, -1.2298301,  1.1979135), nrow = 3)
  # trueOmega2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = truep,
  #                                                     rangeVar = c(0.7, 0.7))$Sigma
  trueOmega2 <- matrix(c(0.7000000, 0.5379098, 0.4837924,
                         0.5379098, 0.7000000, 0.4089374,
                         0.4837924, 0.4089374, 0.7000000), nrow = 3)
  trueOmegaAll <- rbind(trueOmega1, trueOmega2)

  simulatedMVdata <- mvplnDataGenerator(nOccasions = truer,
                                        nResponses = truep,
                                        nUnits = trueN,
                                        mixingProportions = c(0.79, 0.21),
                                        matrixMean = trueMall,
                                        phi = truePhiAll,
                                        omega = trueOmegaAll)

  testthat::expect_that(length(simulatedMVdata), equals(10))
  testthat::expect_that(simulatedMVdata, is_a("mvplnDataGenerator"))
  testthat::expect_that(trunc(simulatedMVdata$units), equals(trueN))
  testthat::expect_that(trunc(simulatedMVdata$occassions), equals(truer))
  testthat::expect_that(trunc(simulatedMVdata$variables), equals(truep))

  expect_type(simulatedMVdata, "list")
  expect_length(simulatedMVdata, 10)
  expect_s3_class(simulatedMVdata, "mvplnDataGenerator")
  expect_identical(simulatedMVdata$units, 1000)
  expect_identical(simulatedMVdata$variables, 3)
  expect_identical(simulatedMVdata$occassions, 2)
  expect_named(simulatedMVdata, c("dataset",
                                  "truemembership",
                                  "units",
                                  "occassions",
                                  "variables",
                                  "mixingProportions",
                                  "means",
                                  "phi",
                                  "omega",
                                  "data2D"))
  expect_output(str(simulatedMVdata), "List of 10")
  expect_vector(simulatedMVdata$truemembership, ptype = double(), size = 1000)
  expect_vector(simulatedMVdata$mixingProportions, c(0.79, 0.21))
  expect_identical(simulatedMVdata$mixingProportions, c(0.79, 0.21))
})


context("Checking for invalid user input")
test_that("Data generate error upon invalid user input", {

  # Generating simulated data
  set.seed(1234)
  trueG <- 2 # number of total G
  truer <- 2 # number of total occasions
  truep <- 3 # number of total responses
  trueN <- 1000 # number of total units

  # M is a r x p matrix
  trueM1 <- matrix(rep(6, (truer*truep)),
                   ncol = truep, nrow = truer,
                   byrow = TRUE)
  trueM2 <- matrix(rep(1, (truer*truep)),
                   ncol = truep, nrow = truer,
                   byrow = TRUE)
  trueMall <- rbind(trueM1, trueM2)

  # Phi is a r x r matrix
  # library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  # truePhi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = truer,
  #                                                   rangeVar = c(1, 1.7))$Sigma
  truePhi1 <- matrix(c(1.435610, -1.105615, -1.105615,  1.426492), nrow = 2)
  truePhi1[1, 1] <- 1 # for identifiability issues

  # truePhi2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                   dim = truer,
  #                                                   rangeVar = c(0.7, 0.7))$Sigma
  truePhi2 <- matrix(c(0.7000000, 0.1727312, 0.1727312, 0.7000000), nrow = 2)
  truePhi2[1, 1] <- 1 # for identifiability issues
  truePhiAll <- rbind(truePhi1, truePhi2)

  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  # trueOmega1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = truep,
  #                                                    rangeVar = c(1, 1.7))$Sigma
  trueOmega1 <- matrix(c(1.4855139,  0.9049113, -0.9063446,
                          0.9049113,  1.3814824, -1.2298301,
                          -0.9063446, -1.2298301,  1.1979135), nrow = 3)

  # trueOmega2 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                     dim = truep,
  #                                                     rangeVar = c(0.7, 0.7))$Sigma
  trueOmega2 <- matrix(c(0.7000000, 0.5379098, 0.4837924,
                          0.5379098, 0.7000000, 0.4089374,
                          0.4837924, 0.4089374, 0.7000000), nrow = 3)

  trueOmegaAll <- rbind(trueOmega1, trueOmega2)


  # nOccasions provided as character
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = "truer",
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = trueOmegaAll))

  # nResponses provided as vector
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = c(truep, truep),
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = trueOmegaAll))

  # nUnits provided as a decimal
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = 1000.25,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = trueOmegaAll))


  # mixingProportions doesn't sum to 1
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.8, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = trueOmegaAll))


  # Generating simulated data - M has incorrect dimension
  M1Testing <- matrix(rep(6, (truer * (truep + 1))),
                       ncol = (truep + 1),
                       nrow = truer, byrow = TRUE)
  M2Testing <- matrix(rep(1, (truer * (truep + 1))),
                       ncol = (truep + 1),
                       nrow = truer, byrow = TRUE)
  MTesting <- rbind(M1Testing, M2Testing)


  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = MTesting,
    phi = truePhiAll,
    omega = trueOmegaAll))


  # Generating simulated data - Phi has incorrect dimension
  # Phi is a r x r matrix
  # library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  # Phi1Testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                      dim = (truer + 1),
  #                                                      rangeVar = c(1, 1.7))$Sigma
  Phi1Testing <- matrix(c(1.3599758, 0.2883504, 0.7627700,
                           0.2883504, 1.4855139, 0.3311862,
                           0.7627700, 0.3311862, 1.3814824), nrow = 3)
  Phi1Testing[1, 1] <- 1 # for identifiability issues

  # Phi2Testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                      dim = (truer + 1),
  #                                                      rangeVar = c(0.7, 0.7))$Sigma
  Phi2Testing <- matrix(c(0.700000, -0.3762910,  0.4978850,
                           -0.376291,  0.7000000, -0.2129134,
                           0.497885, -0.2129134,  0.7000000), nrow = 3)
  Phi2Testing[1, 1] <- 1 # for identifiability issues
  PhiTesting <- rbind(Phi1Testing, Phi2Testing)

  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = PhiTesting,
    omega = trueOmegaAll))

  # Generating simulated data - omega has incorrect dimension
  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  # Omega1Testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                        dim = (truep + 1),
  #                                                        rangeVar = c(1, 1.7))$Sigma
  Omega1Testing <- matrix(c(1.1408736, -0.48972027, -0.11173329,  0.7281382,
                            -0.4897203,  1.18116687,  0.05143694, -0.3539734,
                            -0.1117333,  0.05143694,  1.69450529,  0.8320678,
                             0.7281382, -0.35397344,  0.83206778,  1.5651466), nrow = 4)

  # Omega2Testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
  #                                                        dim = (truep + 1),
  #                                                        rangeVar = c(0.7, 0.7))$Sigma
  Omega2Testing <- matrix(c(0.70000000,  0.05289814,  0.1825991,  0.4078816,
                             0.05289814,  0.70000000, -0.1909522, -0.2348412,
                             0.18259910, -0.19095216,  0.7000000, -0.1719357,
                             0.40788158, -0.23484119, -0.1719357,  0.7000000), nrow = 4)
  Omega_testing <- rbind(Omega1Testing, Omega2Testing)

  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = Omega_testing))


  # Generating simulated data - omega is not a matrix, but a data frame
  # Omega is a p x p matrix
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = truer,
    nResponses = truep,
    nUnits = trueN,
    mixingProportions = c(0.79, 0.21),
    matrixMean = trueMall,
    phi = truePhiAll,
    omega = data.frame(trueOmegaAll)))
})
# [END]
