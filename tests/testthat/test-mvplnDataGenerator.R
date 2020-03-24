context("Checking for data simulation")
library(mixMVPLN)

test_that("Data generation is as expected", {

  # Generating simulated data
  set.seed(1234)
  true_G <- 2 # number of total G
  true_r <- 2 # number of total occasions
  true_p <- 3 # number of total responses
  true_n <- 1000 # number of total units

  # M is a r x p matrix
  true_M1 <- matrix(rep(6, (true_r*true_p)), ncol = true_p, nrow = true_r, byrow = TRUE)
  true_M2 <- matrix(rep(1, (true_r*true_p)), ncol = true_p, nrow = true_r, byrow = TRUE)
  true_M_all <- rbind(true_M1, true_M2)
  # Phi is a r x r matrix
  library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  true_Phi1 <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                    dim=true_r,
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

  expect_that(length(simulatedMVdata), equals(9))
  expect_that(simulatedMVdata, is_a("mvplnDataGenerator"))
  expect_that(trunc(simulatedMVdata$units), equals(true_n))
  expect_that(trunc(simulatedMVdata$occassions), equals(true_r))
  expect_that(trunc(simulatedMVdata$variables), equals(true_p))
})


context("Checking for invalid user input")
test_that("Data generate error upon invalid user input", {

  # nOccasions provided as character
  expect_error(simulatedMVdata <- mvplnDataGenerator(nOccasions = "true_r",
                                                     nResponses = true_p,
                                                     nUnits = true_n,
                                                     mixingProportions = c(0.79, 0.21),
                                                     matrixMean = true_M_all,
                                                     phi = true_Phi_all,
                                                     omega = true_Omega_all))

  # nOccasions provided as character
  expect_error(simulatedMVdata <- mvplnDataGenerator(nOccasions = true_r,
                                                     nResponses = c(true_p, true_p),
                                                     nUnits = true_n,
                                                     mixingProportions = c(0.79, 0.21),
                                                     matrixMean = true_M_all,
                                                     phi = true_Phi_all,
                                                     omega = true_Omega_all))


  # Generating simulated data - mu has incorrect dimension
  trueMu1 <- c(6.5, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2)

  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))


  # Generating simulated data - sigma has incorrect dimension
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
  trueSigma1 <- diag(5) * 2
  trueSigma2 <- diag(5)

  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))


  # mixingProportions does not sum to 1
  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.2),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))


  # Incorrect ImageName format
  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "Yes",
    ImageName = 1234))

})
