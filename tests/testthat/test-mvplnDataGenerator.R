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

  testthat::expect_that(length(simulatedMVdata), equals(9))
  testthat::expect_that(simulatedMVdata, is_a("mvplnDataGenerator"))
  testthat::expect_that(trunc(simulatedMVdata$units), equals(true_n))
  testthat::expect_that(trunc(simulatedMVdata$occassions), equals(true_r))
  testthat::expect_that(trunc(simulatedMVdata$variables), equals(true_p))
})


context("Checking for invalid user input")
test_that("Data generate error upon invalid user input", {

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


  # nOccasions provided as character
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = "true_r",
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = true_Omega_all))

  # nResponses provided as vector
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = c(true_p, true_p),
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = true_Omega_all))

  # nUnits provided as a decimal
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = 1000.25,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = true_Omega_all))


  # mixingProportions doesn't sum to 1
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.8, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = true_Omega_all))


  # Generating simulated data - M has incorrect dimension
  M1_testing <- matrix(rep(6, (true_r*(true_p + 1))),
                       ncol = (true_p + 1),
                       nrow = true_r, byrow = TRUE)
  M2_testing <- matrix(rep(1, (true_r*(true_p + 1))),
                       ncol = (true_p + 1),
                       nrow = true_r, byrow = TRUE)
  M_testing <- rbind(M1_testing, M2_testing)


  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = M_testing,
    phi = true_Phi_all,
    omega = true_Omega_all))


  # Generating simulated data - Phi has incorrect dimension
  # Phi is a r x r matrix
  library(clusterGeneration)
  # Covariance matrix containing variances and covariances between r occasions
  Phi1_testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                       dim = (true_r + 1),
                                                       rangeVar = c(1, 1.7))$Sigma
  Phi1_testing[1, 1] <- 1 # for identifiability issues
  Phi2_testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                       dim = (true_r + 1),
                                                       rangeVar = c(0.7, 0.7))$Sigma
  Phi2_testing[1, 1] <- 1 # for identifiability issues
  Phi_testing <- rbind(Phi1_testing, Phi2_testing)

  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = Phi_testing,
    omega = true_Omega_all))

  # Generating simulated data - omega has incorrect dimension
  # Omega is a p x p matrix
  # Covariance matrix containing variance and covariances of p responses/variables
  Omega1_testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                         dim = (true_p + 1),
                                                         rangeVar = c(1, 1.7))$Sigma
  Omega2_testing <- clusterGeneration::genPositiveDefMat(covMethod = "unifcorrmat",
                                                         dim = (true_p + 1),
                                                         rangeVar = c(0.7, 0.7))$Sigma
  Omega_testing <- rbind(Omega1_testing, Omega2_testing)

  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = Omega_testing))


  # Generating simulated data - omega is not a matrix, but a data frame
  # Omega is a p x p matrix
  testthat::expect_error(simulatedMVdata <- mvplnDataGenerator(
    nOccasions = true_r,
    nResponses = true_p,
    nUnits = true_n,
    mixingProportions = c(0.79, 0.21),
    matrixMean = true_M_all,
    phi = true_Phi_all,
    omega = data.frame(true_Omega_all)))

})
