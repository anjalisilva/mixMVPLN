# **Should have all the function R files and MPLN.stan in the same working directory as this file

# Read all the necessary functions and check that packages needed are present
source("Setup.R")

#####################################  DATA GENERATION  #####################################
# Generating simulated data

set.seed(1)

true_G <- 2 # number of total G
true_r <- 2 # number of total occasions
true_p <- 3 # number of total responses
true_n <- 100 # number of total units

# loading needed packages for generating data
if (!require(clusterGeneration)) install.packages("clusterGeneration") 
if (!require(mvtnorm)) install.packages("mvtnorm") 
if (!require(edgeR)) install.packages("edgeR") 
if (!require(mclust)) install.packages("mclust") 

# Mu is a r x p matrix
true_M1 <- matrix(rep(6,(true_r*true_p)), ncol=true_p, nrow=true_r, byrow=TRUE)
true_M2 <- matrix(rep(1,(true_r*true_p)), ncol=true_p, nrow=true_r, byrow=TRUE)
true_M_all <- rbind(true_M1, true_M2)

# Phi is a r x r matrix
# Covariance matrix containing variances and covariances between r occasions
true_Phi1 <- genPositiveDefMat("unifcorrmat",dim=true_r, rangeVar=c(1,1.7))$Sigma
true_Phi1[1,1] <- 1 # for identifiability issues
true_Phi2 = genPositiveDefMat("unifcorrmat",dim=true_r, rangeVar=c(0.7,0.7))$Sigma
true_Phi2[1,1] = 1 # for identifiability issues
true_Phi_all = rbind(true_Phi1,true_Phi2)

# Omega is a p x p matrix 
# Covariance matrix containing variances and covariances between p responses
true_Omega1 <- genPositiveDefMat("unifcorrmat",dim=true_p, rangeVar=c(1,1.7))$Sigma
true_Omega2 <- genPositiveDefMat("unifcorrmat",dim=true_p, rangeVar=c(0.7,0.7))$Sigma
true_Omega_all <- rbind(true_Omega1,true_Omega2)

simulated_counts <- Datagenerator(i=1, r=true_r, p=true_p, n=true_n, pi_g=c(0.79,0.21), mu=true_M_all, phi=true_Phi_all, omega=true_Omega_all)

#####################################################################################################
# Clustering data for G = 1:3
testing_dataset <- simulated_counts # Assign test dataset using the variable name 'testing_dataset'
clustering_results <- MVPLNClustering(dataset=testing_dataset$dataset, membership=simulated_counts$truemembership, Gmin=1, Gmax=3, n_chains=3, n_iterations=300, init_method="kmeans", n_init_iterations=5)

