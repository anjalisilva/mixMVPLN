# Calling initialization function
initialization_run <- function(r,p,z,dataset,G,mod,norm_factors,n_chains,n_iterations){
  n=nrow(z)
  
  Phi_all_outer = Omega_all_outer = M_all_outer = Sigma_all_outer = list()
  median_mu_outer = median_sigma_outer = median_phi_outer = median_omega_outer = list()
  conv_outer=0
  it_outer=2
  obs=PI = logL = norm_mu_outer = norm_sigma_outer = vector()
  
  
  # initialize Phi; rxr times G
  Phi_all_outer[[1]] = Phi = do.call("rbind", rep(list(diag(r)), G)) 
  
  # initialize Omega; pxp times G
  Omega_all_outer[[1]] = Omega = do.call("rbind", rep(list(diag(p)), G)) 
  
  M = matrix(NA,ncol=p, nrow=r*G) #initialize M = rxp matrix
  Sigma = do.call("rbind", rep(list(diag(r*p)), G)) # sigma (rp by rp)
  for (g in 1:G){
    M[((g-1)*r+1):(g*r),] = matrix(log(mean(dataset[c(which(z[,g]==1)),])),ncol=p, nrow=r)
    Sigma[((g-1)*(r*p)+1):(g*(r*p)),] = (cov(log(dataset[c(which(z[,g]==1)),]+(1/3)))) #diag(r*p) 
  }
  M_all_outer[[1]] = M
  Sigma_all_outer[[1]] = Sigma
  

  for(g in 1:G){
      obs[g] = sum(z[,g]) # number of observations in each group
      PI[g] = obs[g]/n  # obtain probability of each group
  }
    
  stanresults = stanrun(r=r,p=p,dataset=dataset,G=G,n_chains=n_chains,n_iterations=n_iterations,mod=mod,Mu=M_all_outer[[it_outer-1]],Sigma=Sigma_all_outer[[it_outer-1]],norm_factors=norm_factors)
    
  n_iterations = stanresults$n_iterations
    
  # update parameters
  paras = parameter_estimation(r=r, p=p, G=G, z=z, fit=stanresults$fitrstan, n_iterations=n_iterations, n_chains=n_chains)
    
  M_all_outer[[it_outer]] = paras$Mu
  Sigma_all_outer[[it_outer]] = paras$Sigma
  Phi_all_outer[[it_outer]] = paras$Phi
  Omega_all_outer[[it_outer]] = paras$Omega
    
  theta_Stan = paras$theta
    
    vectorized_M=matrix(NA,ncol=r*p,nrow=G)
    for (g in 1:G){
      vectorized_M[g,] = as.vector(t( M_all_outer[[it_outer]][((g-1)*r+1):(g*r),]))
    }
    
    # calculate logL
    logL[it_outer] = calc_loglikelihood(dataset=dataset,z=z,G=G,PI=PI,norm_factors=norm_factors,M=vectorized_M,Sigma=Sigma_all_outer[[it_outer]], theta_Stan=theta_Stan)
    #plot(logL[-1], xlab="iteration", ylab=paste("Initialization logL value for",G) )
    
    z <- zvalue_calculation(PI=PI,dataset=dataset,M=vectorized_M,G=G,Sigma=Sigma_all_outer[[it_outer]],theta_Stan=theta_Stan, norm_factors=norm_factors)
    

  results <- list(finalmu=M_all_outer[[it_outer]]+ do.call("rbind", rep(list(matrix(norm_factors,byrow=TRUE,ncol=p)),G)), 
    finalsigma=Sigma_all_outer[[it_outer]],
    finalphi = Phi_all_outer[[it_outer]],
    finalomega = Omega_all_outer[[it_outer]],
    probaPost =z,
    FinalRstan_iterations = n_iterations,
    loglikelihood = logL)
  
  class(results) <- "MPLN_initialization"
  return(results)
  # Developed by Anjali Silva
}