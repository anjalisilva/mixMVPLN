# clustering function
cluster_mpln <- function(r,p,z,dataset,G,mod,norm_factors,n_chains,n_iterations,initialization){
  #print("Entering cluster_mpln")
  n=nrow(dataset)
  
  Phi_all_outer = Omega_all_outer = M_all_outer = Sigma_all_outer = list()
  median_mu_outer = median_sigma_outer = median_phi_outer = median_omega_outer = list()
  conv_outer=0
  it_outer=2
  obs = PI = logL = norm_mu_outer = norm_sigma_outer = vector()
  
  if (all(is.na(initialization))==TRUE){
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
  }else{ # running after initialization has been done
    M_all_outer[[1]] = M = initialization$finalmu
    Phi_all_outer[[1]] = Phi = initialization$finalphi
    Omega_all_outer[[1]] = Omega = initialization$finalomega
    Sigma_all_outer[[1]] = Sigma = initialization$finalsigma
    z = initialization$probaPost
    n_iterations = initialization$FinalRstan_iterations
  }
  
  
  # start the loop
  while(!conv_outer){
    cat("************** Running for G =",G,"and Iteration =",it_outer,"******************")
    
    obs=apply(z, 2, sum) # number of observations in each group
    PI=sapply(obs, function(x) x/n)  # obtain probability of each group
    
    stanresults = stanrun(r=r,p=p,dataset=dataset,G=G,n_chains=n_chains,n_iterations=n_iterations,mod=mod,Mu=M_all_outer[[it_outer-1]],Sigma=Sigma_all_outer[[it_outer-1]],norm_factors=norm_factors)
    
    # update parameters
    paras = parameter_estimation(r=r, p=p, G=G, z=z, fit=stanresults$fitrstan, n_iterations=stanresults$n_iterations, n_chains=n_chains)
    
    
    M_all_outer[[it_outer]] = paras$Mu
    Sigma_all_outer[[it_outer]] = paras$Sigma
    Phi_all_outer[[it_outer]] = paras$Phi
    Omega_all_outer[[it_outer]] = paras$Omega
    
    theta_Stan = paras$theta
    
    vectorized_M=t(sapply(c(1:G), function(g) ( M_all_outer[[it_outer]][((g-1)*r+1):(g*r),]) ))
    
    logL[it_outer] = calc_loglikelihood(dataset=dataset,z=z,G=G,PI=PI, norm_factors=norm_factors,M=vectorized_M,Sigma=Sigma_all_outer[[it_outer]], theta_Stan=theta_Stan)
    #plot(logL[-1], xlab="iteration", ylab=paste("Initialization logL value for",G) )
    
    
    threshold_outer <- 15
    if(it_outer>(threshold_outer+1)){
      
      if (all(heidel.diag(logL[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1) || it_outer>100){
        programclust <- vector()
        programclust <- map(z)
        
        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[!K])>0){ # J[!K] tells which are empty clusters
          z <- z[,-J[!K]]
          programclust <- map(z)
        }
        
        conv_outer <- 1
      } 
    }
    
    
    if(conv_outer!=1){ # only update until convergence, not after
      # updating z value
      z= zvalue_calculation(PI=PI,dataset=dataset,M=vectorized_M,G=G,Sigma=Sigma_all_outer[[it_outer]],theta_Stan=theta_Stan,norm_factors=norm_factors)
      it_outer=it_outer+1
      n_iterations = n_iterations+10
    }
  } # end of loop
  
  print("Ending cluster_mpln")
  
  results <- list(finalmu=M_all_outer[[it_outer]]+ do.call("rbind", rep(list(matrix(norm_factors,byrow=TRUE,ncol=p)),G)), 
    finalsigma=Sigma_all_outer[[it_outer]],
    finalphi = Phi_all_outer[[it_outer]],
    finalomega = Omega_all_outer[[it_outer]],
    allmu = lapply(M_all_outer, function(x) (x+ do.call("rbind", rep(list(matrix(norm_factors,byrow=TRUE,ncol=p)),G)))),      
    allsigma = Sigma_all_outer, 
    allphi = Phi_all_outer,
    allomega = Omega_all_outer,
    clusterlabels =programclust,
    iterations = it_outer, 
    FinalRstan_iterations = n_iterations,
    proportion = PI, 
    loglikelihood = logL,
    probaPost = z)
  
  class(results) <- "MPLNcluster"
  return(results)
  
}
