
#### functions ####

# initialization function
initialization_function <- function(r,p,dataset,gmodel,mod,norm_factors,n_chains,n_iterations,init_method,n_init_iterations){
  z <- init_runs <- list()
  logL_init <- vector()
  n <- nrow(dataset)
  d <- ncol(dataset)
  for(iterations in 1:n_init_iterations){
    if (init_method=="kmeans" | is.na(init_method)){
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) 
      suppressWarnings(library(mclust))
      z[[iterations]]<-unmap(kmeans(log(dataset+1/3),gmodel)$cluster)
    }else if (init_method=="random"){
      if(gmodel==1){ # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times=n), ncol=gmodel, nrow=n)
      } else { # generating z if g>1
        z_conv=0
        while(!z_conv){ # ensure that dimension of z is same as G (i.e. 
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(rmultinom(n, size = 1, prob=rep(1/gmodel,gmodel))) 
          if(length(which(colSums(z[[iterations]])>0)) ==gmodel){
            z_conv=1
          }
        }
      }
    }else if (init_method=="medoids"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) # loading needed packages
      suppressWarnings(library(mclust))
      
      z[[iterations]]<-unmap(pam(log(dataset+1/3),k=gmodel)$cluster)
    }else if (init_method=="clara"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(clara(log(dataset+1/3),k=gmodel)$cluster)
    }else if (init_method=="fanny"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(fanny(log(dataset+1/3),k=gmodel)$cluster)
    }
    init_runs[[iterations]]=initialization_run(r=r, p=p, z=z[[iterations]], mod=mod, dataset=dataset, G=gmodel, norm_factors=norm_factors, n_chains=n_chains, n_iterations=n_iterations)
    logL_init[iterations] <- unlist(tail((init_runs[[iterations]]$loglikelihood), n=1)) 
  }
  
  initialization <- init_runs[[which(logL_init==max(logL_init, na.rm = TRUE))[1]]]

  return(initialization)
}

# calling initialization function
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
}

# z value calculation function
zvalue_calculation <- function(PI,dataset,M,G,Sigma,theta_Stan,norm_factors){
  # print("Entering zvalue_calculation")
  d <- ncol(dataset)
  n <- nrow(dataset)
  forz=sapply(c(1:G), function(g) sapply(c(1:n), function(i) PI[g]*exp(t(dataset[i,])%*%(theta_Stan[[g]][i,]+norm_factors)-sum(exp(theta_Stan[[g]][i,]+norm_factors))-sum(lfactorial(dataset[i,]))-
      d/2*log(2*pi)-1/2*log(det(Sigma[((g-1)*d+1):(g*d),]))-0.5*t(theta_Stan[[g]][i,]-M[g,])%*%solve(Sigma[((g-1)*d+1):(g*d),])%*%(theta_Stan[[g]][i,]-M[g,])) ) )
  
  if (G==1){
    errorpossible <- Reduce(intersect, list(which(forz==0),which(rowSums(forz)==0)))
    zvalue <- forz/rowSums(forz)
    zvalue[errorpossible,]<-1
  }else {zvalue <- forz/rowSums(forz)}
  
  return(zvalue)
}

# likelihood calculation function
calc_loglikelihood <- function(dataset,z,G,PI,norm_factors,M,Sigma,theta_Stan){
  n <- nrow(dataset)
  like <- matrix(NA, nrow=n, ncol=G)
  d <- ncol(dataset)
  like= sapply(c(1:G), function(g) sapply(c(1:n), function(i) z[i,g] *(log(PI[g]) +
      t(dataset[i,])%*%(theta_Stan[[g]][i,]+norm_factors)-sum(exp(theta_Stan[[g]][i,]+norm_factors))-sum(lfactorial(dataset[i,]))-
      d/2*log(2*pi)-1/2*log(det(Sigma[((g-1)*d+1):(g*d),]))-0.5*t(theta_Stan[[g]][i,]-M[g,])%*%solve(Sigma[((g-1)*d+1):(g*d),])%*%(theta_Stan[[g]][i,]-M[g,])) ) )
  loglike <- sum(rowSums(like))
  return(loglike)
}

# parameter estimation function
parameter_estimation <- function(r, p, G, z, fit, n_iterations, n_chains){
  d = r*p
  n = nrow(z)
  
  M = matrix(NA,ncol=p, nrow=r*G) #initialize M = rxp matrix
  Sigma = do.call("rbind", rep(list(diag(r*p)), G)) # sigma 
  Phi = do.call("rbind", rep(list(diag(r)), G)) 
  Omega = do.call("rbind", rep(list(diag(p)), G)) 
  
  # generating stan values
  E_theta_omega = E_theta_phi = E_theta_phi2 =list()
  
  # update parameters
  theta_mat=lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)), function(o) as.matrix(fit[[g]])[, c(o,sapply(c(1:n), function(i) c(1:(d-1))*n+i)[,o]) ]) )
  
  theta_Stan=lapply(as.list(c(1:G)), function(g) t(sapply(c(1:n), function(o) colMeans(as.matrix(fit[[g]])[,c(o, sapply(c(1:n), function(i) c(1:(d-1))*n+i)[,o]) ]))) )
  
  
  # updating mu
  M=lapply(as.list(c(1:G)), function(g) matrix(data = colSums(z[,g]*theta_Stan[[g]])/sum(z[,g]), ncol=p, nrow=r, byrow = T) )
  M=do.call(rbind,M)
  
  # updating phi
  E_theta_phi=lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)), function(i) lapply(as.list(c(1:nrow(theta_mat[[g]][[i]]))), function(e) ((matrix(theta_mat[[g]][[i]][e,],r,p, byrow = T))-M[((g-1)*r+1):(g*r),])%*%solve(Omega[((g-1)*p+1):(g*p),])%*%t((matrix(theta_mat[[g]][[i]][e,],r,p, byrow = T))-M[((g-1)*r+1):(g*r),]) ) ) )
  E_theta_phi2=lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)), function(i) z[i,g]*Reduce("+",E_theta_phi[[g]][[i]])/((0.5*n_iterations)*n_chains) ) )
  Phi=lapply(as.list(c(1:G)), function(g) (Reduce("+",E_theta_phi2[[g]])/(p*sum(z[,g]))))
  Phi=do.call(rbind,Phi)
  
  # updating omega
  E_theta_omega=lapply(as.list(c(1:G)), function(g) lapply(as.list(c(1:n)), function(i) lapply(as.list(1:nrow(theta_mat[[g]][[i]])), function(e) t((matrix(theta_mat[[g]][[i]][e,],r,p, byrow = T))-M[((g-1)*r+1):(g*r),])%*%solve(Phi[((g-1)*r+1):(g*r),])%*%((matrix(theta_mat[[g]][[i]][e,],r,p, byrow = T))-M[((g-1)*r+1):(g*r),]) ) ) )
  E_theta_omega2=lapply(as.list(c(1:G)) ,function(g) lapply(as.list(c(1:n)), function(i)  z[i,g]*Reduce("+",E_theta_omega[[g]][[i]])/((0.5*n_iterations)*n_chains) ) )
  Omega=lapply(as.list(c(1:G)), function(g) Reduce("+",E_theta_omega2[[g]])/(r*sum(z[,g]))  )
  Omega=do.call(rbind,Omega)
  
  Sigma=lapply(as.list(c(1:G)), function(g) kronecker(Phi[((g-1)*r+1):(g*r),], Omega[((g-1)*p+1):(g*p),]) )
  Sigma=do.call(rbind,Sigma)
  
  results <- list(Mu = M,
    Sigma = Sigma,
    Omega = Omega,
    Phi = Phi,
    theta = theta_Stan)
  class(results) <- "RStan"
  return(results)
  
}

# stan sampling function
stanrun <- function(r,p,dataset,G,n_chains,n_iterations,mod,Mu,Sigma,norm_factors){
  n=nrow(dataset)
  
  fitrstan=list()
  for (g in 1:G){
    data1 = list(d=ncol(dataset), N=nrow(dataset), y=dataset, mu=as.vector(Mu[((g-1)*r+1):(g*r),]), Sigma=Sigma[((g-1)*(r*p)+1):(g*(r*p)),], normfactors=as.vector(norm_factors))
    stanproceed <- 0
    try=1
    
    while (!stanproceed){
      
      fitrstan[[g]]<-sampling(object=mod,
        data=data1,
        iter=n_iterations, chains = n_chains, verbose=FALSE, refresh=-1)
      if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) == TRUE && all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) == TRUE){
        #if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) == TRUE){
        #  print("Rhat met")}
        #if (all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) == TRUE){
        #  print("N_eff met")}
        stanproceed <- 1
      } else if(all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE || all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) != TRUE){
        if(try == 5){ # stop after 5 tries
          stanproceed = 1
        }
        #if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE){
        #  print("Rhat NOT met")}
        #if (all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) != TRUE){
        #  print("N_eff NOT met")}
        n_iterations = n_iterations+100
        try=try+1
      }
    } # close while loop
  }# close g loop
  
  results <- list(fitrstan = fitrstan,
    n_iterations = n_iterations)
  class(results) <- "RStan"
  return(results)
}

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

# calling the clustering function
calling_clustering <- function(n, r, p, d, dataset, Gmin, Gmax, n_chains, n_iterations=NA, init_method=NA, n_init_iterations=NA, norm_factors, mod){
  ptm_inner = proc.time() 
  
  for (gmodel in 1:(Gmax-Gmin+1)){
    
    if(length(1:(Gmax-Gmin+1)) == Gmax){
      clustersize = gmodel
    }else if(length(1:(Gmax-Gmin+1)) < Gmax){
      clustersize = seq(Gmin, Gmax, 1)[gmodel]
    }
    
    if(n_init_iterations!=0){
      initializeruns=initialization_function(r=r, p=p, gmodel=clustersize, mod=mod, dataset=dataset, init_method=init_method, n_init_iterations=n_init_iterations, n_chains=n_chains, n_iterations=n_iterations, norm_factors=norm_factors)
      allruns=cluster_mpln(z=NA,n_iterations=NA,initialization=initializeruns,r=r,p=p,dataset=dataset,G=clustersize,mod=mod,norm_factors=norm_factors,n_chains=n_chains)
    }else if(n_init_iterations == 0){
      allruns=cluster_mpln(z=unmap(kmeans(log(dataset+1/3),clustersize)$cluster),initialization=NA,r=r,p=p,dataset=dataset,G=clustersize, mod=mod,norm_factors=norm_factors,n_chains=n_chains,n_iterations=n_iterations)
    }
  }
  
  final_inner <- proc.time()-ptm_inner
  RESULTS <- list(  gmin = Gmin,
    gmax = Gmax,
    initalization_method = init_method,
    allresults = allruns,
    totaltime = final_inner)
  
  
  class(RESULTS) <- "MPLN"
  return(RESULTS)
}

# parameter calculation function
calculate_parameters <- function(g,dataset){
  p=ncol(dataset[[1]])
  r=nrow(dataset[[1]])
  mu_para <- r*p*g
  sigma_para<- g/2 *( (r*(r+1)) + (p*(p+1)))
  pi_para <- g-1 # because if you have g-1 parameters, you can do 1-these to get the last one
  paratotal <- mu_para+sigma_para+pi_para # total parameters are
  return(paratotal)
}

# BIC function 
BIC_function <- function(ll, k, n, run, gmin, gmax){
  BIC <- -2*ll+ (k* log(n))
  BICmodel <- seq(gmin, gmax, 1)[grep(min(BIC,na.rm = TRUE), BIC)]
  BICmodel_labels <- run[[grep(min(BIC,na.rm = TRUE), BIC)]]$allresults$clusterlabels
  BICMessage <- NA
  
  if (max(BICmodel_labels)!=BICmodel){
    BICmodel <- max(BICmodel_labels)
    BICMessage<-"Spurious or empty cluster resulted."
  }
  
  BICresults <- list(allBICvalues=BIC,
    BICmodelselected=BICmodel,
    BICmodelselected_labels=BICmodel_labels,
    BICMessage=BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}  

# ICL function
ICL_function <- function(bIc, gmax, gmin, run){
  ICL <- vector()
  for (g in 1:(gmax-gmin+1)){
    z <- run[[g]]$allresults$probaPost
    mapz <- mclust::unmap(run[[g]]$allresults$clusterlabels)
    forICL <- function(g){sum(log(z[which(mapz[,g]==1),g]))}
    ICL[g] <- bIc$allBICvalues[g] + sum(sapply(1:ncol(mapz),forICL))
  }
  ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]
  ICLmodel_labels <- run[[grep(min(ICL, na.rm = TRUE), ICL)]]$allresults$clusterlabels
  ICLMessage <- NA
  
  if (max(ICLmodel_labels)!=ICLmodel){
    ICLmodel <- max(ICLmodel_labels)
    ICLMessage<-"Spurious or empty cluster resulted."
  }
  
  ICLresults <- list(allICLvalues=ICL,
    ICLmodelselected=ICLmodel,
    ICLmodelselected_labels=ICLmodel_labels,
    ICLMessage=ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
}

# AIC function
AIC_function <- function(ll, k, run, gmin, gmax){
  AIC <- -2*ll+ 2*k
  AICmodel <- seq(gmin, gmax, 1)[grep(min(AIC,na.rm = TRUE), AIC)]
  AICmodel_labels <- run[[grep(min(AIC,na.rm = TRUE), AIC)]]$allresults$clusterlabels
  AICMessage <- NA
  
  if (max(AICmodel_labels)!=AICmodel){
    AICmodel <- max(AICmodel_labels)
    AICMessage<-"Spurious or empty cluster resulted."
  }
  
  AICresults <- list(allAICvalues=AIC,
    AICmodelselected=AICmodel,
    AICmodelselected_labels=AICmodel_labels,
    AICMessage=AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}  

# AIC3 function
AIC3_function <- function(ll, k, run, gmin, gmax){
  AIC3 <- -2*ll+ 3*k
  AIC3model <- seq(gmin, gmax, 1)[grep(min(AIC3,na.rm = TRUE), AIC3)]
  AIC3model_labels <- run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$allresults$clusterlabels
  AIC3Message <- NA
  
  if (max(AIC3model_labels)!=AIC3model){
    AIC3model <- max(AIC3model_labels)
    AIC3Message<-"Spurious or empty cluster resulted."
  }
  AIC3results <- list(allAIC3values=AIC3,
    AIC3modelselected=AIC3model,
    AIC3modelselected_labels=AIC3model_labels,
    AIC3Message=AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}

# main function
MVPLNClustering <- function(dataset, membership=NA, Gmin, Gmax, n_chains=3, n_iterations=NA, init_method="kmeans", n_init_iterations=5){
  
  ptm <- proc.time() 
  
  if (typeof(unlist(dataset)) != "double" & typeof(unlist(dataset)) != "integer"){
    stop("Dataset type needs to be integer");}
  
  if(n_iterations<40){
    stop("RStan n_iterations argument should be greater than 40");}
  
  if((is.na(n_init_iterations) != TRUE && n_init_iterations == !0) && is.na(init_method) == TRUE){
    stop("Number of initialization iterations specified, but no initialization method selected");}
  
  n = length(dataset)
  p = ncol(dataset[[1]])
  r = nrow(dataset[[1]])
  d = p*r
  
  if(all(is.na(membership)!=TRUE) && length(membership)!=n){
    stop("Length of membership character vector and sample size of dataset should match");}
  
  if(all(is.na(membership)!=TRUE) && all((diff(sort(unique(membership)))==1)!=TRUE) ){
    stop("Cluster memberships in the membership vector are missing a cluster, e.g. 1,3,4,5,6 is missing cluster 2");}
  
  # Changing dataset into a n x rp dataset
  normfactor_dataset = matrix(NA,ncol=d, nrow=n)
  sample_matrix = matrix(c(1:d),nrow=r, byrow=TRUE)
  for (u in 1:n){
    for (e in 1:p){
      for (s in 1:r){
        normfactor_dataset[u,sample_matrix[s,e]] = dataset[[u]][s,e]
      }
    }
  }
  
  # Check if entire row is zero
  if(length(which(apply(normfactor_dataset, 1, function(x) all(x==0))==TRUE))!=0){
    #cat("\nDataset row(s)", c(which(apply(normfactor_dataset, 1, function(x) all(x==0))==TRUE)), "will be removed as this/these contain(s) all zeros")
    if(all(is.na(membership)==FALSE)){membership = membership[-c(which(apply(normfactor_dataset, 1, function(x) all(x==0))==TRUE))]}
    normfactor_dataset = normfactor_dataset[-c(which(apply(normfactor_dataset, 1, function(x) all(x==0))==TRUE)),]
    n = nrow(normfactor_dataset)
  }
  
  if(all(is.na(membership)==TRUE)){
    membership <- "Not provided"}
  
  # loading needed packages
  if (!require(mvtnorm)) install.packages("mvtnorm") 
  
  if (!require(mclust)) install.packages("mclust") 
  
  if (!require(mclust)) install.packages("coda") 
  
  if (!require(mclust)) install.packages("capushe") 
  
  if (!require(edgeR)) install.packages("edgeR") 
  
  if (!require(clusterGeneration)) install.packages("clusterGeneration") 
  
  if (!require(pheatmap)) install.packages("pheatmap") 
  
  if (!require(RColorBrewer)) install.packages("RColorBrewer") 
  
  if (!require(gplots)) install.packages("gplots") 
  
  # Generating normalization factors
  norm_factors <- log(as.vector(calcNormFactors(as.matrix(normfactor_dataset), method = "TMM")))

  # Making RStan model 
  if (!require(rstan)) {
    install.packages('rstan') 
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores()) }
  
  if (!require(Rcpp)) install.packages('Rcpp') 
  
  mod <<- stan_model("MPLN.stan")
  
  if (!require(parallel)) install.packages("parallel") 
  
  # Running code in parallel
  # Calculate the number of cores
  no_cores = detectCores()
  
  # Initiate cluster
  cl = makeCluster(no_cores-1) 
  
  # Doing clusterExport
  clusterExport(cl,c("mod", "testing_dataset", "initialization_function", "initialization_run", "zvalue_calculation", "calc_loglikelihood", "parameter_estimation", "stanrun", "cluster_mpln", "calling_clustering", "calculate_parameters", "BIC_function", "ICL_function", "AIC_function", "AIC3_function"))
  
  # Packages need to be downloaded using clusterEvalQ
  clusterEvalQ(cl, library(rstan))
  clusterEvalQ(cl, library(Rcpp))
  clusterEvalQ(cl, library(mclust))
  clusterEvalQ(cl, library(mvtnorm))
  clusterEvalQ(cl, library(edgeR))
  clusterEvalQ(cl, library(capushe))
  clusterEvalQ(cl, library(clusterGeneration))
  clusterEvalQ(cl, library(coda))
  
  MPLN_parallel = function(g){
    ## ** Never use set.seed(), use clusterSetRNGStream() instead,
    # to set the cluster seed if you want reproducible results
    #clusterSetRNGStream(cl=cl, iseed=g)
    test = calling_clustering(n=n, r=r, p=p, d=d, dataset=normfactor_dataset, Gmin=g, Gmax=g, n_chains=n_chains, n_iterations=n_iterations, init_method=init_method, n_init_iterations=n_init_iterations, norm_factors=norm_factors, mod=mod)
    return(test)
  }
  # empty list to save output
  parallel.Wei_2 = list()
  cat("\nRunning parallel code now.")
  parallel.Wei_2 = clusterMap(cl=cl,fun=MPLN_parallel, g=Gmin:Gmax)
  cat("\nDone parallel code.")
  
  stopCluster(cl)
  
  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- k <- ll <- vector()
  
  for(g in 1:(Gmax-Gmin+1)) {
    # save the final log-likelihood
    ll[g]<-unlist(tail(parallel.Wei_2[[g]]$allresults$loglikelihood, n=1)) 
    
    k[g]<-calculate_parameters(g,dataset)
    
    # starting model selection
    if (g==max(1:(Gmax-Gmin+1))){ 
      bic <- BIC_function(ll=ll,k=k, n=n, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
      icl <- ICL_function(bIc=bic, gmin=Gmin, gmax=Gmax, run=parallel.Wei_2)
      aic <- AIC_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax )
      aic3 <- AIC3_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
    }
  }
  
  
  
  final = proc.time()-ptm
  RESULTS = list(dataset = dataset,
    units = r,
    variables = p,
    occassions = n,
    normalization_factors = norm_factors,
    Gmin = Gmin,
    Gmax = Gmax,
    initalization_method = init_method,
    allresults = parallel.Wei_2,
    loglikelihood = ll, 
    numbofparameters = k,
    truelabels = membership,
    ICL.all = icl,
    BIC.all = bic,
    AIC.all = aic,
    AIC3.all = aic3,
    totaltime = final)
  
  class(RESULTS) = "MVPLN"
  return(RESULTS)
}
