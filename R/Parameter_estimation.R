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
  # Developed by Anjali Silva
}

