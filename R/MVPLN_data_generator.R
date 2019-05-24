
#### function ####
Datagenerator <- function(i, r, p, n, pi_g, mu, phi, omega){
  
  set.seed(i)
  z <- t(rmultinom(n,size=1,pi_g))
  Y = Y2 = norms = X = list()
  
  y <- theta <- n_g <- vector("list", length = length(pi_g)) 
  
  # generating theta
  for (g in 1:length(pi_g)){
    n_g[[g]]<-which(z[,g]==1)
    for (u in 1:length(n_g[[g]]) ){
      # generating V
      theta[[n_g[[g]][u]]]<-matrix(rmvnorm(1, sigma=kronecker(phi[((g-1)*r+1):(g*r),], omega[((g-1)*p+1):(g*p),]) ),r,p)
      
      # adding M1 
      theta[[n_g[[g]][u]]] =  theta[[n_g[[g]][u]]] + mu[((g-1)*r+1):(g*r),] 
    }
  }
  
  # doing exp(theta) to get counts
  for (j in 1:n){
    Y[[j]] = Y2[[j]] = matrix(ncol=p, nrow=r)
    for (i in 1:r){
      for (k in 1:p){
        Y[[j]][i,k]<-rpois(1,exp(theta[[j]][i,k])) 
      }
    }
  }
  
  unlisting_dataset = matrix(NA,ncol=r*p, nrow=n)
  for (u in 1:n){
    unlisting_dataset[u,1:3] = Y[[u]][1,]
    unlisting_dataset[u,4:6] = Y[[u]][2,]
  }
  
  
  # now that data is generated, add norm factors
  norms = matrix(log(as.vector(calcNormFactors(as.matrix(unlisting_dataset), method = "TMM"))), nrow=r, ncol=p, byrow = TRUE) # added one because if the entire column is zero, this gives an error
  
  for (j in 1:n){    
    for (i in 1:r){
      for (k in 1:p){
        Y2[[j]][i,k]<-rpois(1,exp(theta[[j]][i,k]+norms[i,k])) 
      }
    }
  } 
  
  
  results <- list(dataset=Y2,
    truemembership=map(z),
    occassions = n,
    units = r,
    variables = p,
    pi_g = pi_g,
    means = mu,
    phi = phi,
    psi = omega)
  class(results) <- "MVPLN_datagenerator"
  return(results)
  
  # Developed by Anjali Silva
}  

