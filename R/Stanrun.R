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
