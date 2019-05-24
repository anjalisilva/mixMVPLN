# Calling the clustering function
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
  # Developed by Anjali Silva
}