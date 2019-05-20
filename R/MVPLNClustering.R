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
