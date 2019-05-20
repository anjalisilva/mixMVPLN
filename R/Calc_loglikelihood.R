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