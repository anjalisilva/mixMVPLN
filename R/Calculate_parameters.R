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