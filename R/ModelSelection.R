# AIC function
AICFunction <- function(ll,
                        k,
                        run,
                        gmin,
                        gmax,
                        parallel = TRUE) {
  AIC <- -2 * ll + 2 * k
  AICmodel <- seq(gmin, gmax, 1)[grep(min(AIC,na.rm = TRUE), AIC)]

  if(isTRUE(parallel) == "FALSE"){
    # if non parallel run
    AICmodelLabels <- run[[grep(min(AIC,na.rm = TRUE), AIC)]]$clusterlabels
  }else{
    # if parallel run
    AICmodelLabels <- run[[grep(min(AIC,na.rm = TRUE), AIC)]]$allresults$clusterlabels
  }
  AICMessage <- NA

  if (max(AICmodelLabels) != AICmodel) {
    AICmodel <- max(AICmodelLabels)
    AICMessage <- "Spurious or empty cluster resulted."
  }

  AICresults <- list(allAICvalues = AIC,
                     AICmodelselected = AICmodel,
                     AICmodelselectedLabels = AICmodelLabels,
                     AICMessage = AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}

# AIC3 function
AIC3Function <- function(ll,
                         k,
                         run,
                         gmin,
                         gmax,
                         parallel = TRUE) {
  AIC3 <- - 2 * ll + 3 * k
  AIC3model <- seq(gmin, gmax, 1)[grep(min(AIC3, na.rm = TRUE), AIC3)]
  if(isTRUE(parallel) == "FALSE"){
    # if non parallel run
    AIC3modelLabels <- run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$clusterlabels
  } else{
    # if parallel run
    AIC3modelLabels <- run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$allresults$clusterlabels
  }
  AIC3Message <- NA

  if (max(AIC3modelLabels) != AIC3model) {
    AIC3model <- max(AIC3modelLabels)
    AIC3Message <- "Spurious or empty cluster resulted."
  }
  AIC3results <- list(allAIC3values = AIC3,
                      AIC3modelselected = AIC3model,
                      AIC3modelselectedLabels = AIC3modelLabels,
                      AIC3Message = AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}

# BIC function
BICFunction <- function(ll,
                        k,
                        n,
                        run,
                        gmin,
                        gmax,
                        parallel = TRUE) {
  BIC <- -2 * ll + (k * log(n))
  BICmodel <- seq(gmin, gmax, 1)[grep(min(BIC,na.rm = TRUE), BIC)]
  if(isTRUE(parallel) == "FALSE") {
    # if non parallel run
    BICmodelLabels <- run[[grep(min(BIC, na.rm = TRUE),
                                BIC)]]$clusterlabels
  } else {
    # if parallel run
    BICmodelLabels <- run[[grep(min(BIC, na.rm = TRUE),
                                BIC)]]$allresults$clusterlabels
  }
  BICMessage <- NA

  if (max(BICmodelLabels) != BICmodel) {
    BICmodel <- max(BICmodelLabels)
    BICMessage <- "Spurious or empty cluster resulted."
  }

  BICresults <- list(allBICvalues = BIC,
                     BICmodelselected = BICmodel,
                     BICmodelselectedLabels = BICmodelLabels,
                     BICMessage = BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}

# ICL function
ICLFunction <- function(bIc,
                        gmax,
                        gmin,
                        run,
                        parallel = TRUE) {
  ICL <- vector()
  for (g in 1:(gmax - gmin + 1)) {
    if(isTRUE(parallel) == "FALSE") {
      # if non parallel run
      z <- run[[g]]$probaPost
      mapz <- mclust::unmap(run[[g]]$clusterlabels)
    } else {
      # if parallel run
      z <- run[[g]]$allresults$probaPost
      mapz <- mclust::unmap(run[[g]]$allresults$clusterlabels)
    }
    forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}
    ICL[g] <- bIc$allBICvalues[g] + sum(sapply(1:ncol(mapz),forICL))
  }
  ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]


  if(isTRUE(parallel) == "FALSE") {
    # if non parallel run
    ICLmodelLabels <- run[[grep(min(ICL, na.rm = TRUE),
                                ICL)]]$clusterlabels
  } else {
    # if parallel run
    ICLmodelLabels <- run[[grep(min(ICL, na.rm = TRUE),
                                ICL)]]$allresults$clusterlabels
  }

  ICLMessage <- NA

  if (max(ICLmodelLabels) != ICLmodel){
    ICLmodel <- max(ICLmodelLabels)
    ICLMessage<-"Spurious or empty cluster resulted."
  }

  ICLresults <- list(allICLvalues = ICL,
                     ICLmodelselected = ICLmodel,
                     ICLmodelselectedLabels = ICLmodelLabels,
                     ICLMessage = ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
  # Developed by Anjali Silva
}
