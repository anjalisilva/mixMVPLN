
mvplnHybridclus <- function(dataset,
                         membership = "none",
                         gmin,
                         gmax,
                         nChains,
                         nIterations = NA,
                         initMethod = "kmeans",
                         nInitIterations = 0,
                         normalize = "Yes") {

  mvplnVGAclusOutput <-  mvplnVGAclus(dataset = dataset,
                           membership = membership,
                           gmin = gmin,
                           gmax = gmax,
                           initMethod = initMethod,
                           nInitIterations = nInitIterations,
                           normalize = normalize)

  mvplnMCMCclusOutput <- mvplnMCMCclus(dataset = dataset,
                            membership = membership,
                            gmin = gmin,
                            gmax = gmax,
                            nChains = nChains,
                            nIterations = nIterations,
                            initMethod = initMethod,
                            nInitIterations = 0,
                            normalize = normalize,
                            numNodes = numNodes,
                            VGAparameters = mvplnVGAclusOutput)

}
