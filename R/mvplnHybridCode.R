sampleData


mvplnHybridclus(dataset = sampleData$dataset,
                membership = "none",
                gmin = 1,
                gmax = 2,
                nChains = 3,
                nIterations = NA,
                initMethod = "kmeans",
                nInitIterations = 1,
                normalize = "Yes",
                numNodes = NA)


mvplnHybridclus <- function(dataset,
                         membership = "none",
                         gmin,
                         gmax,
                         nChains = 3,
                         nIterations = NA,
                         initMethod = "kmeans",
                         nInitIterations = 0,
                         normalize = "Yes",
                         numNodes = NA) {

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
                            nInitIterations = 0, # must be 0 for hybrid approach
                            normalize = normalize,
                            numNodes = numNodes,
                            VGAparameters = mvplnVGAclusOutput)
  return(mvplnMCMCclusOutput)
}
