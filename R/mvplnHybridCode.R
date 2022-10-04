
mvplnHybridclus <- function(dataset,
                         membership = "none",
                         gmin,
                         gmax,
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



}
