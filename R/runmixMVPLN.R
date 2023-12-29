#' Launch Shiny App For Package mixMVPLN
#'
#' A function that launches the shiny app for this package.
#' The shiny app permit to perform clustering using mixtures
#' of matrix variate Poisson-log normal (MVPLN) via
#' variational Gaussian approximations. Model selection
#' can be done using AIC, AIC3, BIC and ICL.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runMPLNClust()
#' }
#'
#' @author {Anjali Silva, \email{anjali@alumni.uoguelph.ca}}
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Silva, A. et al. (2018). Finite Mixtures of Matrix Variate Poisson-Log Normal Distributions
#' for Three-Way Count Data. \href{https://arxiv.org/abs/1807.08380}{arXiv preprint arXiv:1807.08380}.
#'
#' @export
#' @importFrom shiny runApp
#' @import shinyalert
runmixMVPLN <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "mixMVPLN")
  shinyalert::closeAlert(num=0,id=NULL)
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
