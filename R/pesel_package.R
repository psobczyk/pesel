#' Automatic estimation of number of principal components in PCA
#'
#' @description Automatic estimation of number of principal components in PCA
#' with PEnalized SEmi-integrated Likelihood (PESEL).
#'
#' @details Version: 0.6.0
#' @docType package
#' @name pesel
#' @importFrom stats cov
#' @importFrom graphics plot points
#' @author{ Piotr Sobczyk,
#'          Julie Josse,
#'          Malgorzata Bogdan
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' }
#' @references \emph{Bayesian dimensionality reduction with PCA using penalized semi-integrated likelihood},
#' Piotr Sobczyk, Malgorzata Bogdan, Julie Josse
#'
#' @examples
#' \dontrun{
#' library(varclust)
#' set.seed(1)
#' sim.data <- data.simulation(n = 100, SNR = 1, K=1, numb.vars = 200, max.dim = 5)
#' pesel(sim.data$X, 1, 10, prior = c(1, rep(0, 9)))
#' }
NULL
