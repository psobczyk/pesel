#' pesel
#'
#' @docType package
#' @name pesel
#' @details Version: 0.5.0
#' @importFrom stats cov
#' @importFrom graphics plot points
#' @author Piotr Sobczyk,
#'         Julie Josse,
#'         Malgorzata Bogdan
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#'
#' @examples
#' \dontrun{
#' library(varclust)
#' set.seed(1)
#' sim.data <- data.simulation(n = 100, SNR = 1, K=1, numb.vars = 200, max.dim = 5)
#' pesel(sim.data$X, 1, 10, prior = c(1, rep(0, 9)))
#' }
NULL
