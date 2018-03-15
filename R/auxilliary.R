#' PEnalized SEmi-integrated Likelihood for homogeneous singular values and
#' large number of variables
#'
#' Derived under assumption that number of variables tends to infinity
#' while number of observations is limited.
#'
#' @param X a matrix containing only continuous variables
#' @param minK minimal number of principal components fitted
#' @param maxK maximal number of principal components fitted
#' @export
#' @return numeric vector, PESEL criterion for each k in range [minK, maxK]
pesel_homogeneous <- function(X, minK, maxK){
  d <- dim(X)[1]
  N <- dim(X)[2]
  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  if(any(lambda < 0)){ #numerical error, eigen values of covariance matrix should be positive
    lambda[lambda < 0] = 1e-16 #we change it to something close to zero
  }

  pesel <- NULL
  for(k in minK:maxK){
    m <- d*k - k*(k+1)/2
    v <- sum(lambda[(k+1):d])/(d-k)

    t0 <- -N*d/2*log(2*pi)
    t1 <- -N*k/2*log(mean(lambda[1:k]))
    t2 <- -N*(d-k)/2*log(v)
    t3 <- -N*d/2
    pen <- -(m+d+1+1)/2*log(N)
    pesel <- c(pesel, t0 + t1 + t2 + t3 + pen)
  }
  pesel
}

#' PEnalized SEmi-integrated Likelihood for heterogeneous singular values and
#' large number of variables
#'
#' Derived under assumption that number of variables tends to infinity
#' while number of observations is limited.
#'
#' @param X a matrix containing only continuous variables
#' @param minK minimal number of principal components fitted
#' @param maxK maximal number of principal components fitted
#' @export
#' @return numeric vector, PESEL criterion for each k in range [minK, maxK]
pesel_heterogeneous <- function(X, minK, maxK){
  d <- dim(X)[1]
  N <- dim(X)[2]
  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  if(any(lambda < 0)){ #numerical error, eigen values of covariance matrix should be positive
    lambda[lambda < 0] = 1e-16 #we change it to something close to zero
  }

  pesel <- NULL
  for(k in minK:maxK){
    m <- d*k - k*(k+1)/2
    v <- sum(lambda[(k+1):d])/(d-k)

    t0 <- -N*d/2*log(2*pi)
    t1 <- -N/2*sum(log(lambda[1:k]))
    t2 <- -N*(d-k)/2*log(v)
    t3 <- -N*d/2
    pen <- -(m+d+k+1)/2*log(N)
    pesel <- c(pesel, t0+t1+t2+t3+pen)
  }
  pesel
}
