#' Version of BIC for PCA based on paper by Rajan, Rayner
#'
#' Computes the value of BIC criterion for given data set and
#' number of factors. Assumes uniform distribution of coefficients.
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @return BIC value of BIC criterion
rajan.uniform.BIC <- function(X, k){
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  eig <- eigen(cov(t(X)))
  lambda <- eig$values
  U <- eig$vectors
  Uk <- U[,1:k, drop=F]
  v <- sum(lambda[(k+1):d])/(d-k)

  beta <- max(abs(t(Uk)%*%X))
  w1 <- sqrt(2*v)*(beta - t(Uk)%*%X)
  w2 <- sqrt(2*v)*(-beta - t(Uk)%*%X)

  t0 <- -(k-N)*d/2*log(2*pi*v)
  t1 <- -d*k*log(2*beta)
  t2 <- -d*log(v)
  t3 <- -N*d/2
  t4 <- sum(log((erf(w1)-erf(w2))/2))
  pen <- -(m+d+1+1)/2*log(N)
  t0+t1+t2+t3+t4+pen
}


#' Penalized likelihood criterion for PCA based on paper by Rajan, Rayner
#'
#' Computes the value of BIC criterion for given data set and
#' number of factors. Assumes normal distribution of coefficients.
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @return BIC value of BIC criterion
rajan.BIC <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  v <- sum(lambda[(k+1):d])/(d-k)

  t0 <- -N*d/2*log(2*pi)
  t1 <- -N*k/2*log(mean(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -N*d/2
  pen <- -(m+d+1+1)/2*log(N)
  t0+t1+t2+t3+pen
}


#' Non-penalized likelihood criterion for PCA based on paper by Rajan, Rayner
#'
#' Computes the value of BIC criterion for given data set and
#' number of factors. Assumes normal distribution of coefficients.
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @references Automatic choice of dimensionality for PCA, Thomas P. Minka
#' @return BIC value of BIC criterion
rajan.noBIC <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  v <- sum(lambda[(k+1):d])/(d-k)

  t0 <- -N*d/2*log(2*pi)
  t1 <- -N*k/2*log(mean(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -N*d/2
  t0+t1+t2+t3
}

#' Penalized likelihood for PCA
#'
#' Computes the value of BIC-like criterion for given data set and
#' number of factors. Assumes that number of variables is large
#' compared to number of observations
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @return BIC value of BIC criterion
pca.new.BIC <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  v <- sum(lambda[(k+1):d])/(d-k)

  t0 <- -N*d/2*log(2*pi)
  t1 <- -N/2*sum(log(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -N*d/2
  pen <- -(m+d+k+1)/2*log(N)
  t0+t1+t2+t3+pen
}


#' Penalized likelihood for PCA
#'
#' Computes the value of BIC-like criterion for given data set and
#' number of factors. Assumes that number of variables is large
#' compared to number of observations
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @return BIC value of BIC criterion
pca.new.BIC.fast <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  # lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  lambda <- svd(X, nu = 0, nv = 0)$d^2/(N-1)
  v <- sum(lambda[(k+1):d])/(d-k)

  t0 <- -N*d/2*log(2*pi)
  t1 <- -N/2*sum(log(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -N*d/2
  pen <- -(m+d+k+1)/2*log(N)
  t0+t1+t2+t3+pen
}


#' BIC for PCA, as given by Minka
#'
#' Computes the value of BIC criterion for given data set and
#' number of factors.
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @references Automatic choice of dimensionality for PCA, Thomas P. Minka
#' @return BIC value of BIC criterion
pca.BIC <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- N*k - k*(k+1)/2

  lambda <- eigen(cov(X), only.values = TRUE)$values
  v <- sum(lambda[(k+1):N])/(N-k)

  -d/2*sum(log(lambda[1:k])) -d*(N-k)/2*log(v) -(m+k)/2*log(d)
}


#' Laplace evidence for PCA, as given by Minka
#'
#' Computes the value of Laplace approximation for given data set and
#' number of factors.
#'
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @references Automatic choice of dimensionality for PCA, Thomas P. Minka
#' @return L value of Laplace evidence
pca.Laplace <- function(X, k, alfa=1){
  X <- t(X)
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2

  lambda <- abs(eigen(cov(t(X)), only.values = TRUE)$values)
  v <- sum(lambda[(k+1):d])/(d-k)

  t1 <- -N/2*sum(log(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -k/2*log(N)
  Az <- sapply(1:k, function(i) sum( log(1/lambda[(i+1):d] - 1/lambda[i] ) + log(lambda[i] - lambda[(i+1):d]) + log(N) ))
  if( any(is.nan(Az)) )
    warning(paste("Number of observations ", N, " is to little compared to number of variables ", d,
                  " to perform a meaningful estimation"))
  t4 <- sum(Az)*(-1)/2
  t5 <- log(2*pi)*(m+k)/2
  t6 <- -k*log(2) + sum( lgamma( (d-1:k+1)/2 ) - (d-1:k+1)/2*log(pi) )
  t1+t2+t3+t4+t5+t6
}

#' BIC for subspace clustering
#'
#' Computes the value of BIC criterion for given data set and partition.
#' In each cluster we assume that variables are spanned by few factors.
#' Considering maximum likelihood we get that those factors are in fact
#' principal components. Noise sigma can be computed jointly for all clusters (default),
#' seperately for each cluster or be specified as input.
#'
#'
#' @param X a matrix with only continuous variables
#' @param segmentation a vector, segmentation for which likelihood is computed. Clusters
#'        numbers should be from range [1, numb.clusters]
#' @param dims a vector of integers, dimensions of subspaces. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param max.dim an integer, upper bound for allowed dimension of subspace
#' @param flat.prior boolean, if TRUE (default is FALSE) then flat prior on models is used
#' @keywords internal
#' @return BIC value of BIC criterion
cluster.pca.BIC <- function(X, segmentation, dims, numb.clusters, max.dim, flat.prior = FALSE){
  if(!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]

  formula <- rep(0, numb.clusters)
  for(k in 1:numb.clusters){
    #one cluster
    dim1 = dims[k]
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2]>dim1){
      formula[k] <- pca.new.BIC(Xk, dim1)
    } else{
      formula[k] <- -Inf
    }
  }
  #apriori
  apriori.segmentations <- -p*log(numb.clusters)
  apriori.dimensions <- - log(max.dim)*numb.clusters
  if (flat.prior){
    BIC <- sum(formula)
    return(BIC)
  }
  BIC <- sum(formula) + apriori.segmentations + apriori.dimensions
  return(BIC)
}


#' BIC for subspace clustering
#'
#' Computes the value of BIC criterion for given data set and partition.
#' In each cluster we assume that variables are spanned by few factors.
#' Considering maximum likelihood we get that those factors are in fact
#' principal components. Noise sigma can be computed jointly for all clusters (default),
#' seperately for each cluster or be specified as input.
#'
#'
#' @param X a matrix with only continuous variables
#' @param segmentation a vector, segmentation for which likelihood is computed. Clusters
#'        numbers should be from range [1, numb.clusters]
#' @param dims a vector of integers, dimensions of subspaces. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param adjustment a numeric, percentage of BIC penalty applied
#' @param sigma a numeric, (default is NULL) value of sigma provided by the user
#' @param estimateJointly a boolean, (default value is TRUE) indicating if sigma should be estimated jointly for all clusters
#' @keywords internal
#' @return BIC value of BIC criterion
adjusted.cluster.BIC <- function(X, segmentation, dims, numb.clusters, adjustment = 1, sigma=NULL, estimateJointly=TRUE){
  if(!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if(is.null(sigma) & estimateJointly){
    RES.sigma <- sum(vapply(1:numb.clusters, function(k) {
      max.dim = dims[k]
      Xk = X[,segmentation==k, drop=F];
      if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
        svdSIGNAL= svd(Xk);
        SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*%
          diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*%
          t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim));
        return(sum((Xk - SIGNAL)^2))
      }
      return(0)
    }, 0.9))
    max.dim = max(dims)
    degrees.freedom <- D*p-p-numb.clusters*D*max.dim-p*max.dim+numb.clusters*max.dim^2+numb.clusters*max.dim
    #     print(degrees.freedom)
    #     ps <- as.numeric(table(segmentation))
    #     degrees.freedom <- sum(D*ps-ps-numb.clusters*D*dims-ps*dims+numb.clusters*dims^2+numb.clusters*dims)
    #     print(degrees.freedom)
    sigma <- sqrt(RES.sigma/degrees.freedom)
  }
  likelihoods <- rep(0, numb.clusters)
  penalties <- rep(0, numb.clusters)
  for(k in 1:numb.clusters){
    #one cluster
    max.dim = dims[k]
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*%
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*%
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RESIDUAL = Xk - SIGNAL
      if(!estimateJointly & is.null(sigma)) {
        df <- D*ncol(Xk)-ncol(Xk)-D*max.dim-ncol(Xk)*max.dim+max.dim^2+max.dim
        sigma = sqrt(sum(RESIDUAL^2)/((D-1)*(ncol(Xk)-1)))
        sigma <- sqrt(sum(RESIDUAL^2)/df)
      }
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
      mk <- ncol(Xk)
      penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
      if(!estimateJointly)
        penalties[k] <- penalties[k] + 1/2*log(mk) #for estimating sigma
    } else{
      if(!estimateJointly & is.null(sigma)){
        likelihoods[k] <- -Inf
        penalties[k] <- 0
      } else{
        RESIDUAL = Xk - Xk
        likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
        mk <- max(ncol(Xk),1)
        penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
      }
    }
  }
  BIC <- sum(likelihoods) - adjustment*sum(penalties)
  return(BIC)
}


#' BIC for subspace clustering
#'
#' Computes the value of BIC criterion for given data set and partition.
#' In each cluster we assume that variables are spanned by few factors.
#' Considering maximum likelihood we get that those factors are in fact
#' principal components. Noise sigma can be computed jointly for all clusters (default),
#' seperately for each cluster or be specified as input.
#'
#'
#' @param X a matrix with only continuous variables
#' @param segmentation a vector, segmentation for which likelihood is computed. Clusters
#'        numbers should be from range [1, numb.clusters]
#' @param max.dim an integer, maximum dimension of subspace. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param sigma a numeric, (default is NULL) value of sigma provided by the user
#' @param estimateJointly a boolean, (default value is TRUE) indicating if sigma should be estimated jointly for all clusters
#' @keywords internal
#' @return BIC value of BIC criterion
cluster.BIC <- function(X, segmentation, max.dim, numb.clusters, sigma=NULL, estimateJointly=TRUE){
  if (!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if (is.null(sigma) & estimateJointly){
    RES.sigma <- sum(vapply(1:numb.clusters, function(k) {
      Xk = X[,segmentation==k, drop=F];
      if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
        svdSIGNAL= svd(Xk);
        SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*%
          diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*%
          t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim));
        return(sum((Xk - SIGNAL)^2))
      }
      return(0)
    }, 0.9))
    degrees.freedom <- D*p-p-numb.clusters*D*max.dim-p*max.dim+numb.clusters*max.dim^2+numb.clusters*max.dim
    sigma <- sqrt(RES.sigma/degrees.freedom)
  }
  likelihoods <- rep(0, numb.clusters)
  penalties <- rep(0, numb.clusters)
  for(k in 1:numb.clusters){
    #one cluster
    Xk = X[,segmentation==k, drop=F]
    if (dim(Xk)[2]>max.dim){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*%
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*%
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RESIDUAL = Xk - SIGNAL
      if(!estimateJointly & is.null(sigma)){ sigma = sqrt(sum(RESIDUAL^2)/((D-1)*(ncol(Xk)-1))) }
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
      mk <- ncol(Xk)
      penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
    }
    else{
      RESIDUAL = Xk - Xk
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
      mk <- max(ncol(Xk),1)
      penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
    }
  }
  BIC <- sum(likelihoods) - sum(penalties)
  return(BIC)
}


##' Computing joint sigma for all clusters
##'
##' Computes unbiased noise estimator under assumption that
##' all subspaces are of the same dimension
##'
##' @param X data
##' @param segmentation variable segmentation
##' @param max.dim maximal subspace dimension
##' @param n number of rows
##' @param p number of variables
##' @param numb.clusters total number of clusters
##' @return unbiased noise estimator
##' @keywords internal
getSigma <- function(X, segmentation, max.dim, n, p, numb.clusters){
  RES.sigma=0
  for(k in 1:numb.clusters){
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*%
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*%
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RES.sigma = RES.sigma + sum((Xk - SIGNAL)^2)
    }
  }
  degrees.freedom <- n*p-p-n*max.dim-p*max.dim+max.dim^2+max.dim
  sigma <- sqrt(RES.sigma/degrees.freedom)
  sigma
}

#' Selects subspace closest to given variable (according to BIC)
#'
#' The most similar subspace is choosen based on BIC criterion
#'
#' @param variable variable to be assigned
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @keywords internal
#' @return index number of subspace closest to variable
choose.cluster.BIC <- function(variable, pcas, numberClusters){
  BICs <- NULL
  for(i in 1:numberClusters){
    nparams <- ncol(pcas[[i]])
    n <- length(variable)
    res <- fastLmPure(pcas[[i]], variable, method = 0L)$residuals
    sigma.hat <- sqrt(sum(res^2)/100)
    if (sigma.hat < 1e-15)
      warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-16. It might corrupt the result.")
    loglik <- sum(dnorm(res, 0, sigma.hat, log=T))
    BICs[i] <- loglik - nparams*log(n)/2
  }
  which.max(BICs)
}

#' Choose subspace closest to the given variable
#'
#' The most similar subspace is choosen based on R^2
#'
#' @param variable variable to be assigned
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @return index number of subspace closest to variable
#' @keywords internal
choose.cluster <- function(variable, pcas, numberClusters){
  v1 = var(variable)
  which.max( vapply(1:numberClusters, function(i){
    v2 <- var(fastLmPure(pcas[[i]], variable, method = 0L)$residuals);
    p <- ncol(pcas[[i]]);
    n <- length(variable);
    (v1-v2)/v1
  }, 0.9) )
}

#' Plot mlcc.fit class object
#'
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.mlcc.fit <- function(x,...){
  clusterNumbs <- lapply(x$all.fit, function(y) y$nClusters)
  BICVals <- lapply(x$all.fit, function(y) y$BIC)
  plot.default(clusterNumbs, BICVals, type="b", xaxt="n", ylab="BIC", xlab="Number of clusters")
  axis(side = 1, labels = clusterNumbs, at=clusterNumbs)
}

#' Print mlcc.fit class object
#'
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.mlcc.fit <- function(x,...){
  cat("$nClusters: ", x$nClusters, "\n")
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$subspacesDimensions:\n", unlist(x$subspacesDimensions), "\n")
}

#' Print mlcc.reps.fit class object
#'
#' @param x mlcc.reps.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.mlcc.reps.fit <- function(x,...){
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$basis:\n")
  cat(str(x$basis))
}
