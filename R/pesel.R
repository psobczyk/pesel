#' Automatic estimation of number of principal components in PCA
#' with PEnalized SEmi-integrated Likelihood (PESEL)
#'
#' Underlying assumption is that only small number of principal components,
#' associated with largest singular values, is relevent, while the rest of them
#' is noise. For a given numeric data set, function estimates the number of PCs
#' according to penalized likelihood criterion. Function adjusts the model used
#' to the case when number of variables is larger than the number of
#' observations.
#'
#' Please note that no categorical variables and missing values are allowed.
#'
#'
#' @param X a data frame or a matrix contatining only continuous variables
#' @param npc.min minimal number of principal components, for all the possible
#' number of PCs between npc.min and npc.max criterion is computed
#' @param npc.max maximal number of principal components, if greater than
#' dimensions of X, min(ncol(X), nrow(X))-1 is used, for all the possible
#' number of PCs between npc.min and npc.max criterion is computed
#' @param method name of criterion to be used
#' @param scale a boolean, if TRUE (default value) then data is scaled before
#' applying criterion
#' @export
#' @return number of components
#' @examples
#' \dontrun{
#' library(MetabolAnalyze)
#' data(UrineSpectra)
#' pesel(UrineSpectra[[1]], method = "heterogenous")
#' }
pesel <- function(X, npc.min = 1, npc.max = 10, scale = FALSE,
                      method = c("heterogenous", "homogenous")){
  # preprocessing on X
  # number of components must be smaller than dimensions of X
  n = nrow(X)
  p = ncol(X)
  npc.max = min(npc.max, min(n,p)-1)
  npc.min = max(npc.min, 1)

  method = match.arg(method)

  if(class(X) == "data.frame"){
    X = as.matrix(X)
  }

  if(sum(sapply(X, is.numeric)) < p){
    stop("All the variables have to be numeric")
  }

  missing = which(is.na(X))
  if(length(missing) !=  0){
    stop("There are missing values")
  }

  if(scale)
    X = scale(X)

  vals = numeric(10)
  vals = if(p>n) {
    switch(method,
                 "heterogenous" = pesel_heterogeneous_big_p(X, npc.min, npc.max),
                 "homogenous" = pesel_homogeneous_big_p(X, npc.min, npc.max))
    } else {
      switch(method,
             "heterogenous" = pesel_heterogeneous_big_n(X, npc.min, npc.max),
             "homogenous" = pesel_homogeneous_big_p(t(X), npc.min, npc.max))
    }

  result = NULL
  result$nPCs = npc.min-1+which.max(vals)
  result$vals = vals
  result$npc.min = npc.min
  result$npc.max = npc.max
  class(result) = "pesel.result"
  result
}


#' Plot pesel.result class object
#'
#' @param x pesel.result class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.pesel.result <- function(x,...){
  vals = x$vals - max(x$vals) + 20
  probs = exp(vals)/sum(exp(vals))
  plot(x$npc.min:x$npc.max, probs, xlab = "Number of components",
       ylab = "Posterior probability", main = "PESEL", type = "b")
  points(x$npc.min-1+which.max(vals), max(probs), col = "red")
}

#' Print pesel.result class object
#'
#' @param x pesel.result class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.pesel.result <- function(x,...){
  cat("$nPCs: ", x$nPCs, "\n")
  cat("$vals: value of PeSeL criterion\n")
  cat("$npc.min: ", x$npc.min, "\n")
  cat("$npc.max: ", x$npc.max, "\n")
}
