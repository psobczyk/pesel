% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxilliary.R
\name{pesel_heterogeneous}
\alias{pesel_heterogeneous}
\title{PEnalized SEmi-integrated Likelihood for heterogeneous singular values and
large number of variables}
\usage{
pesel_heterogeneous(X, minK, maxK)
}
\arguments{
\item{X}{a matrix containing only continuous variables}

\item{minK}{minimal number of principal components fitted}

\item{maxK}{maximal number of principal components fitted}
}
\value{
numeric vector, PESEL criterion for each k in range [minK, maxK]
}
\description{
Derived under assumption that number of variables tends to infinity
while number of observations is limited.
}
