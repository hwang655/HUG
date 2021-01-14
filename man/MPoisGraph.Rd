\name{MPoisGraph}
\alias{MPoisGraph}
\title{Fit a dependent Poisson graphical model}
\description{
Fit a dependent Hurdle graphical model via nodewise random effect Hurdle regression.
}
\usage{
MPoisGraph(Y.r, Y.p, B.ini, Omega = NULL, nlambda = 50, lambda.max = 10, lambda.min.ratio = 0.01)
}
\arguments{
  \item{Y.r}{response of nodewise regressions, an nxp matrix}
  \item{Y.p}{predictors of nodewise regressions, a matrix of the same dimension as \code{Y.r}}
  \item{B.ini}{a reasonable initial coefficient matrix for the nodewise Hurdle regression}
  \item{Omega}{(optional) precision matrix of the sample dependence model}
  \item{nlambda}{number of lambda values on grid (default 50)}
  \item{lambda.max}{maximum of the lambda sequence}
  \item{lambda.min.ratio}{ratio between the minimum and the maximum of the lambda sequence}
}
\value{
\item{lambda}{the lambda sequence used in nodewise regressions}
\item{graphs}{a sequence of estimated graphs}
\item{coef.opt}{a pxp matrix of EBIC-selected coefficient estimates of all nodewise regressions}
\item{coef.aic}{a pxp matrix of AIC-selected coefficient estimates of all nodewise regressions}
\item{time}{time spent on each regression in second}
}
\author{Jianyu Liu}

\seealso{GLMGraph, hugeGraph, MHurdGraph}