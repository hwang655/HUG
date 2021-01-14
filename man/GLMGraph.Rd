\name{GLMGraph}
\alias{GLMGraph}
\title{Fit a graphical model nodewise via generalized linear regression}
\description{
Fit a graphical model nodewise via generalized linear regression
}
\usage{
GLMGraph(Y.r, Y.p, Z = NULL, nlambda = 50, family = "gaussian", intercept = T, lambda.min.ratio = 0.01, weight.zero = 1, hetero.penalty = F)
}
\arguments{
  \item{Y.r}{response of nodewise regressions, an nxp matrix}
  \item{Y.p}{predictors of nodewise regressions, a matrix of the same dimension as \code{Y.r}}
  \item{Z}{extra predictor matrix for all nodewise regressions}
  \item{nlambda}{number of lambda values on grid (default 50)}
  \item{family}{"gaussian" or "poisson"}
  \item{intercept}{whether to include an intercept term in each nodewise regressions}
  \item{lambda.min.ratio}{ratio between the minimum and the maximum of the lambda sequence}
  \item{weight.zero}{Weight of zero-response samples in nodewise regressions}
}
\value{
\item{lambda}{the lambda sequence used in nodewise regressions}
\item{graphs}{a sequence of estimated graphs}
\item{coef.opt}{a pxp matrix of EBIC-selected coefficient estimates of all nodewise regressions}
\item{coef.aic}{a pxp matrix of AIC-selected coefficient estimates of all nodewise regressions}
}
\author{Jianyu Liu}

\seealso{
hugeGraph, MPoisGraph, MHurdGraph
}