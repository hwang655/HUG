\name{hugeGraph}
\alias{hugeGraph}
\title{
A wrapper function for huge::huge
}
\description{
A wrapper function for \code{huge::huge} function that works for data with constance columns.
}
\usage{
hugeGraph(Y, ngraph = 50, method = "glasso", ...)
}
\arguments{
  \item{Y}{an nxp matrix containing the data}
  \item{ngraph}{number of lambda values}
  \item{method}{"glasso" or "mb"}
  \item{\dots}{passed to \code{huge}}
}
\value{A \code{huge} object}
\author{Jianyu Liu}

\seealso{
GLMGraph, MPoisGraph, MHurdGraph
}
