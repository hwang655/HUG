\name{rMVP}
\alias{rMVP}
\title{Generate Poisson random vectors with sample dependence}
\description{Generate Poisson random vectors with sample dependence}
\usage{
rMVP(mu, size, g.type, factor.re = 1)
}
\arguments{
  \item{mu}{a p-dimensional mean vector of the underlying Gaussian random vector}
  \item{size}{sample size}
  \item{g.type}{type of the underlying graph, one of "banded", "hub", and "random"}
  \item{factor.re}{inverse scale of the sample dependence}
}
\value{
\item{Y}{an nxp sample matrix, where each row for a sample and each column for a variable}
\item{X}{an nxp matrix of the underlying Gaussian data with sample dependence}
\item{Z}{a list of two Gaussian components}
\item{graph}{the underlying true graph of the random variables}
}
\author{Jianyu Liu}

\seealso{rMVH}
\examples{
Y = rMVP(mu=rep(2, 100), size=80, g.type='hub', factor.re=.5)
}
