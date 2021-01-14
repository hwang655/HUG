\name{rMVH}
\alias{rMVH}
\title{Generate Hurdle random vectors with sample dependence}
\description{
Generate Hurdle random vectors with sample dependence
}
\usage{
rMVH(coef.hurd, mu, size, g.type, factor.re = 1)
}
\arguments{
  \item{coef.hurd}{\eqn{gamma_0} and \eqn{gamma_1} of the Hurdle model}
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

\seealso{rMVP}
\examples{
Y = rMVH(c(0, 0.5), mu=rep(2, 100), size=80, g.type='hub', factor.re=.5)
}
