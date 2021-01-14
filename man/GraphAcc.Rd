\name{GraphAcc}
\alias{GraphAcc}
\title{
Evaluate the accuracy of graph estimation
}
\description{
Evaluate the accuracy of a list of graph estimates.
}
\usage{
GraphAcc(graphs, graph.t)
}
\arguments{
  \item{graphs}{a list of graph estimates in the form of adjacency matrices}
  \item{graph.t}{the benchmark graph in the form of adjacency matrix}
}
\details{No symmetrization will be performed on the graph estimates \code{graphs}; \code{graph.t} will be symmetrized.}
\value{A data.frame of accuracy evaluation results, including true positives (TP), false positives (FP), false negatives (FN), number of nonzero values in each estimated adjacency matrix (DF), number of nonzero values in the benchmark adjacency matrix (DF_T). Each row corresponds to an estimated graph in \code{graphs}.}
\author{Jianyu Liu}
