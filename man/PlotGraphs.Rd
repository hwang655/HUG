\name{PlotGraphs}
\alias{PlotGraphs}
\title{Plot multiple undirected graphs in one plot}
\description{Plot multiple undirected graphs in one plot and display overlapping}
\usage{
PlotGraphs(icovs, graph4layout = NULL, cmpl.layout = NULL, v.keep = NULL, v.name = NULL, v.label = NULL, rmd = TRUE, name.icovs = NULL, legend = TRUE, layout.only = FALSE, conflict.disp = FALSE, ...)
}
\arguments{
  \item{icovs}{a list of graphs in the form of nxn adjacency matrices}
  \item{graph4layout}{(optional) a graph based on which the node positions are determined}
  \item{cmpl.layout}{(optional) an nx2 matrix containing the coordinates of all n nodes}
  \item{v.keep}{(optional) an index array of nodes to be displayed. By default, only nodes with >=1 connection are kept}
  \item{v.name}{(optional) a character array of node names}
  \item{v.label}{(optional) a character array of node symbols}
  \item{rmd}{open a new plotting window if \code{FALSE}. \code{TRUE} by default}
  \item{name.icovs}{(optional) a character array of graph names}
  \item{legend}{whether to plot a legend}
  \item{layout.only}{do not display the graphs if \code{TRUE}. \code{FALSE} by default}
  \item{conflict.disp}{display edges of opposite signs if \code{TRUE}. \code{FALSE} by default}
}
\value{An nx2 matrix containing the coordinates of all n nodes}
\author{Jianyu Liu}
