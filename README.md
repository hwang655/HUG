# HUG

**author**: `Jianyu Liu`, `Haodong Wang`

# Overview

**HUG** is a data analysis package for estimating Poisson and Hurdle graphical models. The proposed graphical model is designed to be appropriately modeling sample dependence and possible zero-inflation in the data. The sample dependence means there is dependence across the units of analysis (e.g. individuals, cells). The zero-inflation refers to the case when the number of zeros is so large that the data do not readily fit standard distributions. It contains two functions.
- _MPoisGraph_, which fit a dependent Poisson graphical model via nodewise random effect Poisson regression.
- _MHurdGraph_, which Fit a dependent Hurdle graphical model via nodewise random effect Hurdle regression.

## Prerequisite

 HUG package requires R version 3.6.0 or higher. It requires the packages ```Matrix, foreach, doParallel, qlcMatrix, fields, mvtnorm, glmnet, huge, igraph, fields, doSNOW```

## Installation

The `HUG` package is currently available using devtools

``` r
# install.packages('devtools')
devtools::install_github("hwang655/HUG")
```
 
## Example
  An example use:
  ```
library(HUG)

mb.dep.pois = MPoisGraph(Y, Y.p, Omega = NULL, B.ini, lambda.min.ratio=5e-3)
mb.dep.hurd = MHurdGraph(Y, Y.p, Omega = NULL, B.ini, lambda.min.ratio=5e-3)
  ```
where 
* ```Y``` is the response of nodewise regressions, e.g. scRNA-seq data matrix with rows for cells and columns for genes, 
* ```Y.p```	is the predictors of nodewise regressions, which has the same dimension as Y.r.
* ```B.ini``` is a reasonable initial coefficient matrix for the nodewise Hurdle regression
* ```lambda.min.ratio``` is ratio between the minimum and the maximum of the lambda sequence
* ```Omega``` is the precision matrix of the sample dependence model
 
The output is a list with five elements

* ```lambda``` is the the lambda sequence used in nodewise regressions
* ```graphs``` is a sequence of estimated graphs
* ```coef.opt``` is a pxp matrix of EBIC-selected coefficient estimates of all nodewise regressions
* ```coef.aic``` is a pxp matrix of AIC-selected coefficient estimates of all nodewise regressions
* ```time``` is the time spent on each regression in second

## Help and Support

For questions, issues or feature requests please reach out to Haodong: <Haodong@ad.unc.edu>.

## Citation
  If you use the software, please cite our paper: Liu, J., Wang, H. Sun, W., & Liu, Y. (2021). Prioritizing Autism Risk Genes using Personalized Graphical Models Estimated from Single Cell RNA-seq Data
                              
## Reference

Liu, J., Wang, H. Sun, W., & Liu, Y. (2021). Prioritizing Autism Risk Genes using Personalized Graphical Models Estimated from Single Cell RNA-seq Data
