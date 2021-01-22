# HUG

**author**: `Jianyu Liu`, `Haodong Wang`

This repository contains code needed to reproduce the article:

Liu, J., Wang, H. Sun, W., & Liu, Y. (2021). Prioritizing Autism Risk Genes using Personalized Graphical Models Estimated from Single Cell RNA-seq Data

## Overview
========

**HUG** is a data analysis package for estimating Poisson and Hurdle graphical models. The proposed graphical model is designed to be appropriately modeling sample dependence and possible zero-inflation in the data. The sample dependence means there is dependence across the units of analysis (e.g. individuals, cells). The zero-inflation refers to the case when the number of zeros is so large that the data do not readily fit standard distributions. The primary algorithm in this package is dependent Poisson graphical model and dependent Hurdle graphical model. For a detailed discussion of the algorithm see [Prioritizing Autism Risk Genes using Personalized Graphical Models Estimated from Single Cell RNA-seq Data](https://arxiv.org/).

## Installation
========

The `HUG` package is currently available using devtools

``` r
# install.packages('devtools')
devtools::install_github("hwang655/HUG")
```

Example
=======

``` r
library(HUG)
library(qlcMatrix)
library(huge)
library(PLNet)
library(doParallel)

cl = makeCluster(4)
registerDoParallel(cl)
getDoParWorkers()

set.seed(41)
n.sim = 500
n = 100


# Dependent Sample Generating
set.seed(31381301)
n = 100
p = 80
dat = rMVP(mu=rep(2, p), size=n, g.type='hub', factor.re=.5)                        # HPLN setting
# dat = rMVH(coef.hurd=c(0, .5), mu=rep(3, p), size=n, g.type='hub', factor.re=.5)  # HHLN setting
g.benchmark = dat$graph
Y = dat$Y
logY = log(as.matrix(Y)+1)
npnY = huge.npn(Y)


# Graph Estimation
glasso.log = hugeGraph(logY, 50, 'glasso', verbose=F)
glasso.npn = hugeGraph(npnY, 50, 'glasso', verbose=F)
mb.pois = GLMGraph(Y, logY, NULL, 50, 'poisson', lambda.min.ratio=1e-3)
# mcdavid = fitHurdle(logY, lambda.min.ratio=.4, parallel=T)

# Dependent graphical model estimation (our methods)
B.ini = matrix(0, p+1, p)
for (j in 1:p) B.ini[-(j+1), j] = as.vector(mb.pois$coef.aic[,j])
mb.dep.pois = MPoisGraph(Y, logY, NULL, B.ini, lambda.min.ratio=5e-3)
mb.dep.hurd = MHurdGraph(Y, logY, NULL, B.ini, lambda.min.ratio=5e-3)

# Evaluation
ests = list(glasso=glasso.log$path, glasso.npn=glasso.npn$path, poisson=mb.pois$graphs,  
            # hurdle=mcdavid$adjMat,
            dep.poisson=mb.dep.pois$graphs, dep.hurdle=mb.dep.hurd$graphs)
acc = lapply(ests, GraphAcc, graph.t=g.benchmark)

# Display results
plot(c(0, 0.1), c(0, 1), type='n', xlab='FPR', ylab='TPR')
for (k in seq(length(acc))) {
    with(acc[[k]], lines(FP/(p*(p-1)-DF_T), TP/DF_T, col=k))
    with(acc[[k]], points(FP/(p*(p-1)-DF_T), TP/DF_T, col=k, pch=k-1))
}
legend('bottomright', col=1:6, lty=1, pch=0:5, legend=names(acc))

```

# Help and Support

Additional documentation, examples and code revisions are coming soon.
For questions, issues or feature requests please reach out to Haodong: <Haodong@ad.unc.edu>.

#### Documentation

The source code is located on github: <https://github.com/hwang655/HUG>. 

#### Contributing

We welcome contributions to make this a stronger package: data examples,
bug fixes, spelling errors, new features, etc. 
