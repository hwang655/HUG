library(HUG)
library(qlcMatrix)
library(huge)
library(PLNet)
library(doParallel)

cl = makeCluster(4)
registerDoParallel(cl)
getDoParWorkers()

# Zero inflation caused by random effects
set.seed(41)
n.sim = 500
n = 100
y.avg = numeric(n.sim)
y.var = numeric(n.sim)
y.zero = numeric(n.sim)
for (i in 1:n.sim) {
    x = runif(1, -1, 2)
    mu = exp(x + rnorm(n))
    y.i = rpois(n, mu)
    y.avg[i] = mean(y.i)
    y.var[i] = var(y.i)
    y.zero[i] = mean(y.i == 0)
}
par(mfrow=c(1, 2))
plot(y.avg, y.var, xlab='sample mean', ylab='sample variance', pch='.', cex=3)
abline(0, 1, col=2, lty=3)
plot(y.avg, y.zero, xlab='sample mean', ylab='zero proportion', pch='.', cex=3)
lines(sort(y.avg), exp(-sort(y.avg)), col=2, lty=3)


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

# Dependent graphical model estimation (our methods)
B.ini = matrix(0, p+1, p)
for (j in 1:p) B.ini[-(j+1), j] = as.vector(mb.pois$coef.aic[,j])
mb.dep.pois = MPoisGraph(Y, logY, NULL, B.ini, lambda.min.ratio=5e-3)
mb.dep.hurd = MHurdGraph(Y, logY, NULL, B.ini, lambda.min.ratio=5e-3)

# Evaluation
ests = list(glasso=glasso.log$path, glasso.npn=glasso.npn$path, poisson=mb.pois$graphs, 
            dep.poisson=mb.dep.pois$graphs, dep.hurdle=mb.dep.hurd$graphs)
acc = lapply(ests, GraphAcc, graph.t=g.benchmark)

# Display results
plot(c(0, 0.1), c(0, 1), type='n', xlab='FPR', ylab='TPR')
for (k in seq(length(acc))) {
    with(acc[[k]], lines(FP/(p*(p-1)-DF_T), TP/DF_T, col=k))
    with(acc[[k]], points(FP/(p*(p-1)-DF_T), TP/DF_T, col=k, pch=k-1))
}
legend('bottomright', col=1:5, lty=1, pch=0:4, legend=names(acc))
