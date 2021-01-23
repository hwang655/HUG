library(HUG)
library(qlcMatrix)
library(huge)
library(doParallel)


library(HurdleNormal)
library(PLNet)

pdata <- function (mu, size, g.type, factor.re = 1,nozero) 
{
  n = size
  p = nozero
  r.type = "block"
  Sigma = diag(n)
  ix2 = round(cumsum(c(0.1, 0.2, 0.3, 0.4)) * n)
  ix1 = c(1, ix2[-length(ix2)] + 1)
  for (k in 1:length(ix1)) Sigma[ix1[k]:ix2[k], ix1[k]:ix2[k]] = 0.8
  diag(Sigma) = 1
  Sigma = Sigma/factor.re
  Omega = solve(Sigma)
  stopifnot(all(eigen(Omega)$values > 0))
  Psi = rPrec(g.type, p)$icov
  Z.f = mvtnorm::rmvnorm(n, mu[1:p], as.matrix(solve(Psi)))
  p = length(mu)
  Z.r = mvtnorm::rmvnorm(p, numeric(n), Sigma)
  Z.f = cbind(Z.f,matrix(0,n,p-nozero))
  X = Z.f + t(Z.r)
  Y = matrix(rpois(n * p, exp(c(X))), n, p)
  Psi = cbind(rbind(Psi,matrix(0,(p-nozero),nozero)),matrix(0,p,(p-nozero)))
  diag(Psi) = 0
  list(Y = Y, X = X, Z = list(fixed = Z.f, random = t(Z.r)), 
       graph = (Psi != 0))
}

hdata <- function (coef.hurd, mu, size, g.type, factor.re = 1,nozero) 
{
  n = size
  p = nozero
  r.type = "block"
  Sigma = diag(n)
  ix2 = round(cumsum(c(0.1, 0.2, 0.3, 0.4)) * n)
  ix1 = c(1, ix2[-length(ix2)] + 1)
  for (k in 1:length(ix1)) Sigma[ix1[k]:ix2[k], ix1[k]:ix2[k]] = 0.8
  diag(Sigma) = 1
  Sigma = Sigma/factor.re
  Omega = solve(Sigma)
  stopifnot(all(eigen(Omega)$values > 0))
  Psi = rPrec(g.type, p)$icov
  Omega = solve(Psi)
  Z.f = mvtnorm::rmvnorm(n, mu[1:p], Omega)
  p = length(mu)
  Z.r = mvtnorm::rmvnorm(p, numeric(n), Sigma)
  Z.f = cbind(Z.f,matrix(0,n,p-nozero))
  X = c(Z.f + t(Z.r))
  Y = matrix((1 + rpois(n * p, exp(X))) * rbinom(n * p, 1, 
                                                 1/(1 + exp(-coef.hurd[1] - coef.hurd[2] * X))), n, p)
  Psi = cbind(rbind(Psi,matrix(0,(p-nozero),nozero)),matrix(0,p,(p-nozero)))
  diag(Psi) = 0
  list(Y = Y, X = matrix(X, n, p), Z = list(fixed = Z.f, random = t(Z.r)), 
       graph = (Psi != 0))
}

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)
getDoParWorkers()
total = 5
# dim.seq=100
dim.seq = c(100,200,400)
# dim.seq = c(50,100,150,200)

acc.cum <- vector(mode = "list", length = length(dim.seq))
for (p.index in  1:length(dim.seq)){
  set.seed(31381301)
  print(paste0("p.index:",p.index))
  for (nn in 1:total){
    print(paste0("rep_time:",nn))
    n = 100
    p = dim.seq[p.index]
    nozero=50
    # Dependent Sample Generating
    dat = pdata(mu=rep(2, p), size=n, g.type='hub', factor.re=.5,nozero = 50)                        # HPLN setting
    # dat = hdata(coef.hurd=c (0, .5), mu=rep(3, p), size=n, g.type='hub', factor.re=.5,50)  # HHLN setting
    g.benchmark = dat$graph
    Y = dat$Y
    logY = log(as.matrix(Y)+1)
    npnY = huge.npn(Y)
  
    # Graph Estimation
    glasso.log = hugeGraph(logY, 50, 'glasso', verbose=F)
    glasso.npn = hugeGraph(npnY, 50, 'glasso', verbose=F)
    mcdavid = fitHurdle(logY, nlambda=50, parallel=T)
    PLN = admmPLN(Y,lambda=0.5,rho = 0.5)$A
    diag(PLN) = 0
    cuts = quantile(abs(PLN[upper.tri(PLN)]),(seq(1,50)+50*9)/(50*10))
    PLN = lapply(cuts[50:1],function(v) Matrix(ifelse(abs(PLN)>v,1,0)))
    mcdavid_c = fitHurdle(conditionalCenter(logY) , lambda.min.ratio=.4, parallel=T)
    
    mb.pois = GLMGraph(Y, logY, Z=NULL, 50, 'poisson', lambda.min.ratio=3e-4)
    
    # Dependent graphical model estimation (our methods)
    B.ini = matrix(0, p+1, p)
    for (j in 1:p) B.ini[-(j+1), j] = as.vector(mb.pois$coef.aic[,j])
    mb.dep.pois = MPoisGraph(Y, logY,Z=NULL, B.ini, lambda.min.ratio=5e-3)
    mb.dep.hurd = MHurdGraph(Y, logY,Z=NULL, B.ini, lambda.min.ratio=5e-3)

    Sigma = diag(n)
    ix2 = round(cumsum(c(0.1, 0.2, 0.3, 0.4)) * n)
    ix1 = c(1, ix2[-length(ix2)] + 1)
    for (k in 1:length(ix1)) Sigma[ix1[k]:ix2[k], ix1[k]:ix2[k]] = 0.8
    diag(Sigma) = 1
    Sigma = Sigma/0.5
    Omega = solve(Sigma)
    mb.true.dep.pois = MPoisGraph(Y, logY,NULL, Omega=Omega, B.ini, lambda.min.ratio=5e-3)
    mb.true.dep.hurd = MHurdGraph(Y, logY,NULL, Omega=Omega, B.ini, lambda.min.ratio=5e-3)
    ests = list(glasso=lapply(glasso.log$path,function(x) x[1:50,1:50]), glasso.npn=lapply(glasso.npn$path,function(x) x[1:50,1:50])
                ,poisson=lapply(mb.pois$graphs,function(x) x[1:50,1:50])
                ,hurdle=lapply(mcdavid$adjMat,function(x) x[1:50,1:50])
                ,hurdle_center=lapply(mcdavid_c$adjMat,function(x) x[1:50,1:50])
                ,PLN=lapply(PLN,function(x) x[1:50,1:50])
                # ,sqrgraph = sqrsolution
                ,dep.poisson=lapply(mb.dep.pois$graphs,function(x) x[1:50,1:50]), dep.hurdle=lapply(mb.dep.hurd$graphs,function(x) x[1:50,1:50])
                ,dep.poisson.true = lapply(mb.true.dep.pois$graphs,function(x) x[1:50,1:50]),dep.hurdle.true=lapply(mb.true.dep.hurd$graphs,function(x) x[1:50,1:50])
                )
    temp = lapply(ests, GraphAcc, graph.t=g.benchmark[1:nozero,1:nozero])
    acc = if (nn == 1) temp else mapply("+", acc, temp , SIMPLIFY = FALSE)
  }
  acc.cum[[p.index]] = lapply(acc,function(x) x/total)
}
save(acc.cum,file = paste0("Output/compare_goldstandard_100-400_poisson"))

pdf(file = paste0("compare_goldstandard_100-400_poisson_2linetype",".pdf"),width=18,height=8)
par(mfcol=c(1,3),oma=c(7, 0, 0, 0), xpd=F)
for (i in 1:length(dim.seq)){
  plot(c(0, 100), c(0, 50), type='n',xlab='# of edges we found',ylab='# of edges we confirmed',cex.lab=1.5,cex.axis=1)
  rindex = rep(list(rep(1,1)),length(acc.cum[[i]]))
  for (j in 1:length(acc.cum[[i]])){
    last = acc.cum[[i]][[j]][rindex[[j]][1],2]
    index = 2
    while(index<80){
      while(((acc.cum[[i]][[j]][index,4]-last)<10)&(index<nrow(acc.cum[[i]][[j]])))
        index=index+1
      rindex[[j]] = c(rindex[[j]],index)
      last = acc.cum[[i]][[j]][index,2]
      index=index+1
    }
  }
  for(k in 1:length(acc.cum[[i]])){
    with(acc.cum[[i]][[k]][seq(1,nrow(acc.cum[[i]][[k]]),3),],lines(DF/2,TP/2,col=rainbow(length(acc.cum[[i]]))[k],lty=(k%%2+1)))
    with(acc.cum[[i]][[k]][seq(1,nrow(acc.cum[[i]][[k]]),3),],points(DF/2,TP/2,col=rainbow(length(acc.cum[[i]]))[k],pch=k-1,cex=1.5))
  }
  title(main = paste0("dim=",dim.seq[i],", type=poisson"),cex.main=2 )
}
par(xpd=NA)
legend(x=-220, y=-10, col=rainbow(length(acc.cum[[i]])), pch=0:(length(acc.cum[[i]])-1), legend=names(acc.cum[[i]]), ncol=5, box.lty=3,cex=2,lty = c((1:length(acc.cum[[i]]))%%2+1))
par(xpd=F)
dev.off()