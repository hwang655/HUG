library(HUG)
library(qlcMatrix)
library(huge)
library(HurdleNormal)
library(PLNet)
library(doParallel)
library(R.matlab)

cl = makeCluster(10)
registerDoParallel(cl)
getDoParWorkers()
set.seed(313813)
total = 20
n = 300
p = 100
time.PLN = time.glasso = time.glasso.npn = time.poisson = time.mcdavid = time.dep.pois = time.dep.hurd = rep(0,total)
for (nn in 1:total){
  
  # Dependent Sample Generating
  dat = rMVP(mu=rep(2, p), size=n, g.type='hub', factor.re=.5)                        # HPLN setting
  #dat = rMVH(coef.hurd=c(0, .5), mu=rep(3, p), size=n, g.type='hub', factor.re=.5)  # HHLN setting
  g.benchmark = dat$graph
  Y = dat$Y
  logY = log(as.matrix(Y)+1)
  ptm = proc.time()
  npnY = huge.npn(Y)
  npntime = proc.time()-ptm
  
  # mcdavid_c = fitHurdle(conditionalCenter(cbind(logY,apply(logY,1,function(x) sum(x!=0)))) , lambda.min.ratio=.4, parallel=T)
  ptm = proc.time()
  PLN = admmPLN(Y,lambda=0.5,rho = 0.5)$A
  diag(PLN) = 0
  cuts = quantile(abs(PLN[upper.tri(PLN)]),(seq(0,49)+50*2)/(50*3+1))
  PLN = lapply(cuts[50:1],function(v) Matrix(ifelse(abs(PLN)>v,1,0)))
  temp = proc.time()-ptm
  time.PLN[nn] = time.PLN+temp[2]+temp[1]

  # Graph Estimation
  ptm = proc.time()
  glasso.log = hugeGraph(logY, ngraph=50,lambda.min.ratio=1e-3, 'glasso', verbose=F)
  temp = proc.time()-ptm
  time.glasso[nn] = temp[2]+temp[1]
  ptm = proc.time()
  glasso.npn = hugeGraph(npnY, ngraph=50, lambda.min.ratio=1e-3,'glasso', verbose=F)
  temp = proc.time()-ptm
  time.glasso.npn[nn] = temp[2]+temp[1]+npntime[2]+npntime[1]
  ptm = proc.time()
  mb.pois = GLMGraph(Y, logY, Z=NULL, 50, 'poisson', lambda.min.ratio=1e-3)
  temp = proc.time()-ptm
  time.poisson[nn] = temp[2]+temp[1]
  ptm = proc.time()
  mcdavid = fitHurdle(logY, lambda.min.ratio=.4, parallel=T)
  temp = proc.time()-ptm
  time.mcdavid[nn] = temp[2]+temp[1]

  # Dependent graphical model estimation (our methods)
  B.ini = matrix(0, p+1, p)
  for (j in 1:p) B.ini[-(j+1), j] = as.vector(mb.pois$coef.aic[,j])
  ptm = proc.time()
  mb.dep.pois = MPoisGraph(Y, logY, Z=NULL,B.ini, lambda.min.ratio=1e-3)
  temp = proc.time()-ptm
  time.dep.pois[nn] =  temp[2]+temp[1]
  ptm = proc.time()
  mb.dep.hurd = MHurdGraph(Y, logY, Z=NULL,B.ini, lambda.min.ratio=1e-3)
  temp = proc.time()-ptm
  time.dep.hurd[nn] = temp[2]+temp[1]
  # sqrsolution = readMat(paste0("result",p,nn,".mat"))$result
  # for (i in 1:50){
  #   sqrsolution[[i]] = sqrsolution[[i]][[1]]
  # }
}
