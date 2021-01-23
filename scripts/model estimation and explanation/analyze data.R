# bsub -n 12 -M 48 -R span[hosts=1] -Ip R
# bsub -o MB.out -q day -n 12 -M 48 -R span[hosts=1] Rscript SkewedCount_Test.r
.libPaths(c(.libPaths(), '../Rlibs'))
library(foreach)
library(doParallel)
library(Matrix)
library(qlcMatrix)
library(fields)
library(mvtnorm)
library(glmnet)
# library(mpath)
library(sgd)
library(grpreg)
library(huge)
library(HUG)
# library(PLNet)
# library(lvnet)

cl = makeCluster(10)
registerDoParallel(cl)
getDoParWorkers()

# rm(list=ls())
source('Functions_Simulate.r')
source('Functions_GLMM.r')
source('Functions_GE.r')
source('Functions_IntMatrixGraph.r')

# Benchmark Graph
load(paste0('Data/pathwayCommons12_catalysis-precedes.rData'))
p=1500

# scRNA-seq
load('Data/Gierahn2017_raw.rData')  # Seq-Well
genes=intersect(colnames(conn.gene),colnames(dat.gierahn))
nz.v=colMeans(as.matrix(dat.gierahn)[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
zeroprop=0.1
genes.now = genes[order(nz.v,decreasing = TRUE)[1:1500]]
load(paste0('Data/Velmeshev/genes_Velmeshev_5387_5cluster.rData'))
{
dat.velmeshev.5387 = data.save
genes=intersect(colnames(conn.gene),colnames(data.save))
nz.v=colMeans(data.save[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
zeroprop=0.1
genes.now.5387 = genes[order(nz.v,decreasing = TRUE)[1:1500]]
}
load(paste0('Data/Velmeshev/genes_Velmeshev_5939_5cluster.rData'))
{
  dat.velmeshev.5939 = data.save
  genes=intersect(colnames(conn.gene),colnames(data.save))
  nz.v=colMeans(data.save[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
  genes.now.5939 = genes[order(nz.v,decreasing = TRUE)[1:1500]]
  }
load(paste0('Data/Velmeshev/genes_Velmeshev_5531_5cluster.rData'))
{
  dat.velmeshev.5531 = data.save
  genes=intersect(colnames(conn.gene),colnames(data.save))
  nz.v=colMeans(data.save[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
  genes.now.5531 = genes[order(nz.v,decreasing = TRUE)[1:1500]]
  }


Y.all = list(Gierahn=dat.gierahn[,genes.now], 
             Velmeshev.5387 = dat.velmeshev.5387[,genes.now.5387], 
             Velmeshev.5939 = dat.velmeshev.5939[,genes.now.5939], 
             Velmeshev.5531 = dat.velmeshev.5531[,genes.now.5531])
names.data = names(Y.all)
cat('Zero proportions: ', round(sapply(Y.all, function(x) mean(x != 0)), 3), '\n')
logY.all = list(log(as.matrix(Y.all$Gierahn)+1),log(Y.all$Velmeshev.5387+1),log(Y.all$Velmeshev.5939+1),log(Y.all$Velmeshev.5531+1))

# EDA
pzero.actual = sapply(Y.all, function(x) colMeans(x == 0))
params.PLN = lapply(Y.all, function(x) matrix(unlist(apply(x, 2, FitPLN)), 2, p))
pzero.PLN = sapply(params.PLN, function(pas) apply(pas, 2, function(pa) dPLN(0, pa[1], pa[2])))

cor.matrix.test <- function(S, n) {
  df = n - 2
  r2 = S[upper.tri(S)]^2
  Fstat = r2 * df / (1 - r2)
  1 - pf(Fstat, 1, df)
}

R.all = sapply(c(1,2,3,4), function(k) {R=t(scale(logY.all[[k]])); R[is.na(R)]=0; R})
names(R.all) = names.data[c(1,2,3,4)]
S.all = sapply(R.all, cor)
pval.all = sapply(S.all, function(S) cor.matrix.test(S, p))

pdf('Figure1.pdf', height=5, width=8); par(mfrow=c(2, 3))
for (k in c(1,2)) {
  hist(colMeans(Y.all[[k]] > 0), main=names.data[k], xlab='Expression Proportion in All Cells')
  plot(pzero.PLN[,k], pzero.actual[,k], xlim=c(0, .7), ylim=c(0, .7), 
       xlab='Predicted Zero-Proportion', ylab='Actual Zero-Proportion', 
       main=names.data[k])
  abline(0, 1, col=2)  
  hist(pval.all[[k]], main=names(R.all)[k], xlab='P-value of Pearson Correlation Test')

}
dev.off()

pdf('Figure2.pdf', height=5, width=8); par(mfrow=c(2, 3))
for (k in c(3,4)) {
  hist(colMeans(Y.all[[k]] > 0), main=names.data[k], xlab='Expression Proportion in All Cells')
  plot(pzero.PLN[,k], pzero.actual[,k], xlim=c(0, .7), ylim=c(0, .7), 
       xlab='Predicted Zero-Proportion', ylab='Actual Zero-Proportion', 
       main=names.data[k])
  abline(0, 1, col=2)  
  hist(pval.all[[k]], main=names(R.all)[k], xlab='P-value of Pearson Correlation Test')
  
}
dev.off()