library(doSNOW)
#library(xtable)
library(HUG)
library(qlcMatrix)
library(huge)
# library(HurdleNormal)
library(doParallel)
library(MASS)
#EvaluateGLMmodels
process.m <- function(matrix){
  g = matrix(0,gene.num,gene.num)
  g1 = matrix[3:(gene.num+1),]
  g[upper.tri(g)] = g1[upper.tri(g1)]
  g1 = rbind(0,g1)
  g[lower.tri(g)] = g1[lower.tri(g1)]
  return(pmax(abs(g),t(abs(g))))
}
  name.people=c("4341","5278", "5531","5958")
  glasso.level=0.7
  gene.num = 1500
  
  rel.index=1
  load('Data/relation_12.rData')
  load(paste0('Data/pathwayCommons12_',rel[rel.index],'.rData'))
for(index.people in 1:4){
  #ImportData
  names.data=c(paste0('genes_Velmeshev_',name.people[index.people],'_5cluster'))
  load(paste0('Data/Velmeshev/',names.data,'.rData'))
  genes=intersect(colnames(conn.gene),colnames(data.save))
  nz.v=colMeans(data.save[,genes]>=0.5)
  genes=intersect(colnames(conn.gene),colnames(data.save))
  nz.v=colMeans(data.save[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
  genes.now = genes[order(nz.v,decreasing = TRUE)[1:gene.num]]
  k=1
  load(file=paste0('Output/dep_hurd_',names.data[k],'_',rel[rel.index],'_glasso',glasso.level,'.rData'))
  load(file=paste0('Output/dep_pois_',names.data[k],'_',rel[rel.index],'_glasso',glasso.level,'.rData'))
  load(file=paste0('Output/pois_',names.data[k],'_',rel[rel.index],'.rData'))
  load(file=paste0('Output/glasso_',names.data[k],'_',rel[rel.index],'.rData'))
  load(file=paste0('Output/glasso_npn_',names.data[k],'_',rel[rel.index],'.rData'))
  g = process.m(mb.dep.hurd$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(500)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_dep.hurd = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_dep.hurd, file = paste0('dep_hurd_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_dep.hurd),'edges.rData'))

  
  g = process.m(mb.dep.pois$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(500)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_dep.pois = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_dep.pois, file = paste0('dep_pois_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_dep.pois),'edges.rData'))

  
  g = process.m(mb.pois$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(500)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_pois = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_pois, file = paste0('pois_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_pois),'edges.rData'))

}
  {
  index.people=5
  names.data=('gierahn')
  load(paste0('Data/Gierahn2017_raw.rData'))
  genes=intersect(colnames(conn.gene),colnames(dat.gierahn))
  nz.v=colMeans(dat.gierahn[,genes]>=0.5)
  genes=intersect(colnames(conn.gene),colnames(dat.gierahn))
  nz.v=colMeans(dat.gierahn[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
  genes.now = genes[order(nz.v,decreasing = TRUE)[1:gene.num]]
  load(file=paste0('Output/dep_hurd_',names.data,'_',rel[rel.index],'_glasso',glasso.level,'.rData'))
  load(file=paste0('Output/dep_pois_',names.data,'_',rel[rel.index],'_glasso',glasso.level,'.rData'))
  load(file=paste0('Output/pois_',names.data,'_',rel[rel.index],'.rData'))
  load(file=paste0('Output/glasso_',names.data,'_',rel[rel.index],'.rData'))
  load(file=paste0('Output/glasso_npn_',names.data,'_',rel[rel.index],'.rData'))
  g = process.m(mb.dep.hurd$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(100)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_dep.hurd = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_dep.hurd, file = paste0('dep_hurd_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_dep.hurd),'edges.rData'))
  
  
  g = process.m(mb.dep.pois$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(100)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_dep.pois = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_dep.pois, file = paste0('dep_pois_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_dep.pois),'edges.rData'))
  
  
  g = process.m(mb.pois$coef.opt)
  g[upper.tri(g)] = 0
  res <- order(g,decreasing = TRUE)[seq_len(100)]
  pos <- arrayInd(res, dim(g), useNames = TRUE)
  # pos = which(g!=0,arr.ind = T)
  opt_graph_pois = cbind(genes.now[pos[,1]],genes.now[pos[,2]])
  save(opt_graph_pois, file = paste0('pois_opt_graphs_',name.people[index.people],'_',rel[rel.index],'_glasso',glasso.level,'_',nrow(opt_graph_pois),'edges.rData'))
  }
  