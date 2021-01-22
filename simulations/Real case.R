library(doSNOW)
library(HUG)
library(qlcMatrix)
library(huge)
# library(HurdleNormal)
library(doParallel)
library(MASS)
cores<-parallel::detectCores()
cl<-makeSOCKcluster(12)
registerDoSNOW(cl)
name.people=c("4341","5278", "5531","5958","gierahn")
rel.index=7
for(rel.index in c(7,1)){
  for (index.people in 5:5){
    #ImportData
    load('Data/relation_12.rData')
    load(paste0('Data/pathwayCommons12_',rel[rel.index],'.rData'))
    names.data=c(paste0('Velmeshev/genes_Velmeshev_',name.people[index.people],'_5cluster'))
    name = paste0('Velmeshev_',name.people[index.people])
    if (name.people[index.people]=="gierahn"){
      load(paste0('Data/Gierahn2017_raw.rData'))
      data.save=dat.gierahn
      name = "Gierahn"
    }
    else{
      load(paste0('Data/',names.data,'.rData'))
    }
    genes=intersect(colnames(conn.gene),colnames(data.save))
    nz.v=colMeans(data.save[,genes]>=0.5)
    genes=intersect(colnames(conn.gene),colnames(data.save))
    nz.v=colMeans(data.save[,genes[!(startsWith(genes,'RPL')|startsWith(genes,'RPS'))]]>=0.5)
    genes.now = genes[order(nz.v,decreasing = TRUE)[1:1500]]
    read_depth=rowSums(data.save)
    dat.Velmeshev=as.matrix(data.save)
    g.benchmark=conn.gene[genes.now,genes.now]
    p=length(genes.now)
    g.benchmark=g.benchmark|t(g.benchmark)
    cat(p,'genes,',sum(g.benchmark)/2,'edges\n')
    k=1
    Y=dat.Velmeshev[,genes.now]
    logY=log((Y+1)/read_depth)
    npnY=huge.npn(Y/read_depth)
    cat(round(mean(Y==0)*100),'%ofobservationsare0\n')
    # mcdavid = fitHurdle(logY, nlambda=50, parallel=T)
    glasso.log=hugeGraph(logY,100,'glasso',lambda.min.ratio=0.5)
    #GraphEstimation
    glasso.npn=hugeGraph(npnY,100,'glasso',lambda.min.ratio=0.05)
    mb.pois=GLMGraph(Y,npnY,Z=log(as.matrix(read_depth,,1)),100,'poisson',lambda.min.ratio=1e-6)
    B.ini=matrix(0,p+2,p)
    for(j in 1:(p))
      B.ini[-(j+2),j]=as.vector(mb.pois$coef.aic[,j])
    R=scale(npnY);R[is.na(R)]=0
    huge.cells=huge(t(R[,]),method='glasso',lambda.min.ratio = 0.7)
    O.cells=as.matrix(huge.cells$icov[[length(huge.cells$icov)]])
    mb.dep.pois=MPoisGraph(Y,npnY,Z=log(as.matrix(read_depth,,1)),B.ini,Omega=O.cells,lambda.max=3,lambda.min.ratio=1e-4,nlambda=100)
    mb.dep.hurd=MHurdGraph(Y,npnY,Z=log(as.matrix(read_depth,,1)),B.ini,Omega=O.cells,lambda.max=3,lambda.min.ratio=1e-4,nlambda=100)
    ests=list(glasso=glasso.log$path,glasso.npn=glasso.npn$path,poisson=mb.pois$graphs,
              dep.poisson=mb.dep.pois$graphs,dep.hurdle=mb.dep.hurd$graphs)
    acc=lapply(ests,GraphAcc,graph.t=g.benchmark)
    pdf(paste0('Eval_Case_',name,'_',rel.index,'.pdf'),height=5,width=5)
    if (rel.index == 7){
      plot(c(0,500),c(0,30),type='n',xlab='# of edges we found',ylab='# of edges we confirmed',main=name,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    }
    if (rel.index == 1){
      plot(c(0,5000),c(0,30),type='n',xlab='# of edges we found',ylab='# of edges we confirmed',main=name,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    }
    rindex = rep(list(rep(1,1)),length(acc))
    for (i in 1:length(acc)){
      last = acc[[i]][rindex[[i]][1],2]
      index = 2
      while(index<50){
        while(((acc[[i]][index,2]-last)<50)&(index<50))
          index=index+1
        rindex[[i]] = c(rindex[[i]],index)
        last = acc[[i]][index,2]
        index=index+1
      }
    }
    for(i in 1:length(acc)){
      with(acc[[i]][rindex[[i]],],lines(DF/2,TP/2,col=i,lty=1))
      with(acc[[i]][rindex[[i]],],points(DF/2,TP/2,col=i,pch=i-1))
    }
    legend("topleft",col=1:length(acc),lty=1,pch=0:(length(acc)-1),legend=names(acc),horiz=F,box.lty=3)
    dev.off()
  }
}