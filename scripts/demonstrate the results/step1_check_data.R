
library(data.table)
library(ggplot2)
library(ggpubr)
# https://briatte.github.io/ggnet/
library(GGally)
library(network)
library(sna)

theme_set(theme_classic())

args = commandArgs(TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use '5278' as default.\n")
  ind = "5278"
}else{
  eval(parse(text=args[[1]]))
}

ind

ASD.list = c("5278", "5531")
control.list = c("4341", "5958")

if(ind %in% ASD.list){
  status = "ASD"
}else if(ind %in% control.list){
  status = "control"
}else{
  stop("unknown ind.\n")
}

# ------------------------------------------------------------------------
# read in SFARI gene list
# ------------------------------------------------------------------------

sf = fread("SFARI-Gene_genes_05-19-2020release_06-06-2020export.csv")
dim(sf)
sf[1:2,]

table(sf$status)
length(unique(sf$`gene-symbol`))
table(sf$`gene-score`, useNA='ifany')
table(sf$`syndromic`, useNA='ifany')

table(sf$`gene-score`, sf$`syndromic`, useNA='ifany')

# ------------------------------------------------------------------------
# read in graph estimation results
# ------------------------------------------------------------------------

path1 = 'edge found/catalysis precedes all ebic edges'
path1 = sprintf("%s/%s_velmeshev_%s", path1, ind, status)

fnms = list.files(path=path1)
dhurdF = grep("^dep_hurd_opt_graphs", fnms, value=TRUE)
dpoisF = grep("^dep_pois_opt_graphs", fnms, value=TRUE)
genesF = grep("genes1500_genes_Velmeshev_", fnms, value=TRUE)

dhurdF
dpoisF
genesF

load(file.path(path1, genesF))
ls()

load(file.path(path1, dhurdF))
ls()

load(file.path(path1, dpoisF))
ls()

length(genes.now)
genes.now[1:5]

dim(opt_graph_dep.hurd)
dim(opt_graph_dep.pois)

opt_graph_dep.hurd[1:2,]
opt_graph_dep.pois[1:2,]

fun1 <- function(v){paste(sort(v), collapse=":")}
pairs.dhurd = apply(opt_graph_dep.hurd, 1, fun1)
pairs.dpois = apply(opt_graph_dep.pois, 1, fun1)

length(unique(pairs.dhurd))
length(unique(pairs.dpois))

length(intersect(pairs.dhurd,pairs.dpois))

# ------------------------------------------------------------------------
# add SFARI information to graph model
# ------------------------------------------------------------------------

d_dhurd = sort(table(c(opt_graph_dep.hurd)), decreasing=TRUE)
d_dhurd[1:10]

length(d_dhurd)
summary(as.numeric(d_dhurd))
d_dhurd = data.frame(gene=names(d_dhurd), d=as.numeric(d_dhurd), 
                     stringsAsFactors = FALSE)

dim(d_dhurd)
d_dhurd[1:2,]

table(sf$`gene-symbol` %in% genes.now)
table(sf$`gene-symbol` %in% d_dhurd$gene)

w2kp = which(sf$`gene-symbol` %in% genes.now)
table(sf$`gene-score`[w2kp], sf$`syndromic`[w2kp], useNA='ifany')

w2kp = which(sf$`gene-symbol` %in% d_dhurd$gene)
table(sf$`gene-score`[w2kp], sf$`syndromic`[w2kp], useNA='ifany')

sf$score = sf$`gene-score`
sf$score[which(sf$syndromic==1)] = 3
table(sf$score, sf$`gene-score`, useNA = 'ifany')
table(sf$score, sf$`syndromic`,  useNA = 'ifany')

d_dhurd$SFARIscore = rep(0, nrow(d_dhurd))
wmat = which(d_dhurd$gene %in% sf$`gene-symbol`)
d_dhurd$SFARIscore[wmat] = sf$score[match(d_dhurd$gene[wmat],sf$`gene-symbol`)]
table(d_dhurd$SFARIscore)
dim(d_dhurd)
d_dhurd[1:2,]

summary(d_dhurd$d[d_dhurd$SFARIscore > 0])
summary(d_dhurd$d[d_dhurd$SFARIscore == 0])

wilcox.test(d_dhurd$d ~ d_dhurd$SFARIscore > 0)
wilcox.test(d_dhurd$d ~ d_dhurd$SFARIscore >=3)

# ------------------------------------------------------------------------
# get adjacency matrix for the subset of genes with any connection
# ------------------------------------------------------------------------

A = diag(nrow=nrow(d_dhurd), ncol=nrow(d_dhurd))
rownames(A) = colnames(A) = d_dhurd$gene

for(i in 1:nrow(opt_graph_dep.hurd)){
  gene1 = opt_graph_dep.hurd[i,1]
  gene2 = opt_graph_dep.hurd[i,2]
  A[gene1, gene2] = 1
  A[gene2, gene1] = 1
}
sum(c(A))

d_dhurd$SFARIscore_numeric = d_dhurd$SFARIscore
d_dhurd$SFARIscore = as.factor(d_dhurd$SFARIscore)

dA = solve(diag(rowSums(A))) %*% A
d_dhurd$SFARIscore_1s = dA %*% d_dhurd$SFARIscore_numeric
d_dhurd$SFARIscore_2s = dA %*% dA %*% d_dhurd$SFARIscore_numeric

cbp1 <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")

p1 = ggviolin(d_dhurd, x = "SFARIscore", y = c("SFARIscore_1s", "SFARIscore_2s"),
         combine = TRUE, color = "SFARIscore", palette = "jco",
         ylab = "SFARIscore after nb aggregation") + 
  scale_colour_manual(values=cbp1)

p3 = ggplot(d_dhurd, aes(x=SFARIscore_1s, y=SFARIscore_2s, 
                         col=SFARIscore, shape=SFARIscore)) + 
  geom_point() + scale_colour_manual(values=cbp1) +
  scale_shape_manual(values=c(16, 17, 15, 1)) + 
  theme(legend.position="top", legend.box = "horizontal")

pdf(sprintf("figures/dhurd_%s_catalysis-precedes_SFAR.pdf", ind), 
    width=6.5, height=3.5)
print(ggarrange(p1, p3, ncol=2, nrow=1))
dev.off()

sc1 = d_dhurd$SFARIscore_1s >=1 | d_dhurd$SFARIscore_2s >=1
dd1 = d_dhurd[which(sc1 & d_dhurd$SFARIscore_numeric == 0),]

sc2 = d_dhurd$SFARIscore_1s >=2 | d_dhurd$SFARIscore_2s >=2
dd2 = d_dhurd[which(sc2 & d_dhurd$SFARIscore_numeric > 0),]

dim(dd1)
dim(dd2)
dd1
dd2

dd1$SFARIscore_2s = round(dd1$SFARIscore_2s,4)
fwrite(dd1[,-4], sep="\t", 
       file=sprintf("output/dhurd_%s_catalysis-precedes_SFAR.tsv", ind))

# ------------------------------------------------------------------------
# try to draw the network, 1 step
# ------------------------------------------------------------------------

for(i in 1:nrow(dd1)){
  dd1.A2 = A[,match(dd1$gene[i], d_dhurd$gene)]
  length(dd1.A2)
  table(dd1.A2)
  w2plot = which(dd1.A2 > 0)
  
  net = network(A[w2plot,w2plot], vertex.attr=list(d_dhurd$gene[w2plot]), 
                vertex.attrnames=list("gene"), 
                directed = FALSE, 
                matrix.type="adjacency")
  set.vertex.attribute(net,"SFARIscore", as.character(d_dhurd$SFARIscore)[w2plot])
  net
  
  cbp1 = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
  
  gn = ggnet2(net, color = "SFARIscore", shape = "SFARIscore", size = 3, alpha = 0.75, 
              edge.alpha = 0.5, label = TRUE, vjust = -1, label.size=3, 
              palette = c("0" = "#E69F00", "1" = "#56B4E9", 
                          "2" = "#009E73", "3" = "#CC79A7"),
              shape.palette=c("0" = 16, "1" = 17, "2" = 15, "3" = 1)) +
    xlim(-0.12, 1.12) + ylim(-0.12,1.12) + 
    theme( axis.text.x = element_blank(), axis.line=element_blank(), 
           axis.text.y = element_blank(), axis.ticks=element_blank())
  
  pdf(sprintf("figures/ggnet_1_%s_%s_%s.pdf", ind, status, dd1$gene[i]), 
      width=3.75, height=3)
  print(gn)
  dev.off()
}

# ------------------------------------------------------------------------
# try to draw the network, 2 step
# ------------------------------------------------------------------------

A2 = (A %*% A > 0)
dim(A2)
A2[1:5,1:5]

for(i in 1:nrow(dd1)){
  dd1.A2 = A2[,match(dd1$gene[i], d_dhurd$gene)]
  length(dd1.A2)
  table(dd1.A2)
  w2plot = which(dd1.A2)
  
  net = network(A[w2plot,w2plot], vertex.attr=list(d_dhurd$gene[w2plot]), 
                vertex.attrnames=list("gene"), 
                directed = FALSE, 
                matrix.type="adjacency")
  set.vertex.attribute(net,"SFARIscore", as.character(d_dhurd$SFARIscore)[w2plot])
  net
  
  cbp1 = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
  
  gn = ggnet2(net, color = "SFARIscore", size = 3, alpha = 0.75, 
         edge.alpha = 0.5, label = TRUE, vjust = -1, label.size=3, 
         palette = c("0" = "#E69F00", "1" = "#56B4E9", 
                     "2" = "#009E73", "3" = "#CC79A7")) +
    xlim(-0.12, 1.12) + ylim(-0.12,1.12) + 
    theme( axis.text.x = element_blank(), axis.line=element_blank(), 
           axis.text.y = element_blank(), axis.ticks=element_blank())
  
  if(length(w2plot) < 10){
    width = 3.75; height=3
  }else if(length(w2plot) < 20){
    width = 5; height=4
  }else{
    width = 10; height=8
  }
  
  pdf(sprintf("figures/ggnet_2_%s_%s_%s.pdf", ind, status, dd1$gene[i]), 
      width=width, height=height)
  print(gn)
  dev.off()
  
}

# ------------------------------------------------------------------------
# try to draw the network for SFAR
# ------------------------------------------------------------------------

w.SFAR = which(d_dhurd$SFARIscore_numeric>0)

net = network(A[w.SFAR,w.SFAR], vertex.attr=list(d_dhurd$gene[w.SFAR]), 
              vertex.attrnames=list("gene"), 
              directed = FALSE, 
              matrix.type="adjacency")
set.vertex.attribute(net,"SFARIscore", as.character(d_dhurd$SFARIscore)[w.SFAR])
net

cbp1 = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

gn = ggnet2(net, color = "SFARIscore", size = 3, alpha = 0.75, 
            edge.alpha = 0.5, label = TRUE, vjust = -1, label.size=2.5, 
            palette = c("0" = "#E69F00", "1" = "#56B4E9", 
                        "2" = "#009E73", "3" = "#CC79A7")) +
  xlim(-0.12, 1.12) + ylim(-0.12,1.12) + 
  theme( axis.text.x = element_blank(), axis.line=element_blank(), 
         axis.text.y = element_blank(), axis.ticks=element_blank())

pdf(sprintf("figures/ggnet_SFAR_genes_%s_%s.pdf", ind, status), 
    width=6, height=4.75)
print(gn)
dev.off()

gc()

sessionInfo()
q(save="no")

