# SummaryGraphs: convert a bunch of symmetric matrices to graphs and provide a comparison
# AsSymmetric    : symmetrize matrix for nodewise estimation

GraphAcc = function(graphs, graph.t) {
    if (!is.list(graphs)) {
        graphs = list(graphs)
    }
    K      = length(graphs)
    p      = nrow(graph.t)
    N      = p * (p-1)
    diag(graph.t) = FALSE
    
    err = foreach(i = 1:K, .combine = c) %do% {
        graph.i = (graphs[[i]] != 0)
        diag(graph.i) = FALSE
        tmp = c(sum(graph.i & graph.t), sum(graph.i), sum(graph.t))
        c(tmp[1], tmp[2]-tmp[1], tmp[3]-tmp[1], tmp[2], tmp[3])
    }
    
    crit   = c("TP", "FP", "FN", "DF", "DF_T")
    err = matrix(err, ncol = length(crit), byrow = TRUE)
    colnames(err) = crit
    try(rownames(err) <- names(graphs))
    err = as.data.frame(err)
    err$Lift = with(err, N * TP / DF_T / DF)
    
    return(err)
}


# Symmetrize Nodewise Graphical Estimation
# Nullize conflicting entries
# Take the bigger entry in the output
AsSymmetric = function(Omega, edge.union = T, method = 'small') {
    p = nrow(Omega)
    Omega.new = Matrix(0, p, p)
    
    ix.nz = which(Omega != 0, arr.ind = T)
    ix1.nz = ix.nz[Omega[ix.nz[,c(2,1)]] == 0, ]
    ix2.nz = ix.nz[Omega[ix.nz[,c(2,1)]] != 0, ]
    nz1.Omega = Omega[ix1.nz]
    nz2.Omega = Omega[ix2.nz]
    nz2.OmegaT = Omega[ix2.nz[,c(2,1)]]
    
    if (method == 'small') {
        nz2.new = ifelse(nz2.Omega < nz2.OmegaT, nz2.Omega, nz2.OmegaT)
    } else if (method == 'large') {
        nz2.new = ifelse(nz2.Omega > nz2.OmegaT, nz2.Omega, nz2.OmegaT)
    } else {
        nz2.new = (nz2.Omega + nz2.OmegaT) / 2
    }
    Omega.new[ix2.nz] = nz2.new
    
    if (edge.union) {
        Omega.new[ix1.nz] = nz1.Omega
        Omega.new[ix1.nz[,c(2,1)]] = nz1.Omega
    }
    
    return(Omega.new)
}

# Generate UNDIRECTED Graphs from Symmetric Matrices
SummaryGraphs = function(icovs){
    eps = 1e-10
    if (!is.list(icovs)){
        if (attr(class(icovs), "package") == "Matrix") {
            icovs = as.matrix(icovs)
        }
        stopifnot(is.matrix(icovs), isSymmetric(icovs))
        icovs = list(icovs)
    }
    K = length(icovs)
    p = nrow(icovs[[1]])
    lower.p = which(lower.tri(diag(p)))
    
    graphs = list()
    sum.all = matrix(0, nrow = 1, ncol = K + 4)
    for (k in 1:K){
        graphs[[k]] = (abs(icovs[[k]]) > eps) + 0
        diag(graphs[[k]]) = 0
        sum.all[1, k] = sum(graphs[[k]][lower.p])
    }
    
    path.union = graphs[[1]]
    icov.nonn = (icovs[[1]] >= 0)
    icov.nonp = (icovs[[1]] <= 0)
    if (K > 1){
        for (k in 2:K){
            path.union = path.union + graphs[[k]]
            icov.nonn = icov.nonn + (icovs[[k]] >= 0)
            icov.nonp = icov.nonp + (icovs[[k]] <= 0)
        }
    }
    icov.diff = !((icov.nonn == K) | (icov.nonp == K))
    graph.union = (path.union > 0) + 0
    graph.common = (path.union > 1)
    graph.diff = icov.diff * path.union
    sum.all[1, K + 1] = sum(graph.union[lower.p])
    sum.all[1, K + 2] = sum(graph.common[lower.p])
    sum.all[1, K + 3] = sum(graph.diff[lower.p])
    
    A = p * (p-1) / 2
    sum.all[1, K + 4] = exp(log(sum.all[1, K+2]) - sum(log(sum.all[1, 1:K])) + (K-1)*log(A))
    
    colnames(sum.all) = c(paste0("df_", 1:K), "df_union", "df_common", "df_diff", "lift")
    sum.all = sum.all[1, c(K + (1:4), 1:K)]
    
    return(list(graph.union = graph.union, graph.common = graph.common,
        graph.diff = graph.diff, graph.each = graphs, sum.all = sum.all))
}

