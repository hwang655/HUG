# Generate a Poisson-logNormal sample with sample dependence
rMVP = function(mu, size, g.type, factor.re=1) {
    n = size
    p = length(mu)
    
    r.type = 'block'
    Sigma = diag(n)
    ix2 = round(cumsum(c(0.1, 0.2, 0.3, 0.4)) * n)
    ix1 = c(1, ix2[-length(ix2)]+1)
    for (k in 1:length(ix1)) Sigma[ix1[k]:ix2[k], ix1[k]:ix2[k]] = 0.8
    diag(Sigma) = 1
    Sigma = Sigma / factor.re
    Omega = solve(Sigma)
    stopifnot(all(eigen(Omega)$values > 0))
    Psi = rPrec(g.type, p)$icov
    
    Z.f = mvtnorm::rmvnorm(n, mu, as.matrix(solve(Psi)))
    Z.r = mvtnorm::rmvnorm(p, numeric(n), Sigma)
    X = Z.f + t(Z.r)
    Y = matrix(rpois(n*p, exp(c(X))), n, p)
    diag(Psi) = 0
    list(Y=Y, X=X, Z=list(fixed=Z.f, random=t(Z.r)), graph=(Psi != 0))
}

# Generate a Hurdle-logNormal sample with sample dependence
rMVH = function(coef.hurd, mu, size, g.type, factor.re=1) {
    n = size
    p = length(mu)
    
    r.type = 'block'
    Sigma = diag(n)
    ix2 = round(cumsum(c(0.1, 0.2, 0.3, 0.4)) * n)
    ix1 = c(1, ix2[-length(ix2)]+1)
    for (k in 1:length(ix1)) Sigma[ix1[k]:ix2[k], ix1[k]:ix2[k]] = 0.8
    diag(Sigma) = 1
    Sigma = Sigma / factor.re
    Omega = solve(Sigma)
    stopifnot(all(eigen(Omega)$values > 0))
    Psi = rPrec(g.type, p)$icov
    
    Z.f = mvtnorm::rmvnorm(n, mu, as.matrix(solve(Psi)))
    Z.r = mvtnorm::rmvnorm(p, numeric(n), Sigma)
    X = c(Z.f + t(Z.r))
    Y = matrix((1 + rpois(n*p, exp(X))) * 
                   rbinom(n*p, 1, 1 / (1+exp(-coef.hurd[1]-coef.hurd[2]*X))), n, p)
    diag(Psi) = 0
    list(Y=Y, X=matrix(X, n, p), Z=list(fixed=Z.f, random=t(Z.r)), graph=(Psi != 0))
}


rPrec = function(graph.type, d) {
    g = diag(d)
    if (graph.type == 'band') {
        # keep negative
        g[cbind(c(1:(d-1), 2:d), c(2:d, 1:(d-1)))] = -0.6
        g[cbind(c(1:(d-2), 3:d), c(3:d, 1:(d-2)))] = -0.3
    } else if (graph.type == 'hub') {
        size.hub = 10
        n.hub = ceiling(d / size.hub)
        for (k in 1:n.hub) {
            # keep negative
            g[(k-1)*size.hub+1, ((k-1)*size.hub+2):min(d, k*size.hub)] = -0.5
            g[((k-1)*size.hub+2):min(d, k*size.hub), (k-1)*size.hub+1] = -0.5
        }
    } else if (graph.type == 'power') {
        g = matrix(FALSE, d, d)
        g[2, 1] = TRUE
        for (j in 3:d) {
            sizes.now = Matrix::colSums(g[1:(j-1), 1:(j-1)]) +
                Matrix::rowSums(g[1:(j-1), 1:(j-1)])
            connect.j = sample(j-1, size=1, prob=sizes.now)
            g[j, connect.j] = TRUE
        }
        # either negative or positive
        g[g] = -runif(sum(g), 0.7, 0.8)
        g[upper.tri(g)] = t(g)[upper.tri(g)]
        diag(g) = 1
        eigval.min = min(0, min(eigen(g)$values))
        diag(g) = diag(g) + (abs(eigval.min) + 0.05)
    } else if (graph.type == 'ER') {
        prob = min(0.05, 5/d)
        ix.sel = sample(d*(d-1)/2, d*(d-1)/2*prob)
        ix.sel = which(upper.tri(diag(d)), arr.ind=T)[ix.sel,]
        # either negative or positive
        g[ix.sel] = -runif(nrow(ix.sel), 0.7, 0.8)
        g[ix.sel[,2:1]] = g[ix.sel]
    } else if (graph.type == 'dense') {
        # actually a latent factor model
        g = matrix(0.5, d, d)
        diag(g) = 1
    } else {
        stop('Invalid graph type!')
    }
    eigval.min = min(0, min(eigen(g)$values))
    diag(g) = diag(g) + (abs(eigval.min) + 0.05)
    list(icov=g, cov=solve(g), graph=Matrix(g-diag(diag(g)) != 0))
}

