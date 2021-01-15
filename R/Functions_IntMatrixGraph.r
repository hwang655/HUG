# Functions to Estimate Graphs based on Matrix-Variate Non-Negative Integer Data
# (Hurdle-Poisson-logNormal model for each variate)

# ----------------------- Poisson-logNormal model inference -----------------------
# Parameter estimation
FitPLN = function(Y, name.data='') {
    # poilog::poilogMLE doesn't work well
    stopifnot(all(Y >= 0))
    var.Y = var(Y)
    mean.Y = mean(Y)
    if (var.Y <= mean.Y) {
        cat(name.data,
            'Sample mean (', round(mean.Y, 4), ') >= sample variance (', round(var.Y, 4), ') \n')
        return(list(mu=log(mean.Y), sigma=0))
    }
    sigma2 = log(1 + (var.Y-mean.Y) / mean.Y^2)
    mu = log(mean.Y) - sigma2 / 2
    list(mu=mu, sigma=sqrt(sigma2))
}
# Retrieve underlying Normal data
TransPLN = function(Y, name.data='') {
    if (is.vector(Y)) return(TransPLN.vec(Y))
    foreach(k=1:ncol(Y), .combine=cbind) %do% TransPLN.vec(Y[,k], paste(name.data, 'Col', k, ':'))
}
TransPLN.vec = function(Y, name.data='') {
    pa = FitPLN(Y, name.data)
    if (pa$sigma == 0) return(ifelse(Y == 0, log(Y+0.1), log(Y)))
    unique.Y = unique(Y)
    trans.Y = try({
        out = numeric(length(Y))
        for (y in unique.Y) {
            out[Y == y] = uniroot(function(x) y-exp(x)-(x-pa$mu)/pa$sigma^2, 
                                  pa$mu + pa$sigma*c(-10, 10))$root
        }
        out
    })
    if (class(trans.Y) == 'try-error') return(log(0.1+Y))
    trans.Y
}
# Goodness of fit test
PLN.test = function(Y, min.count=10, name.data='') {
    stopifnot(all(Y >= 0))
    if (var(Y) == 0) {
        return(1)
    }
    n = length(Y)
    M = 10
    pa = FitPLN(Y, name.data)
    probs = if (pa$sigma==0) dpois(0:M, exp(pa$mu)) else poilog::dpoilog(0:M, pa$mu, pa$sigma)
    probs[length(probs)] = probs[length(probs)] + 1-sum(probs)
    expected = n * probs
    while (length(expected) > 2 && expected[length(expected)] < min.count) {
        expected = c(expected[1:(length(expected)-2)], 
                     expected[length(expected)-1] + expected[length(expected)])
    }
    observed = sapply(1:length(expected), function(k) sum(Y==k-1))
    observed[length(observed)] = n - sum(observed[-length(observed)])
    out = sum((observed-expected)^2 / expected)
    1 - pchisq(out, length(observed)-1-1-(pa$sigma>0))
}
dPLN = function(x, mu, sigma) {
    if (sigma==0) dpois(x, exp(mu)) else poilog::dpoilog(x, mu, sigma)
}

# Estimate the graphs of the Hurdle-Mixture model
# Y is the data matrix (n x p)
# Y2 represent extra explanatory covariates
HUG.Graphs = function(Y, Y2 = NULL, nlambda = 100, maxiter = 1, output.RE = T, Omega = NULL,
                      criterion.RE = 'aic') {
    n    = nrow(Y)
    p    = ncol(Y)
    Z    = (Y != 0) + 0
    stopifnot(all(apply(Z, 2, var) > 0))
    X    = log(Y + 0.1)
    p2   = ifelse(is.null(Y2), 1, 1 + ncol(Y2))
    res  = list(Gammas = list(), Bs = list(), Omegas = list(), U = list(), V = list())
    Omega.now = Omega
    maxiter   = maxiter + is.null(Omega)
    
    for (iter in 1:maxiter) {
        if (is.null(Omega)) {
            output.all = foreach(j = 1:p, .packages = c('glmnet')) %dopar% {
                cat('iter = ', iter, ', j = ', j, '\n')
                A.j   = which(Y[,j] != 0)
                fm1.j = cv.glmnet(X[,-j], Z[,j], family = 'binomial', nlambda = nlambda)
                fm2.j = cv.glmnet(X[A.j, -j], Y[A.j, j] - 1, family = 'poisson', nlambda = nlambda)
                fm.j  = list(gamma = coef(fm1.j), beta = coef(fm2.j))
                fm.j
            }
        } else {
            output.all = foreach(j = 1:p, .packages = c('HUG')) %dopar% {
                source('Functions_Hurdle.r')
                source('Functions_GLMM.r')
                cat('iter = ', iter, ', j = ', j, '\n')
                fm.j = HurdleMix.Reg(X[,-j], Z = Y2, Y[,j], Omega = Omega.now, 
                                     penalty.type = 'weighted', nlambda = nlambda,
                                     method.refit = 'mle', criterion1 = 'ebic')
                fm.j
            }
        }
        # (p2+p) x p sparse matrices
        Gamma.now = foreach(j = 1:p, .combine = cbind) %do% {
            gamma.j = numeric(p2 + p)
            gamma.j[-(p2+j)] = output.all[[j]]$gamma
            Matrix(gamma.j, sparse = T)
        }
        B.now = foreach(j = 1:p, .combine = cbind) %do% {
            b.j = numeric(p)
            b.j[-(p2+j)] = output.all[[j]]$beta
            Matrix(b.j, sparse = T)
        }
        # n x p dense matrices or NULL
        U.now = foreach(j = 1:p, .combine = cbind) %do% {
            output.all[[j]]$u # logistic
        }
        V.now = foreach(j = 1:p, .combine = cbind) %do% {
            output.all[[j]]$v # Poisson
        }
        Omega.now = GraphRE.Hurdle(Y, Y2 = Y2, B = B.now, Gamma = Gamma.now, algo = 'glasso', 
                                   Omega.ini = Omega.now, U.ini = U.now, V.ini = V.now, 
                                   nlambda = 10, lambda.min.ratio = 0.05)
        Omega.now = with(Omega.now, Omega[[which.min(eval[,criterion.RE])]])
        
        res$Gammas[[iter]] = Gamma.now
        res$Bs[[iter]]     = B.now
        res$Omegas[[iter]] = Omega.now
        if (output.RE) {
            res$U[[iter]] = U.now
            res$V[[iter]] = V.now
        }
    }
    return(res)
}

# Estimate the graphs of the Poisson-Mixture model
# Y is the data matrix (n x p)
# Y2 represent extra explanatory covariates
PoisMix.Graphs = function(Y, Y2 = NULL, nlambda = 100, nthres = 20, 
                      gamma.ebic = 0.5, maxiter = 2, eps = 1e-4, 
                      maxiter.IRLS = 100, verbose = F,
                      output.U = T, Omega = NULL) {
    n = nrow(Y)
    p = ncol(Y)
    X = log(Y + 0.1)
    p2 = ifelse(is.null(Y2), 1, 1 + ncol(Y2))
    res = list(Bs = list(), Omegas = list(), Us = list())
    if (is.null(Omega)) {
        Omega.now = sparseMatrix(i = 1:n, j = 1:n, x = 1)
    } else {
        Omega.now = Omega
    }
    for (iter in 1:maxiter) {
        FindU_PoisMix = FindU_PoisMix
        output.all = foreach(j = 1:p) %dopar% {
            library(Matrix)
            library(irlba)
            library(HUG)
            source('Functions_Poisson.r')
            source('Functions_gglasso.r')
            
            cat('iter = ', iter, ', j = ', j, '\n')
            fm.j  = list(beta = numeric(p+p2-1))
            if (var(Z[,j]) != 0) {
                fm.j = PoisMix.Reg(X[,-j], Z = Y2, Y[,j], Omega = Omega.now, 
                               mixture.model = T, nlambda = nlambda, method.fit = 'thres',
                               nthres = nthres, criterion1 = 'bic', criterion2 = 'ebic', 
                               gamma.ebic = gamma.ebic, maxiter.IRLS = maxiter.IRLS)
            }
            fm.j
        }
        # a (p+p2)-by-p sparse matrix
        B.now = foreach(j = 1:p, .combine = cbind) %do% {
            b.j = numeric(p+p2)
            b.j[-(p2 + j)] = output.all[[j]]$beta
            Matrix(b.j, sparse = T)
        }
        # an n-by-p dense matrix
        U.now = foreach(j = 1:p, .combine = cbind) %do% {
            output.all[[j]]$u
        }
        Omega.now = GraphRE.Poisson(Y, Y2, B = B.now, Omega.ini = Omega, nlambda = 50)
        # Omega.now = GraphGaussian(t(U.now))
        res$Bs[[iter]] = B.now
        res$Gammas[[iter]] = Gamma.now
        res$Omegas[[iter]] = Omega.now
        if (output.U) {
            res$Us[[iter]] = U.now
        }
    }
    return(res)
}

# Graph estimation with binary nodes
GraphLogistic = function(Y, Y2 = NULL, nlambda = 100, gamma.ebic = 0.5) {
    n = nrow(Y)
    p = ncol(Y)
    Z = (Y != 0)
    X = ifelse(Y == 0, 0, log(Y + 0.1))
    if (is.null(Y2)) {
        Y2 = matrix(1, n, 1)
    } else {
        stopifnot(nrow(Y2) == n)
        Y2 = cbind(1, Y2)
    }
    p2 = ncol(Y2)
    
    Eval.GLMM = Eval.GLMM
    coef.mat = foreach(j = 1:p, .combine = cbind) %dopar% {
        b.j  = Matrix(0, p2 + p, 1, sparse = T)
        if (var(Z[,j]) != 0) {
            X.j  = cbind(Y2, X[, -j])
            fm.j = glmnet(X.j[,-1], Z[,j], family = 'binomial', nlambda = 100,
                          penalty.factor = c(numeric(p2-1), rep(1, p-1)))
            B.j  = coef(fm.j)
            EBIC.j = Eval.GLMM(Z[,j], X.j, family = 'binomial', B = B.j, 
                               p2 = p2, gamma.ebic = gamma.ebic)[,'ebic']
            EBIC.j[is.na(EBIC.j)] = max(EBIC.j, na.rm = T)
            b.j[-(p2+j)] = B.j[, which.min(EBIC.j)]
        }
        b.j
    }
    return(coef.mat) # (p2 + p) x p matrix
}

# L=huge::huge.generator(n=200,d=200,graph="hub"); X=L$data; nlambda=20; criterion='ric'
# Graphical LASSO with tuning parameter selection
GraphGaussian = function(X, nlambda = 20, criterion = 'ric') {
    fm = huge::huge(X, nlambda = nlambda, method = 'glasso')
    sel = huge::huge.select(fm, criterion = criterion)
    return(sel$opt.icov)
}








