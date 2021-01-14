# Graphical models via GLASSO
hugeGraph = function(Y, ngraph=50, method='glasso', ...) {
    S = corSparse(Y); S[is.na(S)] = 0; diag(S) = 1
    huge(S, nlambda=ngraph, method=method, ...)
}

# Dependent Hurdle graph estimation via LSA
MHurdGraph = function(Y.r, Y.p, Z=NULL, B.ini, Omega=NULL, coef.hurd=NULL, 
                      nlambda=50, lambda.max=NULL, lambda.min.ratio=0.01) {
    n = nrow(Y.r)
    if (is.null(Z)) {
        Z = matrix(1, n, 0)
    }
    Z = as.matrix(Z)
    pZ = ncol(Z)
    p = ncol(Y.r)
    stopifnot(nrow(B.ini) == p+pZ+1)
    logY.r = log(as.matrix(Y.r)+1)
    Z.r = 0 + (Y.r > 0)
    XB = as.matrix(cbind(1, Z, Y.p) %*% B.ini)
    cvec = apply(logY.r - cbind(1, Z, Y.p) %*% B.ini, 2, var)
    if (is.null(Omega)) {
        R = scale(logY.r); R[is.na(R)] = 0
        Omega = with(huge(cor(t(R)), method='glasso'), as.matrix(icov[[length(icov)]]))
    }
    v.ini = log(Y.r+1) - XB
    Ghat1 = XB + v.ini
    Yhat1 = exp(Ghat1)
    Vhat1 = Yhat1
    Ytrans1 = Ghat1 - (Yhat1 + 1 - Y.r) / Vhat1
    if (is.null(coef.hurd)) {
        coef.hurd = matrix(coef(glm(c(Z.r) ~ c(XB), family=binomial)), p, 2, byrow=T)
    } else if (is.vector(coef.hurd)) {
        coef.hurd = matrix(coef.hurd, p, 2, byrow=T)
    }
    Ghat0 = matrix(coef.hurd[,1], n, p, byrow=T) + Ghat1*matrix(coef.hurd[,2], n, p, byrow=T)
    Yhat0 = 1 / (1+exp(-Ghat0))
    Vhat0 = Yhat0 * (1-Yhat0) * matrix(coef.hurd[,2], n, p, byrow=T)^2
    Ytrans0 = Ghat1 - (Yhat0 - (Y.r>0)) / Vhat0 * matrix(coef.hurd[,2], n, p, byrow=T)
    
    Y.tilde = as.matrix(chol(Omega) %*% logY)
    Resid.r = ResidualMatrix(Y.tilde, Z, family='gaussian', intercept=T)
    S = t(scale(Y.tilde)) %*% Resid.r / n * 2
    diag(S) = 0
    lambdaseq = exp(log(ifelse(is.null(lambda.max), max(abs(S), na.rm=T), lambda.max)) +
        seq(0, log(lambda.min.ratio), length.out=nlambda))
    rm(v.ini, Ghat1, Ghat0, Y.tilde, Resid.r, S); gc()
    pb <- txtProgressBar(min=1, max=p, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    Betas = foreach(j=1:p, .packages=c('Matrix', 'glmnet', 'foreach'), .options.snow=opts) %dopar% {
        time.old = Sys.time()
        out = matrix(FALSE, p, nlambda)
        b.opt = numeric(pZ+p)
        b.aic = b.opt
        if (var(Y.r[,j]) != 0) {
            v1.j = Vhat1[,j] * (Y.r[,j]>0)
            v0.j = Vhat0[,j]
            K.j = solve(Omega/cvec[j] + diag(v1.j+v0.j))
            DD1.j = diag(v1.j) - matrix(v1.j, n, n) * K.j * matrix(v1.j, n, n, byrow=T)
            DD0.j = diag(v0.j) - matrix(v0.j, n, n) * K.j * matrix(v0.j, n, n, byrow=T)
            K.j = matrix(v0.j, n, n) * K.j * matrix(v1.j, n, n, byrow=T)
            V.j = DD1.j + DD0.j - K.j - t(K.j)
            ix.keep = which(colMeans(abs(V.j)) > 1e-10) # avoid singularity due to very small v1
            Ytrans.j = solve(V.j[ix.keep, ix.keep], 
                             ((DD0.j-t(K.j))%*%Ytrans0[,j] + (DD1.j-K.j)%*%Ytrans1[,j])[ix.keep])
            rm(DD0.j, DD1.j, K.j); gc()
            Q.j = chol(V.j[ix.keep, ix.keep])
            X.j = Q.j %*% cbind(1, Z, Y.p)[ix.keep, -(j+1+pZ)]
            Y.j = Q.j %*% Ytrans.j
            B = try(coef(glmnet(X.j, Y.j, family='gaussian', lambda=lambdaseq, 
                                penalty.factor=c(numeric(1+pZ), rep(1, p-1)), intercept=F))[-1,])
            # since coef() produces an NA intercept even if intercept=F
            if (class(B) != 'try-error') {
                eval.B = Eval.LSA(Y.j, X.j, B=B, p2=1+pZ, gamma.ebic='auto')
                b.opt = as.vector(B[,which.min(eval.B$ebic)])
                b.aic = as.vector(B[,which.min(eval.B$aic)])
                B = cbind(B, B[,rep(ncol(B), nlambda-ncol(B))]) # in case ncol(B) < nlambda
                out[-j,] = as.matrix(B[-(1:(1+pZ)),]) != 0
            }
        }
        list(B=out, b.opt=b.opt, b.aic=b.aic, time=difftime(Sys.time(), time.old, units='secs'))
    }
    close(pb)
    list(lambda=lambdaseq, 
         graphs=foreach(i=1:nlambda, .packages='Matrix') %dopar% {
             Matrix(foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$B[,i])
         }, 
         coef.opt=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.opt, 
         coef.aic=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.aic, 
         coef.hurd=coef.hurd, 
         time=sapply(Betas, `[[`, 'time'))
}

# Dependent Poisson graph estimation via LSA
MPoisGraph = function(Y.r, Y.p, Z=NULL, B.ini, Omega=NULL, nlambda=50, lambda.max=NULL, lambda.min.ratio=0.01) {
    n = nrow(Y.r)
    if (is.null(Z)) {
        Z = matrix(1, n, 0)
    }
    Z = as.matrix(Z)
    pZ = ncol(Z)
    p = ncol(Y.r)
    stopifnot(nrow(B.ini) == p+pZ+1)
    logY.r = log(as.matrix(Y.r)+1)
    XB = as.matrix(cbind(1, Z, Y.p) %*% B.ini)
    cvec = apply(logY.r - XB, 2, var)
    if (is.null(Omega)) {
        R = scale(logY.r); R[is.na(R)] = 0
        Omega = with(huge(cor(t(R)), method='glasso'), as.matrix(icov[[length(icov)]]))
    }
    v.ini = log(Y.r+1) - XB
    g.hat = XB + v.ini
    Y.hat = exp(XB + v.ini)
    Y.trans = as.matrix(Y.r) / Y.hat - 1 + g.hat
    Sigma = Matrix::solve(Omega)
    
    Y.tilde = as.matrix(chol(Omega) %*% logY)
    Resid.r = ResidualMatrix(Y.tilde, Z, family='gaussian', intercept=T)
    S = t(scale(Y.tilde)) %*% Resid.r / n * 2
    diag(S) = 0
    lambdaseq = exp(log(ifelse(is.null(lambda.max), max(abs(S), na.rm=T), lambda.max)) +
        seq(0, log(lambda.min.ratio), length.out=nlambda))
    rm(v.ini, g.hat, Y.tilde, Resid.r, S); gc()
    
    pb <- txtProgressBar(min=1, max=p, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    Betas = foreach(j=1:p, .packages=c('Matrix', 'glmnet', 'foreach'), .options.snow=opts) %dopar% {
        time.old = Sys.time()
        out = matrix(FALSE, p, nlambda)
        b.opt = c(log(mean(Y.r[,j])), numeric(p+pZ-1))
        b.aic = b.opt
        if (var(Y.r[,j]) != 0) {
            Q.j = chol(solve(cvec[j]*Sigma + diag(1/Y.hat[,j])))
            X.j = Q.j %*% cbind(1, Z, Y.p[,-j])
            Y.j = Q.j %*% Y.trans[,j]
            rm(Q.j); gc()
            B = try(coef(glmnet(X.j, Y.j, family='gaussian', lambda=lambdaseq, 
                                penalty.factor=c(numeric(1+pZ), rep(1, p-1)), intercept=F))[-1,])
            # since coef() produces an NA intercept even if intercept=F
            if (class(B) != 'try-error') {
                eval.B = Eval.LSA(Y.j, X.j, B=B, p2=1+pZ, gamma.ebic='auto')
                b.opt = as.vector(B[,which.min(eval.B$ebic)])
                b.aic = as.vector(B[,which.min(eval.B$aic)])
                B = cbind(B, B[,rep(ncol(B), nlambda-ncol(B))]) # in case ncol(B) < nlambda
                out[-j,] = as.matrix(B[-(1:(1+pZ)),]) != 0
            }
        }
        list(B=out, b.opt=b.opt, b.aic=b.aic, time=difftime(Sys.time(), time.old, units='secs'))
    }
    close(pb)
    list(lambda=lambdaseq, 
         graphs=foreach(i=1:nlambda, .packages='Matrix') %dopar% {
             Matrix(foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$B[,i])
         }, 
         coef.opt=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.opt, 
         coef.aic=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.aic, 
         time=sapply(Betas, `[[`, 'time'))
}


# Graphical models via nodewise GLM
GLMGraph = function(Y.r, Y.p, Z=NULL, nlambda=50, family='gaussian', intercept=T, lambda.max=NULL,
                    lambda.min.ratio=0.01, weight.zero=1, hetero.penalty=F) {
    n = nrow(Y.r)
    if (is.null(Z)) {
        Z = matrix(1, n, 0)
    }
    Z = as.matrix(Z)
    pZ = ncol(Z)
    p = ncol(Y.r)
    if (hetero.penalty) {
        penalty.factor = 1 / ifelse(colMeans(Y.p), colMeans(Y.p), 1)^2
    } else {
        penalty.factor = rep(1, p)
    }
    Resid.r = ResidualMatrix(Y.r, Z, family=family, intercept=intercept)
    lambdaseq = NULL
    if (family %in% c('poisson', 'gaussian')) {
        S = t(scale(Y.p, center=intercept)) %*% Resid.r / nrow(Y.p) * ifelse(family=='gaussian', 2, 1)
        diag(S) = 0
        lambdaseq = exp(log(ifelse(is.null(lambda.max), max(abs(S), na.rm=T), lambda.max)) +
            seq(0, log(lambda.min.ratio), length.out=nlambda))
    } else {
        stop('Unsupported family')
    }
    
    pb <- txtProgressBar(min=1, max=p, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    Betas = foreach(j=1:p, .packages=c('Matrix', 'glmnet', 'foreach'), .options.snow=opts) %dopar% {
        out = matrix(FALSE, p, nlambda)
        b.opt = c(mean(Y.r[,j]), numeric(p+pZ-1))
        b.aic = b.opt
        if (var(Y.r[,j]) != 0) {
            weights.j = ifelse(Y.r[,j]==0, weight.zero, 1)
            B = try(coef(glmnet(cbind(Z, Y.p[,-j]), Y.r[,j], family=family, lambda=lambdaseq, 
                                penalty.factor=c(numeric(pZ), penalty.factor[-j]), intercept=intercept, 
                                weights=weights.j)))
            if (class(B) != 'try-error') {
                eval.B = Eval.GLM(Y.r[,j], cbind(1, Z, Y.p[,-j]), family=family, B=B, 
                                  p2=1+pZ, gamma.ebic='auto', weights=weights.j)
                b.opt = as.vector(B[,which.min(eval.B$ebic)])
                b.aic = as.vector(B[,which.min(eval.B$aic)])
                B = cbind(B, B[,rep(ncol(B), nlambda-ncol(B))]) # in case ncol(B) < nlambda
                out[-j,] = as.matrix(B[-(1:(1+pZ)),]) != 0
            }
        }
        list(B=out, b.opt=b.opt, b.aic=b.aic)
    }
    close(pb)
    list(lambda=lambdaseq, 
         graphs=foreach(i=1:nlambda, .packages='Matrix') %dopar% {
             Matrix(foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$B[,i])
         }, 
         coef.opt=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.opt, 
         coef.aic=foreach(j=1:p, .combine=cbind) %do% Betas[[j]]$b.aic)
}

# Nodewise residual matrix
ResidualMatrix = function(Y, Z, family='gaussian', intercept=T) {
    Y = as.matrix(Y)
    if (is.null(Z) || (ncol(Z)==0)) return(if (intercept) scale(Y, scale=F) else Y)
    Z = as.matrix(Z)
    foreach(j=1:ncol(Y), .combine=cbind) %dopar% {
        y.hat = predict(if (intercept) glm(Y[,j]~Z, family=family) else glm(Y[,j]~Z-1, family=family))
        if (family == 'gaussian') {
            Y[,j] - y.hat
        } else if (family == 'poisson') {
            Y[,j] - exp(y.hat)
        } else {
            stop('Unsupported family')
        }
    }
}

# Graphical models based on truncated correlation matrix
corGraph = function(Y, ngraph=50) {
    S = corSparse(Y); S[is.na(S)] = 0
    cuts = quantile(abs(S[upper.tri(S)]), (seq(0, ngraph-1)+ngraph*9)/(ngraph*10+1))
    return(lapply(cuts[ngraph:1], function(v) Matrix(ifelse(abs(S)>v, 1, 0))))
}

# Evaluate GLM models
# X contains an constant column; first row of B are the intercepts
Eval.GLM = function(y, X, offset=NULL, family=c('gaussian', 'binomial', 'poisson'), B, 
                    p2=1, gamma.ebic=0.5, weights=NULL) {
    family = match.arg(family)
    n = nrow(X)
    p = ncol(X) - p2
    B = as.matrix(B)
    K = ncol(B)
    stopifnot(nrow(B) == ncol(X))
    if (is.null(offset)) { offset = numeric(n) }
    if (gamma.ebic == 'auto') { gamma.ebic = 1-log(p)/log(n)/2 }
    if (is.null(weights)) { weights = rep(1, n) }

    p.nz = colSums(as.matrix(B[-(1:p2),]) != 0)
    logmu = as.matrix(X %*% B) + matrix(offset, n, K, byrow = F)
    mu = exp(logmu)
    if (family == 'binomial') {
        mean.all = mu / (1 + mu)
        var.all = mean.all * (1 - mean.all)
        loglik = as.vector(t(y*weights) %*% logmu) - colSums(log(1+mu)*weights)
    } else if (family == 'gaussian') {
        mean.all = logmu
        var.all = colSums((matrix(y, n, K) - logmu)^2 * weights) / sum(weights)
        loglik = -sum(weights)/2*log(var.all)
    } else if (family == 'poisson') {
        mean.all = mu
        var.all = mu
        loglik = as.vector(t(y*weights) %*% logmu) - colSums(mu*weights)
    } else {
        stop(family, ' family not supported yet.')
    }
    eval.all = matrix(NA, K, 4)
    eval.all[,1] = -loglik
    eval.all[,2] = 2*eval.all[,1] + p.nz*2
    eval.all[,3] = 2*eval.all[,1] + p.nz*log(n)
    eval.all[,4] = eval.all[,3] + 2*gamma.ebic * 
        (lgamma(p+1) - lgamma(p.nz+1) - lgamma(p-p.nz+1))
    colnames(eval.all) = c('loglik', 'aic', 'bic', 'ebic')
    rownames(eval.all) = colnames(B)
    
    return(data.frame(eval.all))
}

Eval.LSA = function(y, X, offset=NULL, B, p2=1, gamma.ebic=0.5, weights=NULL) {
    n = nrow(X)
    p = ncol(X) - p2
    B = as.matrix(B)
    K = ncol(B)
    stopifnot(nrow(B) == ncol(X))
    if (is.null(offset)) { offset = numeric(n) }
    if (gamma.ebic == 'auto') { gamma.ebic = 1-log(p)/log(n)/2 }
    if (is.null(weights)) { weights = rep(1, n) }
    
    p.nz = colSums(as.matrix(B[-(1:p2),]) != 0)
    logmu = as.matrix(X %*% B) + matrix(offset, n, K, byrow = F)
    mu = exp(logmu)
    mean.all = logmu
    loglik = -colSums((matrix(y, n, K) - logmu)^2 * weights) / 2
    eval.all = matrix(NA, K, 4)
    eval.all[,1] = -loglik
    eval.all[,2] = 2*eval.all[,1] + p.nz*2
    eval.all[,3] = 2*eval.all[,1] + p.nz*log(n)
    eval.all[,4] = eval.all[,3] + 2*gamma.ebic * 
        (lgamma(p+1) - lgamma(p.nz+1) - lgamma(p-p.nz+1))
    colnames(eval.all) = c('loglik', 'aic', 'bic', 'ebic')
    rownames(eval.all) = colnames(B)
    
    return(data.frame(eval.all))
}

