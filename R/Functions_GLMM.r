# GLMM.Reg:       Fit a GLM mixed model with given Omega
# GLMM.Reg.IRLS:  Fit a GLM mixed model via IRLS and GMD
# ReFitGLMM:      Run a GLM (mixed) model refitting procedure
# FitGLMM:        Fit a GLM mixed model without penalty
# Eval.GLMM:      Evaluate a bunch of GLM (mixed) models

MHurdReg = function(X, y, b.ini, Omega, coef.hurd=NULL, nlambda=30, 
                    lambda.max=10, lambda.min.ratio=0.01) {
    n = nrow(X)
    p = ncol(X)
    S = t(scale(X, center=T)) %*% (y - mean(y)) / n
    maxlambda = min(max(abs(S), na.rm=T), lambda.max)
    lambdaseq = exp(log(maxlambda) + seq(0, log(lambda.min.ratio), length.out=nlambda))
    z = 0 + (y > 0)
    
    Ghat1 = as.vector(cbind(1, X) %*% b.ini)
    Yhat1 = exp(Ghat1)
    Vhat1 = Yhat1
    Ytrans1 = Ghat1 - (Yhat1 + 1 - y) / Vhat1
    if (is.null(coef.hurd)) {
        coef.hurd = coef(glm(z ~ Ghat1, family=binomial))
    }
    Ghat0 = coef.hurd[1] + Ghat1*coef.hurd[2]
    Yhat0 = 1 / (1+exp(-Ghat0))
    Vhat0 = Yhat0 * (1-Yhat0) * coef.hurd[2]^2
    Ytrans0 = Ghat1 - (Yhat0 - z) / Vhat0 * coef.hurd[2]
    
    B = matrix(0, p+1, nlambda)
    b.opt = numeric(p+1)
    b.aic = b.opt
    if (var(y) != 0) {
        v1 = Vhat1 * z
        v0 = Vhat0
        K = solve(Omega + diag(v1+v0))
        DD1 = diag(v1) - matrix(v1, n, n) * K * matrix(v1, n, n, byrow=T)
        DD0 = diag(v0) - matrix(v0, n, n) * K * matrix(v0, n, n, byrow=T)
        KK = matrix(v0, n, n) * K * matrix(v1, n, n, byrow=T)
        V = DD1 + DD0 - KK - t(KK)
        ix.keep = which(colMeans(abs(V)) > 1e-10) # avoid singularity due to very small v1
        Ytrans = solve(V[ix.keep, ix.keep], 
                       ((DD0-t(KK))%*%Ytrans0 + (DD1-KK)%*%Ytrans1)[ix.keep])
        rm(DD0, DD1, K, KK); gc()
        Q = chol(V[ix.keep, ix.keep])
        X = Q %*% cbind(1, X[ix.keep,])
        Y = Q %*% Ytrans
        B = try(coef(glmnet(X, Y, family='gaussian', lambda=lambdaseq, intercept=F))[-1,])
        if (class(B) != 'try-error') {
            eval.B = Eval.GLM(Y, X, family='gaussian', B=B, p2=1, gamma.ebic='auto')
            b.opt = as.vector(B[,which.min(eval.B$ebic)])
            b.aic = as.vector(B[,which.min(eval.B$aic)])
            B = cbind(B, B[,rep(ncol(B), nlambda-ncol(B))]) # in case ncol(B) < nlambda
        }
    }
    list(B=B, coef.opt=b.opt, coef.aic=b.aic, coef.hurd=coef.hurd, lambda=lambdaseq)
}

HurdReg = function(X, y, b.ini, coef.hurd=NULL, nlambda=30, 
                   lambda.max=10, lambda.min.ratio=0.01) {
    n = nrow(X)
    p = ncol(X)
    S = t(scale(X, center=T)) %*% (y - mean(y)) / n
    maxlambda = min(max(abs(S), na.rm=T), lambda.max)
    lambdaseq = exp(log(maxlambda) + seq(0, log(lambda.min.ratio), length.out=nlambda))
    z = 0 + (y > 0)
    
    Ghat1 = as.vector(cbind(1, X) %*% b.ini)
    Yhat1 = exp(Ghat1)
    Vhat1 = Yhat1
    Ytrans1 = Ghat1 - (Yhat1 + 1 - y) / Vhat1
    if (is.null(coef.hurd)) {
        coef.hurd = coef(glm(z ~ Ghat1, family=binomial))
    }
    Ghat0 = coef.hurd[1] + Ghat1*coef.hurd[2]
    Yhat0 = 1 / (1+exp(-Ghat0))
    Vhat0 = Yhat0 * (1-Yhat0) * coef.hurd[2]^2
    Ytrans0 = Ghat1 - (Yhat0 - z) / Vhat0 * coef.hurd[2]
    
    B = matrix(0, p+1, nlambda)
    b.opt = numeric(p+1)
    b.aic = b.opt
    if (var(y) != 0) {
        v1 = Vhat1 * z
        v0 = Vhat0
        X = rbind(matrix(v1, n, p+1)*cbind(1, X), matrix(v0, n, p+1)*cbind(1, X))
        Y = c(v1*Ytrans1, v0*Ytrans0)
        ix.keep = which(c(v1, v0) > 1e-10) # avoid biased (E)BIC
        X = X[ix.keep,]
        Y = Y[ix.keep]
        B = try(coef(glmnet(X, Y, family='gaussian', lambda=lambdaseq, intercept=F))[-1,])
        if (class(B) != 'try-error') {
            eval.B = Eval.GLM(Y, X, family='gaussian', B=B, p2=1, gamma.ebic='auto')
            b.opt = as.vector(B[,which.min(eval.B$ebic)])
            b.aic = as.vector(B[,which.min(eval.B$aic)])
            B = cbind(B, B[,rep(ncol(B), nlambda-ncol(B))]) # in case ncol(B) < nlambda
        }
    }
    list(B=B, coef.opt=b.opt, coef.aic=b.aic, coef.hurd=coef.hurd, lambda=lambdaseq)
}

# Q.cov is the Cholesky decomposition of the random effect covariance matrix
# i.e. t(Q.cov) %*% Q.cov = Sigma
MPoisReg = function(X, y, Q.cov, lambdaseq=NULL, nlambda=30, lambda.max=100, lambda.min.ratio=0.01) {
    n = nrow(X)
    p = ncol(X)
    S = t(scale(X, center=T)) %*% (y - mean(y)) / n
    if (is.null(lambdaseq)) {
        maxlambda = min(max(abs(S), na.rm=T), lambda.max)
        lambdaseq = exp(log(maxlambda) + seq(log(lambda.min.ratio), 0, length.out=nlambda))
    }
    df.tmp = as.data.frame(cbind(y, X))
    names(df.tmp) = c('y', paste0('X', 1:ncol(X)))
    tQ = t(Q.cov)
    colnames(tQ) = paste0('Q', 1:nrow(Q.cov))
    
    B = matrix(0, p+1, nlambda)
    eval.all = list(aic=rep(NA, nlambda), bic=rep(NA, nlambda), ebic=rep(NA, nlambda))
    phi.all = rep(NA, nlambda)
    sd.all = rep(NA, nlambda)
    control.list = list()
    for (k in 1:length(lambdaseq)) {
        fm.k = glmmLasso(as.formula(paste('y ~', paste(names(df.tmp)[-1], collapse='+'))), 
                         rnd=tQ, data=df.tmp, lambda=lambdaseq[k], 
                         family=poisson(), control=control.list)
        if (class(fm.k) != 'try-error') {
            control.list = list(start=c(fm.k$coefficients, fm.k$ranef), q.start = fm.k$StdDev)
            eval.all$aic[k] = fm.k$aic
            eval.all$bic[k] = fm.k$bic
            B[,k] = fm.k$coefficients
            sd.all[k] = fm.k$StdDev
            phi.all[k] = fm.k$phi
        }
    }
    # glmmLasso(fix, rnd=list(f=~1), ...) is equivalent to glmmLasso(fix, rnd=I) 
    # ... where f=factor(1:p) and I is an identity matrix
    p.nz = colSums(B[-1,] != 0)
    # eval.all$ebic = eval.all$bic + 2*gamma.ebic*(lgamma(p+1)-lgamma(p.nz+1)-lgamma(p-p.nz+1))
    list(B=B, lambda=lambdaseq, phi=phi.all, sd=sd.all, eval=eval.all, 
         coef.opt=B[,which.min(eval.all$bic)], coef.aic=B[,which.min(eval.all$aic)])
}

# # Fit mixed-effect Poisson regression via PQL and proximal gradient descent (ISTA)
# MixedPoisReg = function(X, y, Omega, cval=NULL, lambda=NULL, nlambda=10, learning_rate=1e-4) {
#     n = nrow(X)
#     p = ncol(X)
#     if (is.null(lambda)) {
#         lambda = numeric(nlambda)
#     }
#     nlambda = length(lambda)
#     X = cbind(1, as.matrix(X))
#     y = as.matrix(y)
#     beta_mat = matrix(0, p+1, nlambda)
#     c_vec = numeric(nlambda)
#     beta_now = matrix(0, p+1, 1)
#     v_now = matrix(0, n, 1)
#     
#     maxiter = 100
#     eps = 1e-5
#     for (ilam in 1:nlambda) {
#         for (iter in 1:maxiter) {
#             beta_old = beta_now
#             v_old = v_now
#             resid_now = exp(X %*% beta_now + v_now) - y
#             d_beta = t(X) %*% resid_now / n
#             beta_now = beta_now - learning_rate*d_beta
#             beta_now[-1] = SoftThres(beta_now[-1], learning_rate*lambda[ilam])
#             resid_now = exp(X %*% beta_now + v_now) - y
#             Ov = Omega %*% v_now
#             if (is.null(cval)) {
#                 factor.v = if (any(v_now != 0)) 1/as.numeric(t(v_now)%*%Ov) else 0
#             } else {
#                 factor.v = 1 / n / cval
#             }
#             d_v = resid_now / n + Ov * factor.v
#             v_now = v_now - learning_rate*d_v
#             if ((mean(abs((beta_now-beta_old))) < eps) && (mean(abs((v_now-v_old))) < eps)) {
#                 break
#             }
#         }
#         beta_mat[,ilam] = beta_now
#         c_vec[ilam] = as.numeric(t(v_now)%*%Ov) / n
#     }
#     return(list(beta=beta_mat, cval=c_vec))
# }

SoftThres = function(x, thres) {
    ix.null = (abs(x) <= thres)
    x[ix.null] = 0
    x[!ix.null] = x[!ix.null] - sign(x[!ix.null]) * thres
    return(x)
}



# GLM mixed model fitting/selection/refitting with given Omega
# algo: 'FISTA' - solve the original problem by gradient methods
#       'IRLS'  - solve the LA approximated problem (log|Omega + D| ignored completely) by IRLS
#                 or equivalently, treating U as missing values instead of random effects
# accu.grad: works only when algo = 'FISTA', i.e. solve the original problem
#            TRUE  - use the most accurate gradients computed by LA
#            FALSE - use the gradients of LA approximated problem 
#                    (currently ignoring log|Omega + D| completely, similar to algo = 'IRLS3')
# eps: convergence criterion for FISTA/IRLS iterations
GLMM.Reg = function(X, Z = NULL, y, offset = NULL, extend.X = F, Omega = NULL, 
                    family = c('gaussian', 'binomial', 'poisson'), 
                    algo = 'IRLS', accu.grad = F, PQL = F, nlambda = 100, 
                    lambda = NULL, lambda.max.ratio = 1, 
                    lambda.min.ratio = ifelse(ncol(X) > nrow(X), 0.1, 0.05),
                    method.refit = 'none', nthres = 20, 
                    criterion1 = 'ebic', criterion2 = 'ebic', gamma.ebic = 0.5, 
                    maxiter = ifelse(algo == 'IRLS', 5e2, 1e3), eps = 1e-3, verbose = F,
                    time.limit = Inf, sparse = F, ...) {
    stopifnot(is.matrix(X), !any(is.na(X)))
    family     = match.arg(family)
    y  = drop(y)
    n  = length(y)
    pX = ncol(X)
    if (extend.X) {
        stopifnot(pX %% 2 == 0, algo == 'IRLS')
        pX.0 = pX / 2
        colnames(X) = c(paste0('X', 1:pX.0), paste0('I', 1:pX.0))
    } else {
        pX.0 = pX
        colnames(X) = paste0('X', 1:pX)
    }
    if (max(abs(colMeans(X[,1:pX.0]))) > 1e-10) {
        cat('Conditional centering X...\n')
        CC.X = ConditionalCenter(X[,1:pX.0])
        X[,1:pX.0] = CC.X$X
    }
    sd.X = apply(X, 2, sd)
    X = X / matrix(sd.X, n, pX, byrow = T)
    X[, sd.X == 0] = 0
    sd.X[sd.X == 0] = 1
    if (is.null(Z)) {
        Z = matrix(1, n, 1)
        colnames(Z) = 'Constant'
    } else {
        colnames(Z) = paste0('Z', 1:ncol(Z))
        Z = cbind(Constant = 1, Z)
    }
    pZ      = ncol(Z)
    sd.Z = apply(Z, 2, sd); sd.Z[1] = 1
    Z = Z / matrix(sd.Z, n, pZ, byrow = T)
    Z[, sd.Z == 0] = 0
    
    p       = pX + pZ
    include = 1:pZ
    ZX      = cbind(Z, X)
    if (!is.null(Omega)) {
        is.diag.Omega = isDiagonal(Omega)
        sparse  = sparse || is.diag.Omega
        if (sparse) {
            Omega = as(Omega, 'dgCMatrix')
        } else {
            Omega = as.matrix(Omega)
        }
    }
    timing  = numeric(6)
    names(timing) = c('GMD', 'Newton', 'WholePath', 'Refit', 'Prep1', 'Prep2')
    ybar    = mean(y)
    
    # initialization and determine lambda range
    beta.ini = numeric(p)
    beta.ini[include] = coef(glm(y ~ ZX[,include] - 1, family = family, offset = offset))
    fm.ini = GLMM.Reg.IRLS(ZX, y, offset, F, pZ, Omega, family = family, PQL = PQL, 
                           lambda = 1e10, include = include, b0 = beta.ini, 
                           maxiter = maxiter, eps = eps, verbose = F)
    timing   = timing + fm.ini$timing
    beta.ini = as.vector(fm.ini$B)
    u.ini    = as.vector(fm.ini$U)
    if (is.null(lambda)) {
        if (extend.X) {
            lambdamax = max(sqrt(fm.ini$derv[(1:pX.0)+pZ]^2 + fm.ini$derv[(1:pX.0)+pZ+pX.0]^2))
        } else {
            lambdamax = max(abs(fm.ini$derv))
        }
        lambdamax = lambdamax * lambda.max.ratio
        lambda = exp(seq(log(lambdamax), log(lambdamax)+log(lambda.min.ratio), length.out=nlambda))
    } else {
        lambda = sort(lambda[lambda > 0], decreasing = T)
        lambdamax = max(lambda)
    }
    if (lambdamax > 1e6) { 
        stop('lambdamax = ', lambdamax, ' is too large.\n')
    }
    if (lambdamax < 1e-3) { 
        warning('lambdamax = ', lambdamax, 
                ' may be too small due to overfitting. Try with smaller random effects.\n')
    }
    nlambda = length(lambda)
    
    # compute the solution path
    if (algo == 'IRLS') {
        res = GLMM.Reg.IRLS(ZX, y, offset, extend.X=extend.X, pZ=pZ, Omega=Omega, family=family, 
                            PQL=PQL, lambda=lambda, include=include, b0=beta.ini, u0=u.ini,
                            maxiter=maxiter, eps=eps, verbose=verbose, time.limit=time.limit,
                            sparse=sparse)
        timing =timing + res$timing
    } else {
        res = HUG::GLMM_Reg_GD(ZX, y, Omega, family, 'L1', as.double(lambda), include-1, 
                               b_ini=beta.ini, u_ini=u.ini, accurate=accu.grad, 
                               LS=T, FAST=T, RS=T, maxiter=maxiter, eps=eps, 
                               maxiter_U=200, eps_U=1e-4, verbose=verbose)
    }
    Beta    = res$B
    U       = res$U
    nlambda = ncol(Beta)
    
    # Stage 2. Refit
    beta.opt  = numeric(p)
    u.opt     = numeric(n)
    time.old  = Sys.time()
    eval.all  = Eval.GLMM(y, ZX, offset, family, Omega, Beta, U, pZ, gamma.ebic)
    if (method.refit != 'none') {
        ix.sel    = which.min(eval.all[,criterion1])
        fm2 = ReFitGLMM(ZX, y, offset, Omega, family=family, algo=algo, accu.grad=accu.grad,
                        PQL=PQL,  include=include,
                        beta.ini=as.matrix(Beta[,ix.sel]), u.ini=as.matrix(U[,ix.sel]),
                        method.refit=method.refit, nthres=nthres,
                        criterion=criterion2, gamma.ebic=gamma.ebic,
                        maxiter=maxiter, eps=eps, verbose=verbose)
        beta.opt  = fm2$beta / c(sd.Z, sd.X)
        u.opt     = fm2$u
    }
    timing[4] = timing[4] + as.numeric(Sys.time() - time.old, 'secs')
    
    list(beta=beta.opt, u=u.opt, Beta=Beta / c(sd.Z, sd.X), U=U, lambda=lambda, 
         df=Matrix::colSums(Beta[-(1:pZ), ] != 0), eval=eval.all, timing=timing, family=family)
}


# Fit a sequence of GLM mixed models via IRLS (and GMD)
GLMM.Reg.IRLS = function(ZX, y, offset = NULL, extend.X = F, pZ = length(include), Omega, 
                         family = c('gaussian', 'binomial', 'poisson'), 
                         PQL = F, lambda, include = integer(0), b0 = NULL, 
                         u0 = NULL, maxiter = 500, eps = 1e-3, verbose = F, time.limit = Inf, 
                         sparse = F) {
    family = match.arg(family)
    ZX    = as.matrix(ZX)
    nA    = length(y)
    n     = ifelse(is.null(Omega), nA, nrow(Omega))
    p     = ncol(ZX)
    pX    = p - pZ
    jx.new = 1:p
    if (extend.X) {
        pX.0 = pX / 2
        jx.new[pZ + (1:pX.0)*2 - 1] = pZ + (1:pX.0)
        jx.new[pZ + (1:pX.0)*2]     = pZ + (1:pX.0) + pX.0
    }
    if (is.null(offset)) { offset = numeric(nA) }
    nlambda = length(lambda)
    ybar    = mean(y)
    mean.ZX = colMeans(ZX); mean.ZX[1] = 0
    ZX      = ZX - matrix(mean.ZX, nA, p, byrow = T) # destroys the orthogonality
    X.ext = Matrix(ZX[,jx.new], sparse = F)
    y.ext = numeric(nA)
    Beta  = Matrix(0, nrow = p, ncol = nlambda)
    U     = matrix(NA, n, nlambda)
    if (!is.null(Omega)) {
        is.diag.Omega = isDiagonal(Omega)
        sparse  = sparse || is.diag.Omega
        if (sparse) {
            Omega = as(Omega, 'dgCMatrix')
        } else {
            Omega = as.matrix(Omega)
        }
        if (n > nA) {
            Omega.2 = Matrix::solve(Omega[-(1:nA), -(1:nA)], Omega[-(1:nA), 1:nA])
            Omega.1 = Omega[1:nA, 1:nA] - Omega[1:nA, -(1:nA)] %*% Omega.2
            Omega.2 = as.matrix(Omega.2)
        } else {
            Omega.1 = Omega
            Omega.2 = matrix(0, 0, n)
        }
        Sigma.1 = Matrix::solve(Omega.1)
    } else {
        sparse  = T
        Omega.1 = NULL
        Omega.2 = NULL
        Sigma.1  = NULL
        PQL     = FALSE
    }
    if (is.null(b0)) {
        beta.now = numeric(p)
        beta.now[include] = coef(glm(y ~ ZX[,include] - 1, family = family))
    } else {
        stopifnot(length(b0) == p)
        beta.now = b0
    }
    if (is.null(u0)) {
        u.now = numeric(n)    
    } else {
        stopifnot(length(u0) == n)
        u.now = u0
    }
    Xbeta.now = drop(X.ext %*% beta.now)
    beta.old  = beta.now
    u.old     = u.now
    derv      = rep(NA, p)
    timing    = numeric(6)
    names(timing) = c('GMD', 'Newton', 'WholePath', 'Refit', 'Prep1', 'Prep2')

    # -------------------------- GMD Parameters ----------------------------
    if (extend.X) {
        group  = as.integer(c(1:pZ, rep((1:pX.0)+pZ, each = 2)))
        bs     = as.integer(c(rep(1, pZ), rep(2, pX.0)))
        ix     = as.integer(c(1:pZ, (1:pX.0)*2 + pZ - 1))
        pf     = c(numeric(pZ), rep(1, pX.0))
    } else {
        group  = as.integer(1:p)
        bs     = as.integer(rep(1, p))
        ix     = as.integer(1:p)
        pf     = c(numeric(pZ), rep(1, pX))
    }
    iy     = as.integer(ix + bs - 1)
    bn     = as.integer(max(group))
    pf[intersect(include, 1:bn)] = 0
    pf     = as.double(pf)
    nobs   = as.integer(nA)
    nvars  = as.integer(p)
    vnames = colnames(X.ext)
    gam    = rep(NA, bn)
    maxit  = as.integer(1e5)
    dfmax  = as.integer(min(nA, nvars))
    pmax   = as.integer(dfmax)
    eps.gg = as.double(1e-6)
    flmin  = as.double(1)
    nlam   = as.integer(1)
    intr   = as.integer(0) # X.ext already contains separate constant columns
    
    # -------------------------- Solution Path ----------------------------
    time.oold = Sys.time()
    # Stage 1. Regularized estimation
    for (i.lam in 1:nlambda) {
        if (verbose) {
            cat('lambda =', lambda[i.lam], ',', i.lam, '/', nlambda, '\n')
        }
        diffs.old  = rep(Inf, 2)
        incr.diffs = 0
        
        for (iter1 in 1:maxiter) {
            time.old   = Sys.time()
            converge1  = F
            beta.old   = beta.now
            u.old      = drop(u.now)
            Xbeta.old  = Xbeta.now
            mu.old     = exp(Xbeta.old + offset + u.old[1:nA])
            mean.old   = switch(family, 
                                binomial = mu.old / (1 + mu.old),
                                poisson  = mu.old, 
                                gaussian = Xbeta.old + offset + u.old[1:nA])
            dy.old     = switch(family, 
                                binomial = mean.old * (1 - mean.old),
                                poisson  = mu.old,
                                gaussian = rep(1 / mean((y - Xbeta.old - offset - u.old)^2), nA))
            # Logistic regression may not have a solution
            if (anyNA(dy.old)) break
            if (family == 'gaussian') {
                y.ext.old = y - u.old
            } else {
                if (PQL) {
                    y.ext.old = Xbeta.old + u.old[1:nA] + (y - mean.old) / dy.old
                } else {
                    y.ext.old = Xbeta.old + (y - mean.old) / dy.old
                }
            }
            if (PQL & !is.diag.Omega) {
                W.old     = Matrix(t(chol(Sigma.1 + diag(1 / dy.old))), sparse = sparse)
                sqrtW.old = as.matrix(Matrix::solve(W.old) * sqrt(2))
                X.ext.old = as.matrix(sqrtW.old %*% X.ext)
                y.ext.old = as.vector(sqrtW.old %*% y.ext.old)
            } else {
                if (PQL & is.diag.Omega) {
                    w.old = 2 / (diag(Sigma.1) + 1 / dy.old)
                }
                if (!PQL) {
                    w.old = ifelse(family == 'gaussian', 1, 2) * dy.old # adjust for gglasso
                }
                sqrtw.old = sqrt(w.old)
                X.ext.old = as.matrix(X.ext * sqrtw.old)    # obs i inactive if W[i] is 0
                y.ext.old = y.ext.old * sqrtw.old
            }
            timing[5] = timing[5] + as.numeric(Sys.time() - time.old, 'secs')
            
            # gam[j] is the max eigenvalue of the j-th 2x2 or 1x1 crossprod matrix
            time.old   = Sys.time()
            if (i.lam == 1) {
                if (PQL & !is.diag.Omega & extend.X) {
                    gam = foreach(j = 1:bn, .combine = c) %do% {
                        max(svd(X.ext.old[,ix[j]:iy[j]])$d^2)
                    }
                } else {
                    # Each group of columns in X.ext are orthogonal
                    # ONLY WHEN X IS NOT CENTERED!!!
                    gams = colSums(X.ext.old^2)
                    gam = foreach(j = 1:bn, .combine = c) %do% {
                        max(gams[ix[j]:iy[j]])
                    }
                }
                gam = as.double(gam / nobs)
            }
            timing[6]  = timing[6] + as.numeric(Sys.time() - time.old, 'secs')
            
            # compute beta by GMD/LS
            time.old  = Sys.time()
            if (lambda[i.lam] == 0) {
                output.gg = try(list(beta = matrix(coef(lm(y.ext.old ~ X.ext.old - 1)))))
            } else {
                theta.ini = as.double(c(0, beta.old))
                output.gg = try(gglasso.mls(bn, bs, ix, iy, gam, nobs, nvars, X.ext.old, y.ext.old, 
                                            pf, dfmax, pmax, as.integer(1), flmin, lambda[i.lam], 
                                            theta.ini, eps.gg, maxit, vnames, group, intr, 
                                            as.double(min(time.limit / 2, 300))))
            }
            if ((class(output.gg) == 'try-error') || (ncol(output.gg$beta) < 1)) {
                beta.now = beta.old
                break
            }
            beta.now   = output.gg$beta
            # cat('beta.now = ', paste0(round(beta.now, 4), collapse = ', '), '\n')
            Xbeta.now  = drop(X.ext %*% beta.now)
            timing[1]  = timing[1] + as.numeric(Sys.time() - time.old, 'secs')
            
            # compute u.tilde
            time.old  = Sys.time()
            if (!is.null(Omega)) {
                # prevent too large u (which results in divergence or local optim)
                u.bound = rep(Inf, nA)
                if (family == 'poisson') {
                    u.bound = log(y) - Xbeta.now
                    u.bound[is.na(u.bound) | (u.bound < 0)] = 0
                }
                if (PQL) {
                    if (is.diag.Omega) {
                        u.now = 0.5 * diag(Sigma.1) * sqrtw.old * 
                            (y.ext.old - sqrtw.old * Xbeta.now)
                    } else {
                        u.now = as.vector(0.5 * Sigma.1 %*% 
                                              (t(sqrtW.old) %*% (y.ext.old-sqrtW.old%*%Xbeta.now)))
                    }
                    u.now = ifelse(u.now > u.bound, u.bound, u.now)
                    if (n > nA) { u.now = c(u.now, -as.vector(Omega.2 %*% u.now)) }
                } else {
                    u.now = try(FindU.GLMM(Omega.1, Omega.2, family, y, Xbeta.now,
                                           u.old[1:nA], 200, 1e-4))
                }
            } else {
                u.now = numeric(n)
            }
            # cat('u.now = ', paste0(round(u.now[seq(1, 41, by = 4)], 4), collapse = ', '), '\n')
            if (class(u.now) == 'try-error') {
                u.now = u.old
                break
            }
            timing[2] = timing[2] + as.numeric(Sys.time() - time.old, 'secs')
            
            diffs     = c(sum(abs(beta.now - beta.old)), sum(abs(u.now - u.old))) / 
                (1 + c(sum(abs(beta.old)), sum(abs(u.old))))
            # cat('Iter ', iter1, ', diffs = ', diffs, '\n')
            
            if (all(diffs < eps)) {
                converge1 = T
                break
            }
            incr.diffs = ifelse(all(diffs > diffs.old, na.rm = T), incr.diffs + 1, 0)
            if (incr.diffs >= 5) {
                cat('IRLS diverges after ', iter1, ' iterations.\n')
                break
            }
            diffs.old = diffs
            if (as.numeric(Sys.time() - time.oold, 'secs') > time.limit) { break }
        }
        
        Beta[jx.new, i.lam] = beta.now
        U[,i.lam]    = u.now
        # to determine the max lambda
        if (lambda[i.lam] > 0) {
            derv = t(X.ext.old) %*% (drop(X.ext.old %*% beta.now) - y.ext.old) / nA
        }
        if ((!converge1) & (iter1 == maxiter)) {
            cat('IRLS not converging within ', maxiter, ' iterations.\n')
        } else if (converge1 & verbose) {
            cat('IRLS converges at iteration ', iter1, '.\n')
        }
        if (as.numeric(Sys.time() - time.oold, 'secs') > time.limit) {
            cat('Algorithm stopped after exceeding ', time.limit, ' seconds. ')
            if (nlambda > i.lam) { cat('Neglect lambda < ', lambda[i.lam]) }
            cat('\n')
            Beta  = Matrix(Beta[,1:i.lam], sparse = T)
            U     = as.matrix(U[,1:i.lam])
            nlambda = i.lam
            break
        }
    }
    timing[3] = timing[3] + as.numeric(Sys.time() - time.oold, 'secs')
    
    Beta[1, ]  = Beta[1, ]  - t(Beta[-1, ]) %*% mean.ZX[-1]
    return(list(B = Matrix(Beta, sparse = T), U = U, nlambda = nlambda, 
                timing = timing, derv = derv))
}

# Refit a GLM (mixed) model
ReFitGLMM = function(X, y, offset = NULL, Omega = NULL, 
                     family = c('binomial', 'poisson', 'gaussian'), 
                     algo = c('FISTA', 'IRLS'), accu.grad = T, 
                     PQL = F, include = NULL, beta.ini, u.ini = NULL, 
                     method.refit = c('thres', 'mle'), nthres = 20, 
                     criterion = c('loglik', 'aic', 'bic', 'ebic'), gamma.ebic = 0.5,
                     maxiter = ifelse(algo == 'FISTA', 1e3, 5e2), eps = 1e-3, 
                     verbose = F) {
    family = match.arg(family)
    p      = length(beta.ini)
    nA     = nrow(X)
    if (length(Omega) == 0) { Omega = NULL }
    n      = ifelse(is.null(Omega), nA, nrow(Omega))
    if (!is.null(u.ini)) {
        u.ini = as.matrix(u.ini)
    } else {
        u.ini = matrix(0, n, 1)
    }
    beta.opt = numeric(p)
    eval.opt = rep(Inf, 4)
    names(eval.opt) = c('loglik', 'aic', 'bic', 'ebic')
    
    if (method.refit == "mle") {
        ix.nz = union(include, which(abs(beta.ini) > 1.0e-8))
        fm    = FitGLMM(as.matrix(X[,ix.nz]), y, offset, Omega, family, algo = algo, 
                        accu.grad = accu.grad, PQL = PQL, beta.ini = NULL, 
                        u.ini = NULL, maxiter = maxiter, eps = eps, verbose = verbose)
        beta.opt[ix.nz] = fm$beta
        u.opt = fm$u
        eval.opt = Eval.GLMM(y, X, offset, family, Omega, beta.opt, u.opt, 1, gamma.ebic)
    } else if (method.refit == "thres") {
        beta.nz = abs(beta.ini[abs(beta.ini) > 1.0e-8])
        if (length(beta.nz) == 0) {
            beta.nz = 0
        }
        if (length(beta.nz) > nthres) {
            thres.beta = seq(min(beta.nz), max(beta.nz), length.out = nthres)
        } else {
            thres.beta = sort(beta.nz)
        }
        for (ithres in 1:length(thres.beta)) {
            cat('ithres / nthres: ', ithres, '/', length(thres.beta), ' \n')
            beta.thres = numeric(p)
            ix.ithres  = which(abs(beta.ini) >= thres.beta[ithres])
            ix.ithres  = union(ix.ithres, include)
            beta.now   = numeric(p)
            fm    = FitGLMM(as.matrix(X[,ix.ithres]), y, offset, Omega, family, algo = algo, 
                            accu.grad = accu.grad, PQL = PQL, beta.ini = NULL, 
                            u.ini = NULL, maxiter = maxiter, eps = eps, verbose = verbose)
            beta.now[ix.ithres] = fm$beta
            u.now = fm$u
            eval.now = Eval.GLMM(y, X, offset, family, Omega, beta.now, u.now, 1, gamma.ebic)
            eval.now[is.na(eval.now)] = Inf
            if (eval.now[1, criterion] < eval.opt[criterion]) {
                beta.opt = beta.now
                u.opt    = u.now
                eval.opt = eval.now[1, ]
            }
        }
    }
    return(list(beta = beta.opt, u = u.opt, eval = eval.opt))
}


# Fit a GLM mixed model without penalty
FitGLMM = function(X, y, offset = NULL, Omega, family = c('binomial', 'poisson', 'gaussian'), 
                   algo = c('FISTA', 'IRLS'), accu.grad = T, PQL = F, 
                   beta.ini = NULL, u.ini = NULL, 
                   maxiter = ifelse(algo == 'FISTA', 1e3, 5e2), eps = 1e-3, verbose = F) {
    family = match.arg(family)
    X      = as.matrix(X)
    y      = as.vector(y)
    # Re-scaling
    sd.X   = apply(X, 2, sd)
    if (any(sd.X == 0)) { sd.X[which(sd.X == 0)[1]] = 1 }
    X      = X %*% diag(1 / sd.X)
    X[, sd.X == 0] = 0
    
    if (length(Omega) == 0) {
        res = try(glm(y ~ X - 1, family = family))
        if (class(res)[1] == 'try-error') {
            warning('Cannot fit an unpenalized GLM model. Try with a small penalty.')
            res = glmnet::glmnet(X, y, family = family, intercept = F, offset = offset)
            res$beta = as.vector(res$beta[,ncol(res$beta)])
        } else {
            res$beta = res$coefficients
        }
        res$u = numeric(length(y))
    } else {
        res = try(switch(algo,
                         IRLS = GLMM.Reg.IRLS(X, y, offset, F, ncol(X), Omega, family, 
                                              PQL=PQL, lambda=0, include=1:ncol(X), b0=beta.ini, 
                                              u0=u.ini, maxiter=maxiter, eps=eps, verbose=verbose),
                         FISTA = HUG::GLMM_Reg_GD(X, y, Omega, family, 'none', 0, 
                                                  include=1:ncol(X)-1, 
                                                  b_ini=beta.ini, u_ini=u.ini, 
                                                  accurate=accu.grad, LS=T, FAST=T, RS=T, 
                                                  maxiter=maxiter, eps=eps, maxiter_U=200, 
                                                  eps_U=1e-4, verbose=verbose)))
        if ((class(res) == 'try-error') || anyNA(res$B)) {
            warning('Cannot fit an unpenalized GLMM model. Try with a small penalty.')
            res = switch(algo,
                         IRLS = GLMM.Reg.IRLS(X, y, offset, F, 1, Omega, family, PQL=PQL, 
                                              lambda=1e-10, include=1, b0=beta.ini, u0=u.ini, 
                                              maxiter=maxiter, eps=eps, verbose=verbose),
                         FISTA = HUG::GLMM_Reg_GD(X, y, Omega, family, 'none', 1e-10, include=0, 
                                                  b_ini=beta.ini, u_ini=u.ini, 
                                                  accurate=accu.grad, LS=T, FAST=T, RS=T, 
                                                  maxiter=maxiter, eps=eps, maxiter_U=200, 
                                                  eps_U=1e-4, verbose=verbose))
        }
        res$beta = as.vector(res$B) / sd.X
        res$u    = as.vector(res$U)
    }
    return(res)
}

# Evaluate a bunch of GLM (mixed) models
# X contains an constant column
# The first row of B are the intercepts
Eval.GLMM = function(y, X, offset = NULL, family = c('gaussian', 'binomial', 'poisson'),
                     Omega = NULL, B, U = NULL, p2 = 1, gamma.ebic = 0.5) {
    family = match.arg(family)
    if (length(Omega) == 0) { Omega = NULL }
    nA    = nrow(X)   # possibly < n
    n     = ifelse(is.null(Omega), nA, nrow(Omega))
    p     = ncol(X) - p2
    B     = as.matrix(B)
    K     = ncol(B)
    stopifnot(nrow(B) == ncol(X))
    if (gamma.ebic == 'auto') {
        gamma.ebic = 1 - log(p) / log(nA) / 2
    }
    if (!is.null(Omega)) {
        stopifnot(!is.null(U), ncol(U) == K, nrow(U) == n)
        U = as.matrix(U)
        ldet.Omega = Matrix::determinant(Omega)[[1]]
    }
    if (is.null(offset)) { offset = numeric(nA) }
    
    p.nz  = colSums(as.matrix(B[-(1:p2),]) != 0)
    XB    = as.matrix(X %*% B) + matrix(offset, nA, K, byrow = F)
    if (!is.null(Omega)) {
        logmu = XB + U[1:nA, ]
    } else {
        logmu = XB
    }
    mu    = exp(logmu)
    if (family == 'binomial') {
        mean.all = mu / (1 + mu)
        var.all  = mean.all * (1 - mean.all)
        loglik   = as.vector(t(y) %*% logmu) - colSums(log(1 + mu))
    } else if (family == 'gaussian') {
        mean.all = logmu
        var.all  = colMeans((matrix(y, nA, K, byrow = F) - logmu)^2)
        loglik   = -nA/2*log(var.all)
        var.all  = matrix(1 / var.all, nA, K, byrow = T)
    } else if (family == 'poisson') {
        mean.all = mu
        var.all  = mu
        loglik   = as.vector(t(y) %*% logmu) - colSums(mu)
    } else {
        stop(family, ' family not supported yet.')
    }
    if (!is.null(Omega)) {
        Omega.tmp = Omega
        for (i in 1:K) {
            diag(Omega.tmp)[1:nA] = diag(Omega)[1:nA] + var.all[,i]
            loglik[i] = loglik[i] - as.numeric(t(U[,i]) %*% Omega %*% U[,i]) / 2 + 
                ldet.Omega / 2 - Matrix::determinant(Omega.tmp)[[1]] / 2 
        }
    }
    
    eval.all = matrix(NA, K, 4)
    eval.all[,1] = -loglik
    eval.all[,2] = 2*eval.all[,1] + p.nz*2
    eval.all[,3] = 2*eval.all[,1] + p.nz*log(nA)
    eval.all[,4] = eval.all[,3] + 2*gamma.ebic * 
        (lgamma(p+1) - lgamma(p.nz+1) - lgamma(p-p.nz+1))
    colnames(eval.all) = c('loglik', 'aic', 'bic', 'ebic')
    rownames(eval.all) = colnames(B)
    
    return(data.frame(eval.all))
}

EvalPen = function(G, B, lambda) {
    sqrt(colSums(G^2) + colSums(B^2)) * lambda
}

# For the usage of C++ functions
as.dgC = function(x) as(x, 'dgCMatrix')

FindU.GLMM = function(Omega1, Omega2, family, y, Xbeta, offset = NULL, uA_ini, maxiter, eps) {
    nA  = nrow(Omega1)
    nAc = nrow(Omega2)
    n   = nA + nAc
    if (is.null(offset)) { offset = numeric(nA) }
    stopifnot(length(uA_ini) == nA, length(Xbeta) == nA, length(offset) == nA, 
              length(y) == nA, ncol(Omega1) == nA, ncol(Omega2) == nA)
    
    u       = numeric(n)
    uA      = uA_ini
    S_derv2 = Omega1
    is_conv = F
    # prevent too large u (which results in divergence or local optim)
    u.bound = rep(Inf, nA)
    if (family == 'poisson') {
        u.bound = log(y) - Xbeta
        u.bound[is.na(u.bound) | (u.bound < 0)] = 0
    }
    
    for (iter in 1:maxiter) {
        # Analytical form exists for Gaussian case
        if (family == "gaussian") {
            mu     = Xbeta + uA + offset
            varY   = mean((y - mu)^2)
            diag(S_derv2) = diag(Omega1) + 1 / varY
            uA     = as.vector(Matrix::solve(S_derv2, y - Xbeta)) / varY
            is_conv = T
            break
        }
        if (family == "binomial") {
            mu     = 1 / (1 + exp(-(Xbeta + uA + offset)))
            varY   = mu * (1 - mu)
            S_derv = as.vector(Omega1 %*% uA) + mu - y
            diag(S_derv2) = diag(Omega1) + varY
        } else if (family == "poisson") {
            mu     = exp(Xbeta + uA + offset)
            varY   = mu
            S_derv = as.vector(Omega1 %*% uA) + mu - y
            diag(S_derv2) = diag(Omega1) + varY
        }
        uA_diff = as.vector(Matrix::solve(S_derv2, S_derv))
        uA      = uA - uA_diff
        ix.OOR  = (uA > u.bound)
        if (any(ix.OOR)) {
            warning('Possible outlier: obs ', which.max(abs(S_derv)),
                    ' (extremely influential)')
            uA[ix.OOR] = u.bound[ix.OOR]
        }
        if (sum(abs(uA_diff)) <= eps * (1.0 + sum(abs(uA)))) {
            is_conv = T
            break
        }
    }
    if  (!is_conv) {
        cat("Newton-Raphson for u_0 updating may not converge within ", maxiter, " iterations.\n")
    }
    
    if (nAc == 0) {
        u = uA
    } else {
        u[1:nA] = uA
        u[-(1:nA)] = -as.vector(Omega2 %*% uA)
    }
    return(u)
}

ConditionalCenter = function(X) {
    stopifnot(!anyNA(X))
    X[X == 0] = NA
    mean.X = colMeans(X, na.rm = T)
    mean.X[is.na(mean.X)] = 0 # in case of all-NA X[,j]
    X = X - matrix(mean.X, nrow(X), ncol(X), byrow = T)
    X[is.na(X)] = 0
    return(list(X = X, mean = mean.X))
}

