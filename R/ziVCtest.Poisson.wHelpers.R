ziVCtest.Poisson.Linear =
  function (res, Z, X1, kernel, weights = NULL, method, res.out,
            n.Resampling, r.corr, IsMeta = FALSE){
    if (length(r.corr) > 1 ) {
      stop("SKATO not implemented for ZIP!")
    }
    if (IsMeta) {
      stop("Meta (IsMeta = TRUE) not implemented for ZIP!")

    }
    else if (length(r.corr) == 1) {
      re = ziKMTest.Poisson.Linear(res = res, Z = Z, X1 = X1,
                                   kernel = kernel, weights = weights,
                                   method = method, res.out = res.out,
                                   n.Resampling = n.Resampling, r.corr = r.corr)
    } else {
      stop("method \"optimal\" not implemented for ZIP!")
      # re = SKAT_Optimal_Logistic(res, Z, X1, kernel, weights,
      # pi_1, method, res.out, n.Resampling, r.corr)
    }
    return(re)
  }

ziKMTest.Poisson.Linear = function(res, Z, X1, kernel, weights = NULL,
                                   method = "davies.default", res.out,
                                   n.Resampling, r.corr){
  # for testing:
  # null = ziVCtest_null
  # res = null$res          # (residuals)
  # res.out = null$res.out  # (resampled residuals)
  # res.out = ziVCtest_null$res.out

  # Z = SKAT.example.ZIP$Z
  # X1 = NULL
  # kernel = "linear"
  # weights = NULL
  # method = "davies.default"
  # n.Resampling = 100; r.corr = 0

  # Z = scale(Z,scale = F)

  start = Sys.time()

  if (kernel == "linear.weighted") {
    Z = t(t(Z) * (weights))
  }
  if (r.corr == 1) {
    Z <- cbind(rowSums(Z))
  } else if (r.corr > 0) {
    stop("SKATO not implemented")
    # p.m <- dim(Z)[2]
    # R.M <- diag(rep(1 - r.corr, p.m)) + matrix(rep(r.corr,
    #                                                p.m * p.m), ncol = p.m)
    # L <- chol(R.M, pivot = TRUE)
    # Z <- Z %*% t(L)
  }

  # ziVCtest_null$res[,c("zero","count")] / as.matrix(as.data.frame(resids))

  # Q.Temp = t(res) %*% Z
  # Q = Q.Temp %*% t(Q.Temp)/2
  # Q.res = NULL
  # if (n.Resampling > 0) {
  #   Q.Temp.res = t(res.out) %*% Z
  #   Q.res = rowSums(rbind(Q.Temp.res^2))/2
  # }
  # W.1 = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*% X1) %*% solve(t(X1) %*%
  #                                                              (X1 * pi_1)) %*% (t(X1) %*% (Z * pi_1))

  n = nrow(res)
  q = ncol(Z)

  #SCORE HATS
  res.pi    = res[ , "zero", drop = F  ]
  res.ld    = res[ , "count",  drop = F   ]
  S.pi = c(t(res.pi) %*% Z)
  S.ld = c(t(res.ld) %*% Z)
  # Q.pi = t(S.pi) %*% S.pi
  # Q.ld = t(S.ld) %*% S.ld


  if(length(S.pi) != length(S.ld) || length(res.pi) != length(res.pi) ){
    stop ("Score components do not match dimensionally.")
  }

  # SCORE STARS:
  resids.pi = res.out[ , "zero" , ]
  resids.ld = res.out[ , "count", ]
  S.rsm.pi = t(t(resids.pi) %*% Z)
  S.rsm.ld = t(t(resids.ld) %*% Z)
  # Q.str.subtr.hat.pi = apply(S.rsm.pi-S.pi, 2, function(x) sum(x^2))
  # Q.str.subtr.hat.ld = apply(S.rsm.ld-S.ld, 2, function(x) sum(x^2))


  #variance-standardized combined score:
  scale.pi = sqrt(sum(apply(S.rsm.pi, 1, var)))
  scale.ld = sqrt(sum(apply(S.rsm.ld, 1, var)))
  S.cmb = c((1/scale.pi)*S.pi, (1/scale.ld)*S.ld)
  S.rsm.cmb = rbind((1/scale.pi)*S.rsm.pi,(1/scale.ld)*S.rsm.ld)
  # Q.cmb = sum(S.cmb^2)
  # Q.str.subtr.hat.cmb = apply(S.rsm.cmb-S.cmb, 2, function(x) sum(x^2))

  # results elements
  out = list()

  pval.res.pi = Get_CompQuad_Pval(S.pi, S.rsm.pi, method = method)
  p.val.pi = out$p.value.pi = pval.res.pi$p.val
  out$p.value.pi.rsmp = pval.res.pi$p.val.resamp

  pval.res.ld = Get_CompQuad_Pval(S.ld, S.rsm.ld, method = method)
  p.val.ld = out$p.value.ld = pval.res.ld$p.val
  out$p.value.ld.rsmp = pval.res.ld$p.val.resamp

  pval.res.cmb = Get_CompQuad_Pval(S.cmb, S.rsm.cmb, method = method)
  out$p.value.cmb = pval.res.cmb$p.val
  out$p.value.cmb.rsmp = pval.res.cmb$p.val.resamp


  # Combined p-values
  # Empirical Brown's method
  library(EmpiricalBrownsMethod)
  # (this method uses correlation in these scores to adjust the scaled chisq constant)
  resid.pi.ld = cbind(res.pi,res.ld)
  Browns.result = empiricalBrownsMethod(data_matrix = t(resid.pi.ld),
                                        p_values = c(out$p.value.pi,out$p.value.ld), extra_info = T)
  p.val.B = Browns.result$P_test
  # Note that $P_test often equal to $P_Fisher


  # Min-p: p.val.distribution for p.val.m's
  p.values.m.a = apply(rbind(out$p.value.pi.rsmp, out$p.value.ld.rsmp), 2, min)
  p.values.F.a = apply(rbind(out$p.value.pi.rsmp, out$p.value.ld.rsmp), 2, function(x) log(x[1]) + log(x[2]))

  #as or more extreme in the direction away from the null:
  # the following two lines is replaced by below because of NA's
  # p.val.m = sum(p.values.m.a <= as.vector(min(p.val.1,p.val.2)))/B
  # p.val.F = sum(p.values.F.a <= as.vector(log(p.val.1) + log(p.val.2)))/B
  p.val.m = mean(p.values.m.a <= as.vector(min(p.val.pi,p.val.ld)), na.rm = TRUE)
  p.val.F = mean(p.values.F.a <= as.vector(log(p.val.pi) + log(p.val.ld)), na.rm = TRUE)

  verbose = FALSE
  if(verbose){
    print('p-val time interval')
    print(Sys.time()-start)
    Sys.time()->start;
  }


  # final results list
  ###########################
  re <- list(p.val.pi =  out$p.value.pi,
             p.val.ld =  out$p.value.ld,
             p.val.cmb = out$p.value.cmb,

             p.val.pi.resamp  = out$p.value.pi.rsmp,
             p.val.ld.resamp  = out$p.value.ld.rsmp,
             p.val.cmb.resamp = out$p.value.cmb.rsmp,

             p.val.Fsr = p.val.F,
             p.val.min = p.val.m,

             Test.Type = method,
             Q.pi = pval.res.pi$Q, Q.ld = pval.res.ld$Q,
             Q.resampling.pi = pval.res.pi$Q.resamp,
             Q.resampling.ld = pval.res.ld$Q.resamp,
             param = out$param)
  return(re)
}

Get_CompQuad_Pval = function(S, S.res, method = c('davies.default', 'davies.precise', 'farebrother', 'imhof', 'liu'), ... ){

  method = match.arg(method)

  # p.vals using compquadform:
  S.star.subtr.hat = S.res - S
  S.var = var(t(S.star.subtr.hat))

  Q = t(S) %*% S
  Q.str.subtr.hat = apply( S.star.subtr.hat, 2, function(x) sum(x^2))

  eig.S.var = eigen(S.var, symmetric = TRUE)
  Lds.S.var = eig.S.var$values
  scl.Lda = Lds.S.var[1]

  if(method == 'davies.default'){
    p.rslt =
      CompQuadForm::davies(q=Q, lambda = Lds.S.var)
    p.vals = sapply(Q.str.subtr.hat, function(x)
      CompQuadForm::davies(q=x, lambda = Lds.S.var)$Qq )
    p.val = p.rslt$Qq


  } else if (method == 'davies.precise'){

    p.rslt =
      CompQuadForm::davies(q=Q/scl.Lda, lambda = Lds.S.var/scl.Lda,
                           lim = 5e7, acc = 1e-12)
    p.vals = sapply(Q.str.subtr.hat, function(x)
      CompQuadForm::davies(q=x/scl.Lda, lambda = Lds.S.var/scl.Lda,
                           lim = 5e7, acc = 1e-12)$Qq )
    p.val = p.rslt$Qq


  } else if (method == 'farebrother'){
    p.rslt =
      CompQuadForm::farebrother(q=Q/scl.Lda, lambda = Lds.S.var/scl.Lda,
                                maxit = 100, eps = 1e-12 )

    p.vals = sapply(Q.str.subtr.hat, function(x)
      CompQuadForm::farebrother(q=x/scl.Lda, lambda = Lds.S.var/scl.Lda,
                                maxit = 100, eps = 1e-12 )$Qq )
    p.val = p.rslt$Qq


  } else if (method == 'imhof'){
    p.rslt =
      CompQuadForm::imhof(q=Q/scl.Lda, lambda = Lds.S.var/scl.Lda,
                          epsabs = 1e-4, epsrel = 1e-5, limit = 1e4)
    p.vals = sapply(Q.str.subtr.hat, function(x)
      CompQuadForm::imhof(q=x/scl.Lda, lambda = Lds.S.var/scl.Lda,
                          epsabs = 1e-4, epsrel = 1e-5, limit = 1e4)$Qq   )
      p.val = p.rslt$Qq


  } else if (method == 'liu'){
    p.val = liu(q=Q/scl.Lda, lambda = Lds.S.var/scl.Lda)
    p.vals = sapply(Q.str.subtr.hat, function(x)
      CompQuadForm::liu(q=x/scl.Lda, lambda = Lds.S.var/scl.Lda)$Qq )
    p.val = p.rslt$Qq


  }
 return(list(p.val=p.val, p.val.resamp = p.vals, Q = Q, Q.resamp = Q.str.subtr.hat) )
}

get.pvals.from.scores = function(score.pi,score.ld, scores.pi.resamp, scores.ld.resamp, method = 'davies', verbose){
  res.pi = ziVCtest_null$res[,"count"]
  res.ld = ziVCtest_null$res[,"zero"]
  resids.pi = ziVCtest_null$res.out[,"zero",]
  resids.ld = ziVCtest_null$res.out[,"count",]

  n = nrow(cbind(res.pi,res.ld))
  q = ncol(Z)

  #SCORE HATS:
  S.pi = c(t(res.pi) %*% Z)
  S.ld = c(t(res.ld) %*% Z)
  Q.pi = sum(S.pi^2)
  Q.ld = sum(S.ld^2)


  # SCORE STARS:
  S.res.pi = t(t(res.out[,'zero',]) %*% Z)
  S.res.ld = t(t(res.out[,'count',]) %*% Z)

  Q.str.subtr.hat.pi = apply(S.res.pi-S.pi, 2, function(x) sum(x^2))
  Q.str.subtr.hat.ld = apply(S.res.ld-S.ld, 2, function(x) sum(x^2))

  #standardized combined score:
  scale.pi = sqrt(sum(apply(S.res.pi, 1, var)))
  scale.ld = sqrt(sum(apply(S.res.ld, 1, var)))
  S.cmb = c((1/scale.pi)*S.pi, (1/scale.ld)*S.ld)
  Q.cmb = sum(S.cmb^2)
  S.res.cmb = rbind((1/scale.pi)*S.res.pi,(1/scale.ld)*S.res.ld)
  Q.str.subtr.hat.cmb = apply(S.res.cmb-S.cmb, 2, function(x) sum(x^2))


  if(verbose){
    print('resample scores time interval')
    print(Sys.time()-start)
    Sys.time()->start;
  }

  # p.vals: direct:
  # p.val.1 = sum(Q.str.subtr.hat.pi >= as.vector(Q.pi) )/B
  # p.val.2 = sum(Q.str.subtr.hat.ld >= as.vector(Q.ld) )/B
  # p.val.c.a = sum(Q.str.subtr.hat.cmb >= as.vector(Q.cmb) )/B

  # p.vals: using compquadform:
  star.subtr.hat.1 = S.res.pi - S.pi
  star.subtr.hat.2 = S.res.ld - S.ld
  star.subtr.hat.c = S.res.cmb - S.cmb
  var.S.1 = var(t(star.subtr.hat.1))
  var.S.2 = var(t(star.subtr.hat.2))
  var.S.c = var(t(star.subtr.hat.c))
  E.var.SS.1 = eigen(var.S.1, symmetric = TRUE)
  E.var.SS.2 = eigen(var.S.2, symmetric = TRUE)
  E.var.SS.c = eigen(var.S.c, symmetric = TRUE)
  Lds.var.S.1 = E.var.SS.1$values
  Lds.var.S.2 = E.var.SS.2$values
  Lds.var.S.c = E.var.SS.c$values

  scl.1 = Lds.var.S.1[1]
  scl.2 = Lds.var.S.2[1]
  scl.c = Lds.var.S.c[1]

  if(method == 'davies'){
    p.val.1 = davies(q=Q.pi, lambda = Lds.var.S.1)$Qq
    p.val.2 = davies(q=Q.ld, lambda = Lds.var.S.2)$Qq
    p.val.c.a = davies(q=Q.cmb, lambda = Lds.var.S.c)$Qq

  } else if (method == 'davies.precise'){
    p.val.1 = davies(q=Q.pi/scl.1, lambda = Lds.var.S.1/scl.1, lim = 5e7, acc = 1e-12)$Qq
    p.val.2 = davies(q=Q.ld/scl.2, lambda = Lds.var.S.2/scl.2, lim = 5e7, acc = 1e-12)$Qq
    p.val.c.a = davies(q=Q.cmb/scl.c, lambda = Lds.var.S.c/scl.c, lim = 5e7, acc = 1e-15)$Qq

  } else if (method == 'farebrother'){

    p.val.1 = farebrother(q=Q.pi/scl.1, lambda = Lds.var.S.1/scl.1,
                            maxit = 100, eps = 1e-12 )$Qq
    p.val.2 = farebrother(q=Q.ld/scl.2, lambda = Lds.var.S.2/scl.2,
                            maxit = 100, eps = 1e-12 )$Qq
    p.val.c.a = farebrother(q=Q.cmb/scl.c, lambda = Lds.var.S.c/scl.c,
                            maxit = 100, eps = 1e-15, mode = 1)$Qq

  } else if (method == 'imhof'){

    p.val.1 = imhof(q=Q.pi/scl.1, lambda = Lds.var.S.1/scl.1,  epsabs = 1e-4, epsrel = 1e-5, limit = 1e4)$Qq; p.val.1
    p.val.2 = imhof(q=Q.ld/scl.2, lambda = Lds.var.S.2/scl.2,  epsabs = 1e-4, epsrel = 1e-5, limit = 1e4)$Qq; p.val.2
    p.val.c.a = imhof(q=Q.cmb/scl.c, lambda = Lds.var.S.c/scl.c,  epsabs = 1e-7, epsrel = 1e-7, limit = 1e4)$Qq; p.val.c.a

  } else if (method == 'liu'){

    p.val.1 = liu(q=Q.pi, lambda = Lds.var.S.1)
    p.val.2 = liu(q=Q.ld, lambda = Lds.var.S.2)
    p.val.c.a = liu(q=Q.cmb, lambda = Lds.var.S.c)

  }
  p.val.1; p.val.2; p.val.c.a

  # p.val.1 = imhof(q=100000*unlist(Q.pi), lambda = 100000*Lds.var.S.1)$Qq
  #Note: imhof method gives p-values of ~0.5 when lambda is set to zero.
  # It also does this when lambda is less than 10^-05. Rescaling fixes this.


  # Combined p-values
  # Empirical Brown's method
  # (this method uses correlation in these scores to adjust the scaled chisq constant)
  scores.1.2.a = as.matrix(as.data.frame(scores.a[c('resid.pi','resid.lda')]))
  # cor(scores.1.2.a); plot(scores.a$resid.lda,scores.a$resid.pi)
  p.val.B = empiricalBrownsMethod(data_matrix = t(scores.1.2.a),
                                    p_values = c(p.val.1,p.val.2), extra_info = T)$P_test


  # Min-p: p.val.distribution for p.val.m's
  # p.values.1.a = sapply(Q.str.subtr.hat.pi, function(x) sum(as.vector( Q.str.subtr.hat.pi ) >= x )/B)
  # p.values.2.a = sapply(Q.str.subtr.hat.ld, function(x) sum(as.vector( Q.str.subtr.hat.ld ) >= x )/B)
  p.values.1.a = sapply(Q.str.subtr.hat.pi, function(x) davies(q=x, lambda = Lds.var.S.1)$Qq )
  p.values.2.a = sapply(Q.str.subtr.hat.ld, function(x) davies(q=x, lambda = Lds.var.S.2)$Qq )
  p.values.m.a = apply(rbind(p.values.1.a,p.values.2.a), 2, min)
  p.values.F.a = apply(rbind(p.values.1.a,p.values.2.a),2, function(x) log(x[1]) + log(x[2]))
  #as or more extreme in the direction away from the null:
  # p.val.m = sum(p.values.m.a <= as.vector(min(p.val.1,p.val.2)))/B
  # p.val.F = sum(p.values.F.a <= as.vector(log(p.val.1) + log(p.val.2)))/B
  p.val.m = mean(p.values.m.a <= as.vector(min(p.val.1,p.val.2)), na.rm = TRUE)
  p.val.F = mean(p.values.F.a <= as.vector(log(p.val.1) + log(p.val.2)), na.rm = TRUE)


  if(verbose){
    print(paste('object.size(score.stars.res.a)',object.size(score.stars.res.a) ))
    print('p-val time interval')
    print(Sys.time()-start)
    Sys.time()->start;
  }

  return(list(
    p.vals = list(
      p.pi.VC = p.val.1,
      p.ld.VC = p.val.2,
      p.cmb.VC = p.val.c.a,
      p.min.VC = p.val.m,
      p.Fsr.VC = p.val.F,
      p.Brn.VC = p.val.B
    ),
    p.stars = list(
      p.stars.pi = p.values.1.a,
      p.stars.ld = p.values.2.a,
      p.stars.min = p.values.m.a,
      p.stars.Fsr = p.values.F.a
    ),
    beta0.hats = list(
      beta0.hat.a = beta0.hat.a,
      beta0.hat.pi.a = beta0.hat.pi.a,
      beta0.hat.lda.a = beta0.hat.lda.a
    )
  ))
}



