ziVCtest.Poisson.Other =
  function (res, Z, X1, kernel, weights = NULL, method, res.out,
            n.Resampling, r.corr, IsMeta = FALSE){
    if (length(r.corr) > 1 ) {
      stop("SKATO not implemented for ZIP!")
    }
    if (IsMeta) {
      stop("Meta (IsMeta = TRUE) not implemented for ZIP!")

    }
    else if (length(r.corr) == 1) {
      re = ziKMTest.Poisson.Other(res = res, Z = Z, X1 = X1,
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


ziKMTest.Poisson.Other = function(res, Z, X1, kernel, weights = NULL,
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

  if (r.corr == 1) {
    stop("SKATO not implemented")
    # Z <- cbind(rowSums(Z))
  } else if (r.corr > 0) {
    stop("SKATO not implemented")
    # p.m <- dim(Z)[2]
    # R.M <- diag(rep(1 - r.corr, p.m)) + matrix(rep(r.corr,
    #                                                p.m * p.m), ncol = p.m)
    # L <- chol(R.M, pivot = TRUE)
    # Z <- Z %*% t(L)
  }


  start = Sys.time()

  n <- nrow(Z)
  q = ncol(Z)


  if (class(kernel) == "matrix") {
    K = kernel
  } else {
    if (kernel == "linear.weighted") {
      Z.w = t(t(Z) * (weights))
      K = Z.w %*% t(Z.w)
    } else {
      K = SKAT:::lskmTest.GetKernel(Z, kernel, weights, n, q)
    }
  }

  #SCORE HATS
  res.pi    = res[ , "zero", drop = F  ]
  res.ld    = res[ , "count",  drop = F   ]
  # S.pi = c(t(res.pi) %*% Z)
  # S.ld = c(t(res.ld) %*% Z)
  # Q.pi = t(S.pi) %*% S.pi
  # Q.ld = t(S.ld) %*% S.ld
  Q.hat.pi = t(res.pi) %*% K %*% res.pi
  Q.hat.pi = t(res.ld) %*% K %*% res.ld


  # SCORE STARS:
  resids.pi = res.out[ , "zero" , ]
  resids.ld = res.out[ , "count", ]
  # S.rsm.pi = t(t(resids.pi) %*% Z)
  # S.rsm.ld = t(t(resids.ld) %*% Z)
  # Q.str.subtr.hat.pi = apply(S.rsm.pi-S.pi, 2, function(x) sum(x^2))
  # Q.str.subtr.hat.ld = apply(S.rsm.ld-S.ld, 2, function(x) sum(x^2))
  Q.star.pi = Q.star.ld = NULL
  if (n.Resampling > 0) {
    Q.star.pi = Q.star.ld = rep(NA, n.Resampling)
    for (i in 1:n.Resampling) {
      Q.star.pi[i] = t(resids.pi[, i]) %*% K %*% resids.pi[, i]
      Q.star.ld[i] = t(resids.ld[, i]) %*% K %*% resids.ld[, i]
    }
  }


  if(length(S.pi) != length(S.ld) || length(res.pi) != length(res.pi) ){
    stop ("Score components do not match dimensionally.")
  }


  #variance-standardized combined score:
  # TODO how to scale when the half matrix is not available due to Kernel
  Z = Z.eign = SKAT:::Get_Matrix_Square.1(K)
  S.pi = c(t(res.pi) %*% Z)
  S.ld = c(t(res.ld) %*% Z)
  S.rsm.pi = t(t(resids.pi) %*% Z)
  S.rsm.ld = t(t(resids.ld) %*% Z)

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

    # the following two lines is replaced by below because of potential NA's
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



SKAT:::SKAT.linear.Linear
#   SKAT:::SKAT.linear.Other =
#   function (res, Z, X1, kernel, weights = NULL, s2, method, res.out,
#             n.Resampling)
#   {
#     n <- nrow(Z)
#     m = ncol(Z)
#     if (class(kernel) == "matrix") {
#       K = kernel
#     }
#     else {
#       K = lskmTest.GetKernel(Z, kernel, weights, n, m)
#     }
#     Q = t(res) %*% K %*% res/(2 * s2)
#     Q.res = NULL
#     if (n.Resampling > 0) {
#       Q.res <- rep(0, n.Resampling)
#       for (i in 1:n.Resampling) {
#         Q.res[i] = t(res.out[, i]) %*% K %*% res.out[, i]/(2 *
#                                                              s2)
#       }
#     }
#     W = K - X1 %*% solve(t(X1) %*% X1) %*% (t(X1) %*% K)
#     if (method == "davies") {
#       W1 = W - (W %*% X1) %*% solve(t(X1) %*% X1) %*% t(X1)
#     }
#     if (method == "liu") {
#       out <- Get_Liu_PVal(Q, W, Q.res)
#       pval.zero.msg = NULL
#     }
#     else if (method == "liu.mod") {
#       out <- Get_Liu_PVal.MOD(Q, W, Q.res)
#       pval.zero.msg = NULL
#     }
#     else if (method == "davies") {
#       out <- Get_Davies_PVal(Q, W1, Q.res)
#       pval.zero.msg = out$pval.zero.msg
#     }
#     else {
#       stop("Invalid Method!")
#     }
#     re <- list(p.value = out$p.value, p.value.resampling = out$p.value.resampling,
#                Test.Type = method, Q = Q, param = out$param, pval.zero.msg = pval.zero.msg)
#     return(re)
#   }
#   <environment: namespace:SKAT>
