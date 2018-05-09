
ziVCtest.Hurdle.Linear =
  function (obj.res, Z,
            kernel, weights = NULL,
            dosages = FALSE,
            method = c("default.hurdle"),
            r.corr, IsMeta = FALSE) {


    if (length(r.corr) > 1 && dim(Z)[2] == 1) {
      r.corr = 0
    }
    if (IsMeta) {
      stop("Meta (IsMeta = TRUE) not implemented!")
      # re = SKAT_RunFrom_MetaSKAT(res = res, Z = Z, X1 = X1,
      #                            kernel = kernel, weights = weights, pi_1 = pi_1,
      #                            out_type = "Hurdle", method = method, res.out = res.out,
      #                            n.Resampling = n.Resampling, r.corr = r.corr)
    }
    else if (length(r.corr) == 1) {
      re = ziKMtest.Hurdle.Linear(obj.res = obj.res, Z = Z,
                                  kernel = kernel, weights = weights,
                                  dosages = dosages,
                                  method = method, r.corr = r.corr)

    }
    else {
      stop("method \"optimal\" not implemented!")
      re = SKAT_Optimal_Logistic(res, Z, X1, kernel, weights,
                                 pi_1, method, res.out, n.Resampling, r.corr)
    }
    return(re)
  }

ziKMtest.Hurdle.Linear =
  function (obj.res, Z,
            kernel, weights,
            dosages,
            method, r.corr,
            p.Brown = TRUE){


    # ziKMTest.Poisson.Linear = function(res, Z, X1, kernel, weights = NULL,
    #                                    method = "davies.default", res.out,
    #                                    n.Resampling, r.corr)
    if(method != "default.hurdle"){
      warning("For Hurdle, only implemented method is \"default.hurdle\":
              method.C = \"davies\", method.bin = \"Hybrid\", method.D = \"SKAT\" ")
    }

    method.C = "davies"
    method.D = "SKAT"
    method.bin = "Hybrid"

    # Using adjusted Binary testing method: for  n<2000 applies less conservative inference by kurtosis adjustment
    S.D.adj = SKATBinary(Z, obj.res$obj.D, kernel = kernel, is_dosage = dosages,
                         method.bin = method.bin, method = method.D)
    print("SKAT binary adjustment method:")
    print(paste(S.D.adj$method.bin))

    p.bn = S.D.adj$p.value
    p.bn.dist = S.D.adj$p.value.resampling
    Q.bn.dist = S.D.adj$Q.resampling


    # Continuous trait testing method
    class(obj.res$obj.C) = "SKAT_NULL_Model_Hurdle"
    # Need to specify class SKAT_NULL_Model_Hurdle to enter vcHurdle_With_NullModel
    # in order to run modified function and return Q.resampling for out.type = "C"
    S.C = VCtest(Z, obj.res$obj.C, kernel = kernel, is_dosage = dosages,
                 method = method.C)
    p.C = S.C$p.value
    p.C.dist = S.C$p.value.resampling
    Q.C.dist = S.C$Q.resampling

    # Combined Continuous and Binary
    p.min.bn = mean( apply(cbind(p.C.dist,p.bn.dist),1, min) <  min(p.C,p.bn) )
    p.Fsr.bn.rsmp = mean( apply(cbind(p.C.dist,p.bn.dist),1, function(x) sum(log(x)))  < sum(log(c(p.C,p.bn))))
    p.Fsr.bn = pchisq(q = -2*sum(log(c(p.C,p.bn))),df=4,lower.tail = FALSE)

    # Try to account for possible correlation in test statistics due to adjusting for covariates
    if(p.Brown){
      corr.data = sqrt(cbind(Q.C.dist,Q.bn.dist))
      # corr.data = sqrt(cbind(matrix(Q.C.dist,nrow = n),matrix(Q.bn.dist,nrow = n)))
      # corr.data = rmvnorm(sigma = matrix(c(1,-0.9,-0.9,1),ncol = 2), n = 500)
      # p.Brown and p.Fisher identical under negative correlation

      p.Brn.sqrtQ = empiricalBrownsMethod(data_matrix = t(corr.data),
                                          p_values = c(p.C,p.bn), extra_info = TRUE)

      print("dim corr.data")
      print(dim(corr.data))
      print("Brown's method correlation:")
      print(paste(cor(cbind(Q.C.dist,Q.bn.dist))))
      print("p.Fisher vs p.Brown:")
      print(paste(c(p.Fsr.bn,p.Brn.sqrtQ$P_test)))
    } else {
      p.Brn.sqrtQ = NA
    }


    # final results list
    ###########################
    re <- list('p.C' =      p.C,
               'p.D.Adj' =  p.bn,
               'p.Adj.Mthd' = S.D.adj$Test.Type,

               'p.C.resamp' =     p.C.dist,
               'p.D.Adj.resamp' = p.bn.dist,

               'Q.C' = S.C$Q,
               'Q.C' = S.D.adj$Q,
               'Q.C.resamp' =     Q.C.dist,
               'Q.D.Adj.resamp' = Q.bn.dist,

               'p.min' =      p.min.bn,
               'p.Fsr.rsmp' = p.Fsr.bn.rsmp,
               'p.Fsr.indp' = p.Fsr.bn,
               'p.Brown'=     p.Brn.sqrtQ,

               Test.Type = paste("method.C: ",method.C,
                                 ", method.D: ",method.D,
                                 ", method.Adj: ",method.bin, sep = ""),

               param = list("C" = S.C$param, "D" = S.D.adj$param)
    )
    return(re)

    # return(list('p.Cont' = p.C, 'pDcht.Adj' =  p.bn, 'p.Adj.Mthd' = S.D.adj$Test.Type,
    # 'p.min' = p.min.bn, 'p.Fsr.rsmp' = p.Fsr.bn.rsmp, 'p.Fsr.indp' = p.Fsr.bn, 'p.Brown'=p.Brn.sqrtQ))
  }




# Some extra functions to run SKAT on the null model and output Q.resampling
# in order to obtain a test of both hurdle model parameters under the global null.


vcHurdle_With_NullModel =
  # replaces SKAT::SKAT_With_NullModel() function. Calls modified functions below.
  function (Z, obj.res, kernel = "linear.weighted", method = "davies",
            weights.beta = c(1, 25), weights = NULL, impute.method = "fixed",
            r.corr = 0, is_check_genotype = TRUE, is_dosage = FALSE,
            missing_cutoff = 0.15, max_maf = 1, estimate_MAF = 1, SetID = NULL,
            out.z = NULL) {
    n <- dim(Z)[1]
    m <- dim(Z)[2]
    out.method <- SKAT:::SKAT_Check_Method(method, r.corr, m = m, n = n)
    method = out.method$method
    r.corr = out.method$r.corr
    IsMeta = out.method$IsMeta
    SKAT:::SKAT_Check_RCorr(kernel, r.corr)
    if (is.null(obj.res$n.all)) {
      obj.res$n.all = n
    }
    if (is.null(out.z)) {
      out.z <- SKAT:::SKAT_MAIN_Check_Z(Z, obj.res$n.all, obj.res$id_include,
                                        SetID, weights, weights.beta, impute.method, is_check_genotype,
                                        is_dosage, missing_cutoff, max_maf = max_maf, estimate_MAF = estimate_MAF)
    }
    if (out.z$return == 1) {
      out.z$param$n.marker <- m
      return(out.z)
    }
    if (length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1) {
      r.corr = 0
      method = "davies"
    }
    else if (length(r.corr) > 1 && sum(abs(out.z$Z.test - out.z$Z.test[,
                                                                       1])) == 0) {
      r.corr = 0
      method = "davies"
      msg <- sprintf("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!")
      warning(msg, call. = FALSE)
    }
    if (obj.res$out_type == "C") {
      if (kernel == "linear" || kernel == "linear.weighted") {
        re = vcHurdle.linear.Linear(obj.res$res, out.z$Z.test,
                                    obj.res$X1, kernel, out.z$weights, obj.res$s2,
                                    method, obj.res$res.out, obj.res$n.Resampling,
                                    r.corr = r.corr, IsMeta = IsMeta)
      }
      else {
        warning("Hurdle model not yet implemented for nonlinear Kernel") #TODO

        re = vcHurdle.linear.Other(obj.res$res, out.z$Z.test,
                                   obj.res$X1, kernel, out.z$weights, obj.res$s2,
                                   method, obj.res$res.out, obj.res$n.Resampling)
      }
    }
    else if (obj.res$out_type == "D") {
      if (kernel == "linear" || kernel == "linear.weighted") {
        re = SKAT:::SKAT.logistic.Linear(obj.res$res, out.z$Z.test,
                                         obj.res$X1, kernel, out.z$weights, obj.res$pi_1,
                                         method, obj.res$res.out, obj.res$n.Resampling,
                                         r.corr = r.corr, IsMeta = IsMeta)
      }
      else {
        re = SKAT:::SKAT.logistic.Other(obj.res$res, out.z$Z.test,
                                        obj.res$X1, kernel, out.z$weights, obj.res$pi_1,
                                        method, obj.res$res.out, obj.res$n.Resampling)
      }
    }
    re$param$n.marker <- m
    re$param$n.marker.test <- dim(out.z$Z.test)[2]
    return(re)
  }

vcHurdle.linear.Linear =
  function (res, Z, X1, kernel, weights = NULL, s2, method, res.out,
            n.Resampling, r.corr, IsMeta = FALSE) {

    if (length(r.corr) > 1 && dim(Z)[2] == 1) {
      r.corr = 0
    }
    if (IsMeta) {
      warning("Hurdle model not implemented for MetaSKAT.")
      re = SKAT:::SKAT_RunFrom_MetaSKAT(res = res, Z = Z, X1 = X1,
                                        kernel = kernel, weights = weights, s2 = s2, out_type = "C",
                                        method = method, res.out = res.out, n.Resampling = n.Resampling,
                                        r.corr = r.corr)
    }
    else if (length(r.corr) == 1) {
      re = vcHurdle.KMTest.linear.Linear(res, Z, X1, kernel, weights, #removed unnecessary r.corr argument
                                         s2, method, res.out, n.Resampling)
    }
    else {
      re = SKAT:::SKAT_Optimal_Linear(res, Z, X1, kernel, weights,
                                      s2, method, res.out, n.Resampling, r.corr)
    }
    return(re)
  }

vcHurdle.KMTest.linear.Linear =
  function (res, Z, X1, kernel, weights = NULL, s2, method, res.out,
            n.Resampling) {
    print("KMTest.linear.Linear.MG")
    print(paste("kernel",kernel))
    n = nrow(Z)
    m = ncol(Z)
    if (class(kernel) == "matrix") {
      K = kernel
    }
    else {
      K = vcHurdle.lskmTest.GetKernel(Z, kernel, weights, n, m) #weirdly this does not by default have the linear kernel option
    }
    Q = t(res) %*% K %*% res/(2 * s2)
    Q.res = NULL
    if (n.Resampling > 0) {
      Q.res <- rep(0, n.Resampling)
      for (i in 1:n.Resampling) {
        Q.res[i] = t(res.out[, i]) %*% K %*% res.out[, i]/(2 * s2)
      }
    }
    W = K - X1 %*% solve(t(X1) %*% X1) %*% (t(X1) %*% K)
    if (method == "davies") {
      W1 = W - (W %*% X1) %*% solve(t(X1) %*% X1) %*% t(X1)
    }
    if (method == "liu") {
      out <- SKAT:::Get_Liu_PVal(Q, W, Q.res)
      pval.zero.msg = NULL
    }
    else if (method == "liu.mod") {
      out <- SKAT:::Get_Liu_PVal.MOD(Q, W, Q.res)
      pval.zero.msg = NULL
    }
    else if (method == "davies") {
      out <- SKAT:::Get_Davies_PVal(Q, W1, Q.res)
      pval.zero.msg = out$pval.zero.msg
    }
    else {
      stop("Invalid Method!")
    }
    re <- list(p.value = out$p.value, p.value.resampling = out$p.value.resampling,
               Test.Type = method, Q = Q, Q.resampling = Q.res, param = out$param, pval.zero.msg = pval.zero.msg) #ADDED Q.resampling
    return(re)
  }

vcHurdle.lskmTest.GetKernel = # add linear kernel
  function (Z, kernel, weights, n, m) {
    if (kernel == "linear") { # for some reason this was missing before
      K = (Z %*% t(Z) )
    }
    if (kernel == "quadratic") {
      K = (Z %*% t(Z) + 1)^2
    }
    if (kernel == "IBS") {
      K = SKAT:::call_Kernel_IBS(Z, n, m)
    }
    if (kernel == "IBS.weighted") {
      K = SKAT:::call_Kernel_IBS_Weight(Z, n, m, weights)
    }
    if (kernel == "2wayIX") {
      K = SKAT:::call_Kernel_2wayIX(Z, n, m)
    }
    if (kernel == "IBS.weighted_OLD") {
      if (is.null(weights)) {
        qs = apply(Z, 2, mean)/(2)
        weights = 1/sqrt(qs)
      }
      else {
        weights <- weights^2
      }
      K1 = matrix(nrow = n, ncol = n)
      for (i in 1:n) {
        K1[i, ] = apply(abs(t(Z) - Z[i, ]) * weights, 2,
                        sum)
      }
      K = 1 - (K1)/(2 * sum(weights))
    }
    if (kernel == "IBS_OLD") {
      K1 = matrix(nrow = n, ncol = n)
      for (i in 1:n) {
        K1[i, ] = apply(abs(t(Z) - Z[i, ]), 2, sum)
      }
      K = (2 * m - K1)/(2 * m)
    }
    if (kernel == "2wayIX_OLD") {
      K = 1 + Z %*% t(Z)
      N1 = matrix(nrow = n, ncol = n)
      for (i in 1:n) {
        for (j in i:n) {
          N1[j, i] = N1[i, j] = SKAT:::K1_Help(Z[i, ], Z[j, ])
        }
      }
      K = K + N1
    }
    return(K)
  }


