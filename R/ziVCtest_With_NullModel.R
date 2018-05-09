ziVCtest_With_NullModel = function (Z, obj.res, kernel = "linear.weighted", method = "davies",
                                    weights.beta = c(1, 25), weights = NULL, impute.method = "fixed",
                                    r.corr = 0, is_check_genotype = TRUE, is_dosage = FALSE,
                                    missing_cutoff = 0.15, max_maf = 1, estimate_MAF = 1, SetID = NULL,
                                    out.z = NULL) {
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  out.method <- VC_Check_Method(method, r.corr, m = m, n = n)
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
  else if (length(r.corr) > 1 && sum(abs(out.z$Z.test - out.z$Z.test[, 1])) == 0) {
    r.corr = 0
    method = "davies"
    msg <- sprintf("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!")
    warning(msg, call. = FALSE)
  }
  if (obj.res$out_type == "C") {
    if (kernel == "linear" || kernel == "linear.weighted") {
      print("entering SKAT.linear.Linear")

      re = SKAT:::SKAT.linear.Linear(obj.res$res, out.z$Z.test,
                                     obj.res$X1, kernel, out.z$weights, obj.res$s2,
                                     method, obj.res$res.out, obj.res$n.Resampling,
                                     r.corr = r.corr, IsMeta = IsMeta)
    } else {
      re = SKAT:::SKAT.linear.Other(obj.res$res, out.z$Z.test,
                                    obj.res$X1, kernel, out.z$weights, obj.res$s2,
                                    method, obj.res$res.out, obj.res$n.Resampling)
    }
  } else if (obj.res$out_type == "D") {
    if (kernel == "linear" || kernel == "linear.weighted") {
      print("entering SKAT.logistic.Linear")

      re = SKAT:::SKAT.logistic.Linear(obj.res$res, out.z$Z.test,
                                       obj.res$X1, kernel, out.z$weights, obj.res$pi_1,
                                       method, obj.res$res.out, obj.res$n.Resampling,
                                       r.corr = r.corr, IsMeta = IsMeta)
    } else {
      re = SKAT:::SKAT.logistic.Other(obj.res$res, out.z$Z.test,
                                      obj.res$X1, kernel, out.z$weights, obj.res$pi_1,
                                      method, obj.res$res.out, obj.res$n.Resampling)
    }
  } else if (obj.res$out_type == "ZIP") {
    if (kernel == "linear" || kernel == "linear.weighted") {
      print("using ziVCtest.Poisson.Linear")

      re = ziVCtest.Poisson.Linear(res = obj.res$res, Z = out.z$Z.test, X1 = obj.res$X1,
                                   kernel = kernel, weights = out.z$weights,
                                   method = method, res.out = obj.res$res.out,
                                   n.Resampling = obj.res$n.Resampling,
                                   r.corr = r.corr, IsMeta = IsMeta)
    } else {
      re = ziVCtest.Poisson.Other(res = obj.res$res, Z = out.z$Z.test, X1 = obj.res$X1,
                                   kernel = kernel, weights = out.z$weights,
                                   method = method, res.out = obj.res$res.out,
                                   n.Resampling = obj.res$n.Resampling,
                                   r.corr = r.corr, IsMeta = IsMeta)
      # re = ziVCtest.Poisson.Other(obj.res, out.z$Z.test,
                                  # kernel, out.z$weights,
                                  # method, obj.res$res.out, obj.res$n.Resampling)
    }
  } else if (obj.res$out_type == "ZINB") {
    stop("ZINB not yet implemented.") #TODO

  } else if (obj.res$out_type == "Hurdle") {
    if (kernel == "linear" || kernel == "linear.weighted") {
      print("using ziVCtest.Hurdle.Linear")
      re = ziVCtest.Hurdle.Linear(obj.res, out.z$Z.test,
                                  kernel, out.z$weights,
                                  dosages = is_dosage,
                                  method,
                                  r.corr = r.corr, IsMeta = IsMeta)
    } else {
      stop("Hurdle not yet implemented for nonlinear Kernel") #TODO
      re = ziVCtest.Hurdle.other(obj.res, out.z$Z.test,
                                 kernel, out.z$weights,
                                 dosages = is_dosage,
                                 method)
    }
  }
  re$param$n.marker <- m
  re$param$n.marker.test <- dim(out.z$Z.test)[2]
  return(re)
}


# SKAT:::SKAT_Check_Method
VC_Check_Method = function(method, r.corr, n = NULL, m = NULL)
{
  IsMeta = FALSE
  if (
    # TODO: clean up methods for hurdle and ZI models
    method != "default.hurdle" &&
    method != "davies.default" && method != "davies.precise" &&
    method != "imhof" && method != "farebrother" &&
    method != "liu" && method != "davies" &&
    method != "liu.mod" &&
    method != "optimal" && method != "optimal.moment" &&
    method != "optimal.mod" && method != "adjust" && method !=
    "optimal.adj" && method != "optimal.moment.adj" && method !=
    "SKAT" && method != "SKATO" && method != "Burden" &&
    method != "SKATO.m" && method != "davies.M" && method !=
    "optimal.adj.M") {
    stop("Invalid method!")
  }
  if (method == "davies.M") {
    IsMeta = TRUE
    method = "davies"
  }
  else if (method == "optimal.adj.M") {
    IsMeta = TRUE
    method = "optimal.adj"
  }
  if (method == "SKAT") {
    method = "davies"
    r.corr = 0
  }
  else if (method == "SKATO") {
    method = "optimal.adj"
  }
  else if (method == "SKATO.m") {
    method = "optimal.moment.adj"
  }
  else if (method == "Burden") {
    method = "davies"
    r.corr = 1
  }
  if ((method == "optimal" || method == "optimal.moment") &&
      length(r.corr) == 1) {
    r.corr = (0:10)/10
  }
  else if ((method == "optimal.mod" || method == "optimal.adj" ||
            method == "optimal.moment.adj") && length(r.corr) == 1) {
    r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
  }
  if (method == "optimal") {
    method = "davies"
  }
  else if (method == "optimal.moment") {
    method = "liu.mod"
  }
  if (!is.null(n) && !is.null(m) && length(r.corr) > 1) {
    if (m/n < 1 && n > 5000) {
      IsMeta = TRUE
    }
  }
  re <- list(method = method, r.corr = r.corr, IsMeta = IsMeta)
  return(re)
}
