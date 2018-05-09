
Get_ziVCtest_Residuals.Hurdle =
  function (formula_1, formula_2, data,
            n.Resampling, type.Resampling,
            id_include)
  {
    # for testing
    # kernel = "linear"
    # dosages = FALSE

    # From input to PaperHurdleSKAT
    # Y,X,G,
    # formula_null,
    # B = 1000
    # kernel = "linear", dosages = TRUE,
    # method = c('perturbation', 'bootstrap'),
    # method.bin = "Hybrid",
    # p.Brown = TRUE


    Y.name = strsplit( as.character(formula_1), "~") [[2]]
    Y = (data[Y.name])

    lin.pred.1.name = strsplit( as.character(formula_1), "~") [[3]]
    lin.pred.2.name = strsplit( as.character(formula_2), "~") [[3]]
    X.1.names = unlist(strsplit(as.character(lin.pred.1.name), " \\+ "))
    X.2.names = unlist(strsplit(as.character(lin.pred.2.name), " \\+ "))
    # X.names = unlist(strsplit(as.character(lin.pred.1.name), " \\+ "))

    # Binary Outcome and Covariates
    SubSet.C = as.matrix(as.data.frame(Y)) > 0
    Y.B = 1*SubSet.C
    colnames(Y.B) = "Y.B"
    # Continuous Outcome and Covariates given Y > 0
    Y.C = as.matrix(as.data.frame(Y))
    Y.C[!SubSet.C, ] = NA
    colnames(Y.C) = "Y.C"
    # Zero-inflated outcome
    Y.ZIP = as.matrix(Y)
    colnames(Y.ZIP) = "Y.ZIP"


    # cov.names = names(data.null)[grep('X.', names(data.null))]
    cov.names = union(X.1.names,X.2.names)
    if( all( cov.names %in% names(data)) ){
      data.X = data[cov.names]
    } else {
      data.X = as.data.frame(data)[cov.names]
    }
    # names(data.X)
    data.null = data.frame(Y.B=Y.B, Y.C=Y.C, data.X); names(data.null)

    formula.B = as.formula(paste("Y.B ~", paste(X.1.names, collapse  = " + ")))
    formula.C = as.formula(paste("Y.C ~", paste(X.2.names, collapse  = " + ")))
    # formula.ZIP = as.formula(paste("Y.ZIP ~", paste( paste(cov.names, collapse  = " + "), "|",
    #                                                  paste(cov.names, collapse  = " + ")
    # )))

    # From SKAT get residuals
    # mod = lm(y.c ~ X, data = data)
    # X1.old <- model.matrix(formula, data = data)
    # X1.new <- SKAT:::Get_SKAT_Residuals.Get_X1(X1)
    # s2 = summary(mod)$sigma^2
    # res = mod$resid
    # n1 <- length(res)
    # res.out <- NULL


    # Functions if starting from scratch:
    # binVCfn = function(wgts) {
    #   ZImodel =
    #     suppressWarnings(
    #       glm(formula = formula.B, data = data.null, weights = wgts,
    #                           family = binomial(link = "logit")
    #           )
    #     )
    #   ##
    #   ZImodel$score = "TBD"
    #   return(ZImodel)
    # }
    #
    # linVCfn = function(wgts) {
    #   ZImodel =
    #       lm(formula = formula.C, data = data.null, weights = wgts )
    #   ZImodel$score = "TBD"
    #   return(ZImodel)
    # }


    # ziVCnullmdl = ziVCfn(wgts = 1)
    # res = ziVCnullmdl$score.resids
    # mapply(function(b1,b2,w) scores(X,Y.a,G,b1,b2,w), beta.pi.stars.a, beta.lda.stars.a, weights.mat )
    # res.out = array(res, dim = c(dim(res),n.Resampling))

    # data.boot = as.data.frame(data)
    # n = nrow(data.boot)

    #
    #     if (n.Resampling > 0) {
    #       if (type.Resampling == "permutation") {
    #         # res.out <- res %x% t(rep(1, n.Resampling))
    #         res.out = apply(res.out, c(2,3), sample)
    #         # dim(res.out)
    #         # res.out[1:5,1:2,1:2]
    #       }
    #       else if (type.Resampling == "bootstrap") {
    #         # res.out <- matrix(rnorm(n1 * n.Resampling, mean = 0,
    #         #                         sd = sqrt(s2)), ncol = n.Resampling)
    #         # X1_inv <- solve(t(X1) %*% X1)
    #         # res.out <- res.out - (X1 %*% X1_inv) %*% (t(X1) %*%
    #         #                                             res.out)
    #
    #         for(b in 1:n.Resampling){
    #           # mapply(function(b1,b2,w) scores(X,Y.a,G,b1,b2,w), beta.pi.stars.a, beta.lda.stars.a, weights.mat )
    #           ZImodel = zeroinflVC(formula = formulaZI,
    #                                data = data.boot[sample(1:n,replace = TRUE),],
    #                                weights = 1, #subset, na.action, offset,
    #                                dist = dist, #"poisson", "negbin", "geometric"
    #                                link = link, #"logit", "probit", "cloglog", "cauchit", "log"
    #                                # control = zeroinfl.control(...),
    #                                model = TRUE, y = TRUE, x = FALSE)
    #           res.out[,,b] = ZImodel$score.resids
    #         }
    #
    #       }
    #       else if ( grepl("perturbation", type.Resampling)) {
    #         ## for SKAT package: perturbing the residuals
    #         # res.out <- matrix(rnorm(n1 * n.Resampling, mean = 0,
    #         #                         sd = 1), ncol = n.Resampling)
    #         # res.out <- res.out * res
    #         # stop("Error: Perturbation is no more provided!")
    #
    #         ## for ziVC test: perturbation resampling via refitting the model
    #         if(grepl("exp", type.Resampling)){
    #           weights.mat = as.data.frame(replicate(n.Resampling, rexp(n, rate=1)))
    #         } else if(grepl("beta", type.Resampling)){
    #           weights.mat = as.data.frame(replicate(n.Resampling, 4*rbeta(n, 1/2,3/2 )))
    #         }
    #
    #         res.out.list = lapply(as.data.frame(weights.mat), function(w) ziVCfn(w)$score.resids )
    #         res.out = array(unlist(res.out.list), c(n,2,n.Resampling))
    #         # all.equal(c(res.out.list[[11]]), c(res.out[,,11]))
    #         # dim(res.out)
    #         # res.out[1:5,1:2,1:10]; res[1,1]
    #         # hist(res.out[1,1,1:n.Resampling], xlim = c(-10,20), breaks = -20:40)
    #
    #       }
    #       else {
    #         stop("Error: Wrong resampling method!")
    #       }
    #     }
    #     dimnames(res.out) = list(paste(id_include), c("count","zero"),paste("b",1:n.Resampling,sep = ""))
    #
    #
    # if(asymp.stdz){
    #   res = n ^ -0.5 * res
    #   res.out = res.out / rep( sqrt(colSums(weights.mat)), each = dim(res.out)[1] * dim(res.out)[2] )
    # }


    if( grepl('perturbation', type.Resampling) ){
      # Dichotomous Trait: out_type = "D"
      # browser()
      obj.D = SKAT_Null_Model.ziHurdle(formula.B, data = data.null, out_type="D",
                                       n.Resampling = n.Resampling, type.Resampling = type.Resampling,
                                       wgts.Resampling = NULL,
                                       Adjustment = TRUE)


      wgts.mat = obj.D$wgts.Resampling
      # Continuous Trait: out_type = "C"
      obj.C = SKAT_Null_Model.ziHurdle(formula.C, data = data.null, out_type="C",
                                       n.Resampling = n.Resampling, type.Resampling = type.Resampling,
                                       wgts.Resampling = wgts.mat,
                                       Adjustment = FALSE)

      # w.C = wgts.mat[match(names(res), rownames(obj)),] where obj and res are internal variables below
      # identical(w.C[,2],obj.C$wgts.Resampling[,2])

      print("names obj.D")
      print(names(obj.D))

      print("names obj.C")
      print(names(obj.C))
      # browser()

      print("Class obj.D on creation")
      print(class(obj.D))
      print("Class obj.C on creation")
      print(class(obj.C))

    } else {
      print("Resampling method not implimented for Hurdle model VC test.")
    }

    return(list(
      obj.D = obj.D,
      obj.C = obj.C,
      res.D = obj.D$res, res.C = obj.C$res,
      mdl.D = obj.D$mod, mdl.C = obj.C$mod,
      X1.D = obj.D$X1, X1.C = obj.C$X1,
      res.out.D = obj.D$res.out, res.out.C = obj.C$res.out,
      weights.mat = wgts.mat, out_type = "Hurdle",
      n.Resampling = n.Resampling, type.Resampling = type.Resampling,
      id_include = id_include))
  }


#########################

SKAT_Null_Model.ziHurdle = function (formula, data = NULL, out_type = "C",
                                     n.Resampling = 0,
                                     type.Resampling = "perturbation",
                                     wgts.Resampling = NULL,
                                     Adjustment = TRUE) {
  print("SKAT_Null_Model.ziHurdle")
  # print(paste("type.Resampling", type.Resampling))
  # print(paste("n.Resampling", n.Resampling))

  SKAT_MAIN_Check_OutType.ziHurdle(out_type)
  obj1 <- model.frame(formula, na.action = na.omit, data)
  obj2 <- model.frame(formula, na.action = na.pass, data)
  n <- dim(obj2)[1]
  n1 <- dim(obj1)[1]
  id_include <- SKAT:::SKAT_Null_Model_Get_Includes(obj1, obj2)
  n1 <- length(id_include)
  if (n1 < 2000 && out_type == "D" && Adjustment) {
    MSG <- sprintf("Sample size (non-missing y and X) = %d, which is < 2000. The small sample adjustment is applied!\n",
                   n)
    cat(MSG)
    n.Resampling.kurtosis = 10000
    re <- SKAT_Null_Model_MomentAdjust.ziHurdle(formula, data,
                                                n.Resampling = n.Resampling,
                                                type.Resampling = type.Resampling,
                                                is_kurtosis_adj = TRUE,
                                                n.Resampling.kurtosis = n.Resampling.kurtosis)
    return(re)
  }
  if (n - n1 > 0) {
    MSG <- sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                   n - n1)
    warning(MSG, call. = FALSE)
  }
  if (out_type == "C") {
    re <- Get_SKAT_Residuals.linear.ziHurdle(formula, data,
                                             n.Resampling, type.Resampling, wgts.Resampling,
                                             id_include)

  } else if(out_type == "D") {
    re <- Get_SKAT_Residuals.logistic.ziHurdle(formula, data,
                                               n.Resampling, type.Resampling, wgts.Resampling,
                                               id_include)
  }
  class(re) <- "SKAT_NULL_Model"
  re$n.all <- n
  return(re)
}
# SKAT::SKAT_Null_Model
# SKAT:::SKAT_MAIN_Check_OutType
# assignInNamespace("SKAT_Null_Model", value = SKAT_Null_Model.ziHurdle, ns ="SKAT") # envir = "package:SKAT"
# reassignInPackage(name = "SKAT_Null_Model", pkgName = "SKAT", value = SKAT_Null_Model.ziHurdle, keepOld = FALSE)
# rm(SKAT_Null_Model, envir = .GlobalEnv)
# environment(SKAT_Null_Model)
# For testing
# obj.C = SKAT_Null_Model(Y.C ~ X.C, out_type="C", n.Resampling=500)
# S.C = SKAT(G.C, obj.C, kernel = kernel, is_dosage = dosages)

SKAT_MAIN_Check_OutType.ziHurdle =
  function(out_type)
  {
    if (out_type != "C" && out_type != "D" && out_type != "ZIP") {
      stop("Invalid out_type!. Please use either
           \"C\" for the continous outcome or
           \"D\" for the dichotomous outcome or
           \"ZIP\" for the zero-inflated poisson outcome.")
    }
  }


#########################
SKAT_Null_Model_MomentAdjust.ziHurdle =
  function (formula, data = NULL,
            n.Resampling = 0, type.Resampling = "perturbation", wgts.Resampling=NULL,
            is_kurtosis_adj = TRUE, n.Resampling.kurtosis = 10000){
    print("SKAT_Null_Model_MomentAdjust.ziHurdle")

    obj1 <- model.frame(formula, na.action = na.omit, data)
    obj2 <- model.frame(formula, na.action = na.pass, data)
    n <- dim(obj2)[1]
    n1 <- dim(obj1)[1]
    id_include <- SKAT:::SKAT_Null_Model_Get_Includes(obj1, obj2)
    if (n - n1 > 0) {
      MSG <- sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                     n - n1)
      warning(MSG, call. = FALSE)
    }
    re1 <- Get_SKAT_Residuals.logistic.ziHurdle(formula, data,
                                                n.Resampling, type.Resampling, wgts.Resampling=NULL,
                                                id_include)
    re1$n.all <- n
    re2 <- NULL
    if (is_kurtosis_adj == TRUE) {
      re2 <- Get_SKAT_Residuals.logistic.ziHurdle(formula, data,
                                                  n.Resampling.kurtosis, type.Resampling, wgts.Resampling=NULL,
                                                  id_include)
    }
    class(re1) <- "SKAT_NULL_Model"
    re <- list(re1 = re1, re2 = re2, is_kurtosis_adj = is_kurtosis_adj,
               type = "binary")
    class(re) <- "SKAT_NULL_Model_ADJ"
    return(re)
  }

#########################
Get_SKAT_Residuals.logistic.ziHurdle =
  function (formula, data,
            n.Resampling, type.Resampling, wgts.Resampling=NULL,
            id_include=NA) {
    print("Get_SKAT_Residuals.logistic.ziHurdle")
    # print(paste("n.Resampling", n.Resampling))

    obj <- model.frame(formula, na.action = na.pass, data)
    n <- dim(obj)[1]

    # mod = lm(formula, data)
    X1 <- model.matrix(formula, data = data)
    X1 <- SKAT:::Get_SKAT_Residuals.Get_X1(X1)

    glmfit = glm(formula, data = data, family = "binomial")
    betas = glmfit$coef
    mu = glmfit$fitted.values
    eta = glmfit$linear.predictors
    n.case = sum(glmfit$y)
    pi_1 = mu * (1 - mu)
    res = glmfit$y - mu
    n1 <- length(res)



    res.out <- NULL
    if (n.Resampling > 0) {
      if(type.Resampling == "perturbation") {
        # reinstated depricated perturbation method in order to unify pertrubation weights between
        # Continuous and Dichotomous Analyses
        # not really necessary except to account for covariates.
        # (p-values indpendent between binary and continuous models with intercept-only models)
        if(!is.null(wgts.Resampling[1])){
          wgt.res = wgts.Resampling
        } else {
          wgt.res = matrix(rnorm(n * n.Resampling, mean = 0,
                                 sd = 1), ncol = n.Resampling)
        }
        wgt.res = wgt.res[match(names(res), rownames(obj)),]
        # print("head perturbation weights:")
        # print(cbind(names(res)[1:10],head(wgt.res[,1],10)))
        res.out <- wgt.res * res
      } else {
        if (is.null(res.out)) {
          stop("Hurdle model only implemented with perturbation resampling!")
        }
      }
    }
    return(list(res = res, X1 = X1, res.out = res.out, out_type = "D", mod = glmfit,
                n.Resampling = n.Resampling, type.Resampling = type.Resampling, wgts.Resampling = wgt.res,
                id_include = id_include, mu = mu, pi_1 = pi_1))
  }


#########################
Get_SKAT_Residuals.linear.ziHurdle =
  function (formula, data,
            n.Resampling, type.Resampling, wgts.Resampling=NULL,
            id_include=NA){
    print("Get_SKAT_Residuals.linear.ziHurdle")
    # print(paste("n.Resampling", n.Resampling))

    mod = lm(formula, data = data)
    X1 <- model.matrix(formula, data = data)
    X1 <- SKAT:::Get_SKAT_Residuals.Get_X1(X1)
    s2 = summary(mod)$sigma^2
    res = mod$resid
    n1 <- length(res)

    obj <- model.frame(formula, na.action = na.pass, data)
    n <- dim(obj)[1]

    res.out <- NULL
    if (n.Resampling > 0) {
      if (type.Resampling == "perturbation") {
        if(!is.null(wgts.Resampling[1])){
          wgt.res = wgts.Resampling
        } else {
          wgt.res = matrix(rnorm(n * n.Resampling, mean = 0,
                                 sd = 1), ncol = n.Resampling)
        }
        wgt.res = wgt.res[match(names(res), rownames(obj)),]

        #reinstated depricated perturbation method: matches subject-specific perturbation weights
        res.out <- wgt.res* res
        warning("Perturbation method is not recommended!")
      }
      else {
        stop("Error: Wrong resampling method!")
      }
    }
    return(list(res = res, X1 = X1, res.out = res.out, out_type = "C", mod = mod,
                n.Resampling = n.Resampling, type.Resampling = type.Resampling, wgts.Resampling = wgt.res,
                id_include = id_include, s2 = s2))
  }



