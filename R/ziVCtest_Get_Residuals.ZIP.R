#############
# SKAT:::Get_SKAT_Residuals.linear

Get_ziVCtest_Residuals.ZIP =
  function (formula_1, formula_2, data,
            n.Resampling, type.Resampling, wgtDstr.Resampling,
            id_include,
            dist, link, asymp.stdz = TRUE)
  {

    # mod = lm(y.c ~ X, data = data)
    # X1.old <- model.matrix(formula, data = data)
    # X1.new <- SKAT:::Get_SKAT_Residuals.Get_X1(X1)
    # s2 = summary(mod)$sigma^2
    # res = mod$resid
    # n1 <- length(res)
    # res.out <- NULL

    formulaZI = as.formula(
      paste(
        as.character(formula_1)[[2]],
        paste(c(as.character(formula_1)[[3]],as.character(formula_2)[[3]]), collapse = " | "),
        sep = " ~ ")
    )

    # betas = ZImodel(formulaZI, data, wgts = 1, dist = dist, link = link)

    # formula = y.zi ~ X | X
    # ziVCmdl = zeroinflVC(formula = formulaZI, data = data, weights = weights, #subset, na.action, offset,
    #                      dist = c("poisson"), #, "negbin", "geometric"
    #                      link = c("logit"),  #, "probit", "cloglog", "cauchit", "log"
    #                      # control = zeroinfl.control(...),
    #                      model = TRUE, y = TRUE, x = FALSE)

    ziVCfn = function(wgts) {
      ZImodel =
        suppressWarnings( zeroinflVC(formula = formulaZI, data = data, weights = wgts, #subset, na.action, offset,
                 dist = dist, #, "negbin", "geometric"
                 link = link,  #, "probit", "cloglog", "cauchit", "log"
                 # control = zeroinfl.control(...),
                 model = TRUE, y = TRUE, x = FALSE)
        )
      return(ZImodel)
    }

    ziVCnullmdl = ziVCfn(wgts = 1)
    # names(ziVCnullmdl)
    res = ziVCnullmdl$score.resids
    # mapply(function(b1,b2,w) scores(X,Y.a,G,b1,b2,w), beta.pi.stars.a, beta.lda.stars.a, weights.mat )

    res.out = array(res, dim = c(dim(res),n.Resampling))

    data.boot = as.data.frame(data)
    n = nrow(data.boot)

    # cor(ziVCmdl$score.resids,data$Z)
    # cor(mod$residuals, data$Z)


    if (n.Resampling > 0) {
      if (type.Resampling == "permutation") {
        # res.out <- res %x% t(rep(1, n.Resampling))
        res.out = apply(res.out, c(2,3), sample)
        # dim(res.out)
        # res.out[1:5,1:2,1:2]
      } else if (type.Resampling == "bootstrap") {
        # res.out <- matrix(rnorm(n1 * n.Resampling, mean = 0,
        #                         sd = sqrt(s2)), ncol = n.Resampling)
        # X1_inv <- solve(t(X1) %*% X1)
        # res.out <- res.out - (X1 %*% X1_inv) %*% (t(X1) %*%
        #                                             res.out)

        for(b in 1:n.Resampling){
          # mapply(function(b1,b2,w) scores(X,Y.a,G,b1,b2,w), beta.pi.stars.a, beta.lda.stars.a, weights.mat )
            ZImodel = zeroinflVC(formula = formulaZI,
                                 data = data.boot[sample(1:n,replace = TRUE),],
                                 weights = 1, #subset, na.action, offset,
                                 dist = dist, #"poisson", "negbin", "geometric"
                                 link = link, #"logit", "probit", "cloglog", "cauchit", "log"
                                 # control = zeroinfl.control(...),
                                 model = TRUE, y = TRUE, x = FALSE)
            res.out[,,b] = ZImodel$score.resids
        }

      } else if ( grepl("perturbation", type.Resampling)) {
        ## for SKAT package: perturbing the residuals
        # res.out <- matrix(rnorm(n1 * n.Resampling, mean = 0,
        #                         sd = 1), ncol = n.Resampling)
        # res.out <- res.out * res
        # stop("Error: Perturbation is no more provided!")

        ## for ziVC test: perturbation resampling via refitting the model
        if( wgtDstr.Resampling == "exp"){
          weights.mat = as.data.frame(replicate(n.Resampling, rexp(n, rate=1)))
        } else if( wgtDstr.Resampling == "beta" ){
          weights.mat = as.data.frame(replicate(n.Resampling, 4*rbeta(n, 1/2,3/2 )))
        }

        res.out.list = lapply(as.data.frame(weights.mat), function(w) ziVCfn(w)$score.resids )
        res.out = array(unlist(res.out.list), c(n,2,n.Resampling))
        # all.equal(c(res.out.list[[11]]), c(res.out[,,11]))
        # dim(res.out)
        # res.out[1:5,1:2,1:10]; res[1,1]
        # hist(res.out[1,1,1:n.Resampling], xlim = c(-10,20), breaks = -20:40)

      } else {
        stop("Error: Wrong resampling method!")
      }
    }
    dimnames(res.out) = list(paste(id_include), c("count","zero"),paste("b",1:n.Resampling,sep = ""))

    if(asymp.stdz){
      res = n ^ -0.5 * res
      res.out = res.out / rep( sqrt(colSums(weights.mat)), each = dim(res.out)[1] * dim(res.out)[2] )
    }

    return(list(res = res, X1 = NULL, mdl = ziVCnullmdl,
                res.out = res.out, weights.mat = weights.mat, out_type = "ZIP",
                n.Resampling = n.Resampling, type.Resampling = type.Resampling, wgtDstr.Resampling = wgtDstr.Resampling,
                id_include = id_include))
  }
