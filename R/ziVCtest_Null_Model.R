#' Function for fitting the Null Model (under assumption of no genetic effect).
#'
#' Compute model parameters and residuals. Obtain distribution of residuals via resampling
#' for downstream SNP-set testing and p-values.
#'
#' @param formula a formula object, a symbolic description of the Null Model to be fitted.
#' The response is a zero-inflated variable. If desired, separate models can be set for the zero and mean responses instead (see below).
#' @param formula_1 a formula object describing the zero model. Set both \code{formula_1} and \code{formula_2}, otherwise, specifying \code{formula} will set both models equal.
#' @param formula_2 a formula object describing the mean model. Set both \code{formula_1} and \code{formula_2}, otherwise, specifying \code{formula} will set both models equal.
#' @param data a data frame containing the variables in the model.
#' @param out_type the outcome type. options are "ZIP" (Zero Inflated Poisson), "ZINB" (Zero Inflated Negative Binomial), "Hurdle" (Hurdle Model for zero-inflated continuous outcome).
#' @param n.Resampling a numeric value indicating the number of resampling iterations to use for estimating the distribution of the residuals under the Null Model (default=1000).
#' @param type.Resampling resampling methods "perturbation", "bootstrap", or "permutation" (default="perturbation"). Not all options are available for all models.
#' @param wgtDstr.Resampling the distribution of the resampling weights if the \code{type.Resampling} is set to "perturbation". Options: exponential ("exp), beta ("beta"), normal ("normal").
#'
#'
#' @return This function returns an object that has model parameters and residuals for the Null Model of no association between genetic variables and outcome phenotypes. The resampled distribution of residuals under the null hence only needs to be calculated once and reused for each SNP-set. Association testing is then carried out using the \code{VCtest} function.
#'
#' @seealso \code{\link{VCtest}} \code{\link[SKAT]{SKAT_Null_Model}}
#' @examples
#' #Modified SKAT example
#' data(SKAT.example)
#' SKAT.example.ZI = SKAT.example
#'
#' set.seed(123)
#' SKAT.example.ZI$y.zi = ifelse(y.b, rpois(length(SKAT.example$y.b), round(exp(SKAT.example$y.c),0)), 0)
#' hist(SKAT.example.ZI$y.zi, xlim = c(0,20), breaks = 0:60, right = FALSE)
#'
#' ziVCtest_null = ziVCtest_Null_Model(formula = y.zi ~ X,
#'                                     data = SKAT.example.ZI,
#'                                     out_type = "ZIP",
#'                                     n.Resampling = 500,
#'                                     type.Resampling = "perturbation",
#'                                     wgtDstr.Resampling = "exp"
#' )
#'
#' result = VCtest(Z=SKAT.example.ZI$Z, obj=ziVCtest_null,
#'                 kernel = "linear", method = "davies.default")
#'
#' @export
ziVCtest_Null_Model = function (formula, formula_1=NULL, formula_2=NULL,
                                data = NULL, out_type = "ZIP",
                                n.Resampling = 1000,
                                type.Resampling = c("perturbation",
                                                    "bootstrap",
                                                    "permutation"),
                                wgtDstr.Resampling = c("exp","beta", "normal"))
{
  # for testing:
  # data = SKAT.example.ZIP
  # formula = y.zi~X; formula_1=NULL; formula_2=NULL
  # SKAT:::SKAT_MAIN_Check_OutType(out_type)

  type.Resampling <- match.arg(type.Resampling)
  wgtDstr.Resampling <- match.arg(wgtDstr.Resampling)

  ziVCtest_MAIN_Check_OutType = function (out_type) {
    if (out_type != "ZIP" && out_type != "Hurdle") {
      stop("Invalid out_type!. Please use either \"ZIP\" or \"Hurdle\" to model the zero-inflated outcome.")
    }
  }
  ziVCtest_MAIN_Check_OutType(out_type)

  # set formula:
  if( !is.null(formula) ){
    if ( !is.null(formula_1) || !is.null(formula_2) ) {
      warning( "formula_1 (zeros) and formula_2 (mean) are ignored when \"formula\" is specified ")
    }
    formula_1 = formula_2 = formula
    obj1 <- model.frame(formula, na.action = na.omit, data)
    obj2 <- model.frame(formula, na.action = na.pass, data)

  } else {
    if ( is.null(formula_1) || is.null(formula_2) ) {
      stop( "When \"formula\" not specified, must specify both \"formula_1\" (zeros) and \"formula_2\" (mean)!")
    }
    # for testing:
    # formula_1 = y.zi ~ X
    # formula_2 = y.zi ~ W+X^2
    formula_12 = as.formula(
      paste(
        as.character(formula_1)[[2]],
            paste(c(as.character(formula_1)[[3]],as.character(formula_2)[[3]]), collapse = " + "),
        sep = " ~ ")
    )
    obj1 <- model.frame(formula_12, na.action = na.omit, data)
    obj2 <- model.frame(formula_12, na.action = na.pass, data)

  }

  n <- dim(obj2)[1]
  n1 <- dim(obj1)[1]
  id_include <- SKAT:::SKAT_Null_Model_Get_Includes(obj1, obj2)
  n1 <- length(id_include)

  if (n1 < 2000 && out_type == "D" && Adjustment) {
    MSG <- sprintf("Sample size (non-missing y and X) = %d, which is < 2000. The small sample adjustment is applied!\n",
                   n)
    cat(MSG)
    n.Resampling.kurtosis = 10000
    re <- SKAT_Null_Model_MomentAdjust(formula, data, n.Resampling,
                                       type.Resampling = type.Resampling, is_kurtosis_adj = TRUE,
                                       n.Resampling.kurtosis = n.Resampling.kurtosis)
    return(re)
  }

  if (n - n1 > 0) {
    MSG <- sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                   n - n1)
    warning(MSG, call. = FALSE)
  }

  if (out_type == "C") {
    re <- SKAT:::Get_SKAT_Residuals.linear(formula, data, n.Resampling,
                                    type.Resampling, id_include)
  } else if (out_type == "D") {
    re <- SKAT:::Get_SKAT_Residuals.logistic(formula, data, n.Resampling,
                                      type.Resampling, id_include)

  } else if (out_type == "ZIP") {
    if(wgtDstr.Resampling == "normal"){
      warning("\"wgtDstr.Resampling\" must have mean and variance of 1 for ZIP model: choose \"exp\" or \"beta\".")}
print(wgtDstr.Resampling)
    re <- Get_ziVCtest_Residuals.ZIP(formula_1, formula_2, data,
                                     n.Resampling, type.Resampling, wgtDstr.Resampling = wgtDstr.Resampling,
                                     id_include,
                                     dist = "poisson", link = "logit", asymp.stdz = TRUE)

  } else if (out_type == "ZINB") {
    if(wgtDstr.Resampling == "normal"){
      warning("\"wgtDstr.Resampling\" must have mean and variance of 1 for ZINB model: choose \"exp\" or \"beta\".")}
    # TODO check that dist = "negbin" works properly
    re <- Get_ziVCtest_Residuals.ZIP(formula_1, formula_2, data,
                                     n.Resampling, type.Resampling, wgtDstr.Resampling = wgtDstr.Resampling,
                                     id_include,
                                     dist = "negbin", link = "logit", asymp.stdz = TRUE)

  } else if (out_type == "Hurdle") {
    warning("out_type \"Hurdle\" code is \"beta\" use at own risk.")
    if(wgtDstr.Resampling != "normal"){
      warning("\"wgtDstr.Resampling\" is automatically set to \"normal\" for Hurdle model")}
    re <- Get_ziVCtest_Residuals.Hurdle(formula_1, formula_2, data,
                                        n.Resampling, type.Resampling,
                                        id_include)
  }

  class(re) <- "ziVC_NULL_Model"
  re$n.all <- n
  return(re)
}


