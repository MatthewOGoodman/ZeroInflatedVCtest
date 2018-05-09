#pval.zero.msg
#(only when p.value=0) text message that shows how small the p.value is. ex. "Pvalue < 1.000000e-60" when the p.value is smaller than 10^{-60}
#param
#estimated parameters of each method.
#param$Is_Converged
#(only with method="davies") an indicator of the convergence (1=convergence, 0=non-convergence). When 0 (not converged), "liu" method will be used to compute p-values.
#param$n.marker
#a number of SNPs in the genotype matrix.
#param$n.marker.test
#a number of SNPs used for the test. It can be different from param$n.marker when some markers are monomorphic or have higher missing rates than the missing_cutoff.


#' Test association between a SNP-set
#' and a zero-inflated phenotype using kernel regression.
#'
#' Use pre-computed null model to test association between a zero-inflated phenotype and a given set of SNPs/genes.
#'
#' @param Z a numeric genotype matrix with each row as a different individual
#'   and each column as a separate gene/snp. Each genotype should be coded as 0,
#'   1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele
#'   and a is a minor allele. Missing genotypes will be imputed by the simple
#'   Hardy-Weinberg equilibrium (HWE) based imputation.
#' @param obj an output object of the ziVCtest_Null_Model function.
#' @param kernel a type of kernel (default= "linear.weighted"). See detail
#'   section.
#' @param method a method to compute the p-value (default= "davies.default").
#' using the quadratic form of the test statistic (using R package CompQuadForm)
#'   "davies.default" represents an exact method that computes the p-value by inverting the
#'   characteristic function of the quadratic form,
#'   "davies.precise" includes precision parameters lim = 5e7, acc = 1e-12 (see CompQuadForm::davies)
#'   "farebrother" represents Farebrother's approximation as a single ChiSquare distribution (Algorithm AS204)
#'   "imhof" represents Imhof's approximation using an expansion of the distribution function,
#'   "liu" represents an approximation method that matches the first 3 moments
#'   )
#' @param weights.beta a numeric vector of parameters for the beta weights for
#'   the weighted kernels. If you want to use your own weights, please use the
#'   "weights" parameter. It will be ignored if "weights" parameter is not null.
#' @param weights a numeric vector of weights for the weighted kernels. It is
#'   \eqn{\sqrt{w}} in the SKAT paper. So if you want to use the Madsen and Browning
#'   (2009) weight, you should set each element of weights as \eqn{1/\sqrt{p(1-p)}}, not
#'   1/ p(1-p). When it is NULL, the beta weight with the "weights.beta"
#'   parameter is used.
#' @param impute.method a method to impute missing genotypes (default= "fixed").
#'   "bestguess" imputes missing genotypes as most likely values (0,1,2),
#'   "random" imputes missing genotypes by generating binomial(2,p) random
#'   variables (p is the MAF), and "fixed" imputes missing genotypes
#'   by assigning the mean genotype values (2p).
#' @param r.corr the \eqn{\rho} parameter for the compound symmetric correlation
#'   structure kernels (default= 0). If you give a vector value, SKAT will
#'   conduct the optimal test. It will be ignored if method="optimal" or
#'   method="optimal.adj". See details.
#' @param is_check_genotype a logical value indicating whether to check the
#'   validity of the genotype matrix Z (default= TRUE). If Z has non-SNP data,
#'   please set it FALSE, otherwise you will get an error message. If it is
#'   FALSE and you use weighted kernels, the weights should be given through the
#'   "weights" parameter.
#' @param is_dosage a logical value indicating whether the matrix Z is a dosage
#'   matrix. If it is TRUE, SKAT will ignore "is_check_genotype".
#' @param missing_cutoff a cutoff of the missing rates of SNPs (default=0.15).
#'   Any SNPs with missing rates higher than the cutoff will be excluded from
#'   the analysis.
#' @param max_maf a cutoff of the maximum minor allele frequencies (MAF)
#'   (default=1, no cutoff). Any SNPs with MAF > cutoff will be excluded from
#'   the analysis.
#' @param estimate_MAF a numeric value indicating how to estimate MAFs for the
#'   weight calculation and the missing genotype imputation. If estimate_MAF=1
#'   (default), SKAT uses all samples to estimate MAFs. If estimate_MAF=2, only
#'   samples with non-missing phenotypes and covariates are used to estimate
#'   MAFs.
#' @param SSD.INFO an SSD_INFO object returned from Open_SSD.
#' @param SetID a character value of Set ID. A set ID of each set can be found
#'   from SetInfo object in SSD.INFO.
#' @param SetIndex a numeric value of Set index. A set index of each set can be
#'   found from SetInfo object in SSD.INFO.
#' @param obj.SNPWeight an output object of Read_SNP_WeightFile (default=NULL).
#'   If NULL, the beta weight with the "weights.beta" parameter will be used.
#'
#' @return A list containing values from among the following (depending on
#'   class(obj), where obj is the output of the Null_Model function)
#' \item{p.val.pi}{ p-value from a test of null genetic effect on the excess
#'   zero parameter in the ZIP or ZINB model}
#' \item{p.val.ld}{  p-value from a test of null genetic effect on the mean
#'   parameter in the ZIP or ZINB model}
#' \item{p.val.cmb}{  p-value for the inverse-variance combined test of the
#'   global null for both the pi and lambda parameters}
#' \item{p.val.pi.resamp}{  resampled p-values for the pi parameter}
#' \item{p.val.ld.resamp}{  resampled p-values for the lambda parameter}
#' \item{p.val.cmb.resamp}{  resampled p-values for the combined test}
#' \item{p.val.Fsr}{  p-value for the Fisher method using the resampled p-value
#'   distributions}
#' \item{p.val.min}{  p-value for the min-p method using the resampled p-value
#'   distributions}
#' \item{Test.Type}{  parameter indicating the Test.Type input}
#' \item{Q.pi}{  the test statistic for pi}
#' \item{Q.ld}{  the test statistic for lambda}
#' \item{Q.resampling.pi}{  the resampled test statistics for pi}
#' \item{Q.resampling.ld}{  the resampled test statistics for lambda}
#' \item{param}{  some parameters from the CompQuadForm calculation which obtains
#'   the p-value from the distribution of the quadratic form.}
#' \item{p.C}{  p.value test of null genetic effect on mean parameter of non-zero
#'   continuous outcome in hurdle model}
#' \item{p.D.Adj}{  p.value for test of null genetic effect on dichotomized zero
#'   vs non-zero outcome in hurdle model}
#' \item{p.Fsr.rsmp}{  p-value from Fisher method of combining p.C and p.D, using
#'   resampling}
#' \item{p.Fsr.indp}{  p-value from Fisher method of combining p.C and p.D,
#'   assuming independence}
#' \item{p.Brown.P_test}{  p-value from Brown's method of combining p.C and p.D
#'   using correlation of test statistics}
#' \item{p.Brown.P_Fisher}{  Fisher's p-value under independence as caluculated
#'   by Brown's method software}
#' \item{p.Brown.Scale_Factor_C}{  correction factor in Brown's method}
#' \item{p.Brown.DF}{  degrees of freedom after correcction in Brown's method}
#'
#'
#' @seealso \code{\link{ziVCtest_Null_Model}}; \code{\link[SKAT]{SKAT}}
#'
#' @export
VCtest = function (Z, obj, kernel = "linear.weighted", method = "davies",
          weights.beta = c(1, 25), weights = NULL, impute.method = "fixed",
          r.corr = 0, is_check_genotype = TRUE, is_dosage = FALSE,
          missing_cutoff = 0.15, max_maf = 1, estimate_MAF = 1) {

  if (kernel != "linear" && kernel != "linear.weighted") {
    if (class(obj) == "SKAT_NULL_Model_ADJ") {
      msg <- sprintf("The small sample adjustment only can be applied for linear and linear.weighted kernel in the current version of SKAT! No adjustment is applied")
      warning(msg, call. = FALSE)
      obj <- obj$re1
    }
  }
  if (class(obj) == "SKAT_NULL_Model_EMMAX") {
    re = SKAT:::SKAT_emmaX(Z, obj, kernel = kernel, method = method,
                    weights.beta = weights.beta, weights = weights, impute.method = impute.method,
                    r.corr = r.corr, is_check_genotype = is_check_genotype,
                    is_dosage = is_dosage, missing_cutoff = missing_cutoff,
                    max_maf = max_maf, estimate_MAF = estimate_MAF)
  } else if (class(obj) == "SKAT_NULL_Model_ADJ") {
    re <- SKAT:::SKAT_With_NullModel_ADJ(Z, obj, kernel = kernel,
                                  method = method, weights.beta = weights.beta, weights = weights,
                                  impute.method = impute.method, r.corr = r.corr, is_check_genotype = is_check_genotype,
                                  is_dosage = is_dosage, missing_cutoff = missing_cutoff,
                                  max_maf = max_maf, estimate_MAF = estimate_MAF)
  } else if (class(obj) == "SKAT_NULL_Model") {
    re <- SKAT:::SKAT_With_NullModel(Z, obj, kernel = kernel, method = method,
                              weights.beta = weights.beta, weights = weights, impute.method = impute.method,
                              r.corr = r.corr, is_check_genotype = is_check_genotype,
                              is_dosage = is_dosage, missing_cutoff = missing_cutoff,
                              max_maf = max_maf, estimate_MAF = estimate_MAF)
  } else if (class(obj) == "ziVC_NULL_Model") {
    re <- ziVCtest_With_NullModel(Z, obj, kernel = kernel, method = method,
                                     weights.beta = weights.beta, weights = weights, impute.method = impute.method,
                                     r.corr = r.corr, is_check_genotype = is_check_genotype,
                                     is_dosage = is_dosage, missing_cutoff = missing_cutoff,
                                     max_maf = max_maf, estimate_MAF = estimate_MAF)
  }
  else if (class(obj) == "SKAT_NULL_Model_Hurdle") {
    print("entering vcHurdle_With_NullModel")
    re <- vcHurdle_With_NullModel(Z, obj, kernel = kernel, method = method,
                                  weights.beta = weights.beta, weights = weights, impute.method = impute.method,
                                  r.corr = r.corr, is_check_genotype = is_check_genotype,
                                  is_dosage = is_dosage, missing_cutoff = missing_cutoff,
                                  max_maf = max_maf, estimate_MAF = estimate_MAF)
  } else {
    stop("Please run SKAT_NULL_Model or ziVICtest_Null_Model first!")
  }
  class(re) <- "SKAT_OUT"
  return(re)
}

