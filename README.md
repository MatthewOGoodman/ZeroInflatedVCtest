# ZeroInflatedVCtest

Variance-component score tests for SNP-set association with zero-inflated outcomes, extending the sequence kernel association test (SKAT).

Code accompanying:

> Goodman M, Chibnik L, Cai T. **Variance components genetic association test for zero-inflated count outcomes.**
> *Genetic Epidemiology* 2019;43(1):82–101. doi:[10.1002/gepi.22162](https://doi.org/10.1002/gepi.22162)

For a set of variants, the package aggregates single-variant score statistics into a SNP-set statistic and returns p-values for the zero part and the non-zero part of the outcome, plus an omnibus test of association. Supported null models are zero-inflated Poisson (`"ZIP"`), zero-inflated negative binomial (`"ZINB"`) and a hurdle model for zero-inflated continuous outcomes (`"Hurdle"`). The null distribution of the residuals is estimated once by resampling and reused for every SNP set.

## Install

```r
# install.packages("remotes")
remotes::install_github("MatthewOGoodman/ZeroInflatedVCtest", build_vignettes = TRUE)
```

Depends on `SKAT`, `pscl`, `CompQuadForm`, `MASS` and the Bioconductor package `EmpiricalBrownsMethod`; the command above installs all of them.

## Usage

```r
library(ZeroInflatedVCtest)

# Fit the null model once (no genetic effect); raise n.Resampling for precision
null_obj <- ziVCtest_Null_Model(formula = Y ~ X, data = dat, out_type = "ZIP",
                                n.Resampling = 1000, type.Resampling = "perturbation")

# Test each SNP set: G is an n x p genotype matrix
res <- VCtest(G, null_obj, kernel = "linear", method = "davies.precise")
res$p.val.pi   # zero-inflation component
res$p.val.ld   # count (mean) component
```

The vignette (`vignette("ziVCvignette")`) runs a full simulated example under the null and the alternative.

## Licence

GPL (>= 2), as declared in `DESCRIPTION`.
