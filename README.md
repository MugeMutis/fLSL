# fLSL: Functional Spatial-Index Logit Regression

R implementation of the methods in *Functional Spatial-Index Logit Regression:
Inversion-Free Linearized and Iterated GMM Estimation*.

In the model, a scalar-on-function index is propagated across spatial units by
a spatial autoregressive multiplier and linked to the success probability of a
binary response through a logistic function,
`p_i = Lambda((S_rho m)_i / s_i(rho))`, `S_rho = (I - rho W)^{-1}`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("MugeMutis/fLSL")
```

## Main functions

| Function | Purpose | Paper |
|---|---|---|
| `FLSL()` | FPCA, one-step estimator and inversion-free Wald test of rho = 0 | Sec. 3.2, Cor. 3.1 |
| `IFLSL_LS()` | Projected GMM by the safeguarded Gauss-Newton iteration | Sec. 3.3, Alg. 1 |
| `vcov_theta()`, `wald_rho()` | Covariance and Wald tests, conditional independence or spatial HAC | Thm. 3.1, Sec. 3.5 |
| `hac_omega()` | Spatial HAC moment covariance (Bartlett kernel, sparse sum) | eq. (3.11) |
| `functional_effects()` | Direct, indirect and total effect curves with pointwise bands | Sec. 2.2, 3.6 |
| `identification_r2()`, `moran_test()`, `residual_correlogram()` | Diagnostics | Sec. 5 |
| `spatial_benchmark()`, `NLGMM()`, `FLR()` | QMLE/NLS, direct nonlinear GMM and nonspatial benchmarks | Sec. 3.1, 4.1 |
| `dgp()` | Main, exact and dependent-response designs of the experiments | Sec. 4 |

## Example

```r
library(fLSL)
set.seed(1)
dat <- dgp(n = 500, rho = 0.25)                  # main design of Section 4.1

fit1 <- FLSL(dat$y, dat$x, dat$w, gp = dat$gp)   # one-step estimator and test
fit1

fit2 <- IFLSL_LS(flsl_fit = fit1)                # iterated estimator
wald_rho(fit2)

# spatial HAC inference (coordinates in the units of the bandwidth)
wald_rho(fit1, type = "hac", coords = dat$coords, bandwidth = 0.2)

# functional effect curves with pointwise bands
eff <- functional_effects(fit2)
eff
plot(eff)

# fixed-location training/test split: withheld responses may be NA
test <- sample(500, 100); y_train <- dat$y; y_train[test] <- NA
fit_tr <- IFLSL_LS(flsl_fit = FLSL(y_train, dat$x, dat$w, gp = dat$gp,
                                    train = setdiff(1:500, test)))
p_test <- fitted_probability(fit_tr)[test]
```

The FLSL estimator and its test use only sparse products with `W`; the
iterated estimator, the benchmarks and the effect multipliers evaluate the
dense spatial multiplier and are intended for samples in which an n x n
inverse is feasible.
