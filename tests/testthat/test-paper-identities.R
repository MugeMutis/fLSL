# Checks of the implementation against identities stated in the paper.

set.seed(1)
dat <- dgp(n = 200, rho = 0.3)
fit <- FLSL(dat$y, dat$x, dat$w, gp = dat$gp)

test_that("FLSL is the one-step estimator (3.6)", {
  tr <- fit$train
  Z <- fit$Z[tr, ]; G <- fit$G[tr, ]
  PZ <- Z %*% solve(crossprod(Z), t(Z))
  step <- solve(t(G) %*% PZ %*% G, t(G) %*% PZ %*% fit$u0)
  expect_equal(unname(fit$theta_hat), unname(fit$theta0 + as.numeric(step)), tolerance = 1e-8)
})

test_that("the Jacobian at rho = 0 is R0 [B, W B delta] (3.5)", {
  o <- fLSL:::.flsl_general_gradient(fit$alpha0, fit$b0, 0, fit$Xi, as.matrix(fit$W))
  expect_equal(unname(o$G), unname(fit$G), tolerance = 1e-10)
})

test_that("the analytic Jacobian (Appendix B) matches numerical derivatives", {
  th <- c(0.1, fit$b0, 0.4); K <- fit$ncomp; Wd <- as.matrix(fit$W)
  o <- fLSL:::.flsl_general_gradient(th[1], th[1 + seq_len(K)], th[K + 2], fit$Xi, Wd)
  num <- sapply(seq_along(th), function(j) {
    e <- numeric(length(th)); e[j] <- 1e-6
    p1 <- fLSL:::.flsl_probability((th + e)[1], (th + e)[1 + seq_len(K)], (th + e)[K + 2], fit$Xi, Wd)$P
    p0 <- fLSL:::.flsl_probability((th - e)[1], (th - e)[1 + seq_len(K)], (th - e)[K + 2], fit$Xi, Wd)$P
    (p1 - p0) / 2e-6
  })
  expect_equal(unname(o$G), unname(num), tolerance = 1e-6)
})

test_that("vcov_theta reproduces the inversion-free V^(0) of Corollary 3.1", {
  vc <- vcov_theta(fit)
  expect_equal(unname(vc$V), unname(fit$wald$V), tolerance = 1e-10)
  w <- wald_rho(fit)
  expect_equal(w$statistic, fit$wald$statistic, tolerance = 1e-10)
})

test_that("HAC with a bandwidth below the smallest distance uses squared residuals only", {
  U <- fit$Z * fit$u0
  Om <- hac_omega(U, dat$coords, bandwidth = 1e-8)
  expect_equal(as.numeric(Om), as.numeric(crossprod(U) / nrow(U)), tolerance = 1e-10)
  expect_equal(attr(Om, "mean_neighbours"), 0)
})

test_that("IFLSL-LS converges and its total multiplier satisfies (2.12)", {
  it <- IFLSL_LS(flsl_fit = fit)
  expect_true(it$converged)
  th <- it$theta_hat; K <- it$ncomp; Wd <- as.matrix(fit$W)
  cm <- effect_multipliers(th, fit$Xi, Wd)
  pr <- fLSL:::.flsl_probability(th[1], th[1 + seq_len(K)], th[K + 2], fit$Xi, Wd)
  expect_equal(unname(cm[["total"]]),
               mean(pr$P * (1 - pr$P) / pr$sigma) / (1 - th[[K + 2]]), tolerance = 1e-8)
  expect_equal(predict(it), fitted_probability(it), tolerance = 1e-10)
  ef <- functional_effects(it)
  expect_equal(nrow(ef$curves), 4 * length(dat$gp))
})

test_that("training subsets: all units reproduce the full fit; NA responses are allowed", {
  f_all <- FLSL(dat$y, dat$x, dat$w, gp = dat$gp, train = rep(TRUE, 200))
  expect_equal(f_all$theta_hat, fit$theta_hat, tolerance = 1e-10)
  y_na <- dat$y; y_na[1:40] <- NA
  f_tr <- FLSL(y_na, dat$x, dat$w, gp = dat$gp, train = 41:200)
  expect_equal(f_tr$n_train, 160L)
  it <- IFLSL_LS(flsl_fit = f_tr)
  expect_length(fitted_probability(it), 200)
})

test_that("the exact design selects K = 4 and the main design uses beta = 2 sin(2 pi t)", {
  set.seed(2)
  ex <- dgp(n = 300, rho = 0, design = "exact")
  expect_equal(FLSL(ex$y, ex$x, ex$w, gp = ex$gp)$ncomp, 4L)
  expect_equal(dat$beta, 2 * sin(2 * pi * dat$gp))
})
