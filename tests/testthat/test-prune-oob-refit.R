# prune_components() with select_by = "oob": R's coefficients are fitted on the out-of-bag (OOB) components, for every
# candidate's score, for the final model and for a given rule and tau; the default `select` is "1se"

oobr_data <- function(n = 500, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
  y <- x[, 1] + 2 * x[, 2] * x[, 3] + 0.5 * x[, 4] + rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

oobr_fit <- function(d, max_interaction = 3, ntrees = 20, splits = 60) {
  set.seed(1)
  rpf(d$x, d$y, max_interaction = max_interaction, ntrees = ntrees, splits = splits)
}

# The OOB sums of purify_oob(), with the components named, and the number of OOB families per row
oobr_sums <- function(fit, x) {
  o <- fit$fit$purify_oob(0L, 1L)
  nm <- randomPlantedForest:::prune_set_names(lapply(o$sets, as.integer), colnames(x))
  colnames(o$sum0) <- colnames(o$sum1) <- nm
  o$m0 <- o$count[, 1]
  o$m1 <- o$count[, 2]
  o$m <- o$m0 + o$m1
  o
}

# Least squares of y on the OOB components of `kept`, over the rows with OOB families
oob_lstsq <- function(fit, x, y, kept) {
  o <- oobr_sums(fit, x)
  use <- o$m > 0
  F_oob <- (o$sum0 + o$sum1)[use, kept, drop = FALSE] / o$m[use]
  stats::lm.fit(cbind(1, F_oob), y[use])$coefficients
}

test_that("R chosen out of bag: the chosen tau's components, coefficients fitted on the OOB components", {
  d <- oobr_data()
  fit <- oobr_fit(d)
  # the defaults choose R here
  p <- prune_components(fit, d$x, d$y, seed = 3)
  expect_identical(p$rule, "R")
  expect_identical(p$select_by, "oob")
  expect_identical(p$select, "1se")
  # the true components
  expect_setequal(names(p$beta), c("x1", "x4", "x2:x3"))
  # least squares on the OOB components
  co <- oob_lstsq(fit, d$x, d$y, names(p$beta))
  expect_equal(unname(c(p$intercept, p$beta)), unname(co), tolerance = 1e-10)
  expect_equal(p$components$beta[p$components$kept], unname(p$beta[p$components$component[p$components$kept]]))

  # R given at the chosen tau is the same model: the same components, the same out-of-bag coefficients
  g <- prune_components(fit, d$x, d$y, rule = "R", tau = p$tau)
  expect_identical(g$select_by, "oob")
  expect_null(g$oob)
  expect_null(g$fold_id)
  expect_identical(g$components, p$components)
  expect_identical(g$sets, p$sets)
  expect_identical(g$beta, p$beta)
  expect_identical(g$intercept, p$intercept)
  expect_identical(predict(g, d$x), predict(p, d$x))
  expect_output(print(g), "fitted out of bag")
  expect_output(print(g), "Given, not chosen")

  # in sample (select_by = "cv"): least squares on the forest's own components at the training rows; the OOB fit
  # shrinks the noisier pair
  gi <- prune_components(fit, d$x, d$y, rule = "R", tau = p$tau, select_by = "cv")
  expect_null(gi$select_by)
  expect_identical(names(gi$beta), names(p$beta))
  own <- predict_components(fit, d$x)$m
  co_own <- stats::lm.fit(cbind(1, as.matrix(own[, names(gi$beta), with = FALSE])), d$y)$coefficients
  expect_equal(unname(c(gi$intercept, gi$beta)), unname(co_own), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(p$beta), unname(gi$beta))))
  expect_lt(p$beta[["x2:x3"]], gi$beta[["x2:x3"]])

  # the prediction is the forest's components with these coefficients
  new <- d$x[1:7, ]
  M <- predict_components(fit, new)$m
  manual <- p$intercept + as.numeric(as.matrix(M[, names(p$beta), with = FALSE]) %*% p$beta)
  expect_equal(predict(p, new)$.pred, manual, tolerance = 1e-10)
  pc <- predict_components(p, new)
  expect_equal(rowSums(pc$m) + pc$intercept, manual, tolerance = 1e-10)

  expect_output(print(p), "Rule R \\(drop \\+ refit one coefficient per component, fitted out of bag\\)")
  expect_no_match(format(gi), "fitted out of bag")
})

test_that("R's OOB scores and final fit use the rows with OOB families, R fitted outside each fold", {
  # with 10 families a few rows are in every bootstrap sample, and some have OOB families in one half only
  d <- oobr_data(n = 300, seed = 2)
  fit <- oobr_fit(d, max_interaction = 2, ntrees = 10, splits = 40)
  fold_id <- rep_len(1:4, nrow(d$x))
  p <- prune_components(fit, d$x, d$y, rule = "R", tau = c(0, 0.01), fold_id = fold_id, seed = 1, select = "min")
  o <- oobr_sums(fit, d$x)
  use <- o$m > 0
  both <- (o$m0 > 0 & o$m1 > 0)[use]
  expect_true(any(!use))
  expect_true(any(!both))
  expect_equal(attr(p$oob, "n_rows"), sum(use))

  share <- stats::setNames(p$components$share, p$components$component)
  y <- d$y[use]
  fid <- fold_id[use]
  m0 <- o$m0[use]
  m1 <- o$m1[use]
  m <- o$m[use]
  for (tau in c(0, 0.01)) {
    kept <- names(share)[share >= tau]
    S0 <- o$sum0[use, kept, drop = FALSE]
    S1 <- o$sum1[use, kept, drop = FALSE]
    F_oob <- (S0 + S1) / m
    g <- g0 <- g1 <- numeric(sum(use))
    for (k in 1:4) {
      te <- fid == k
      cf <- stats::lm.fit(cbind(1, F_oob[!te, , drop = FALSE]), y[!te])$coefficients
      g[te] <- cbind(1, F_oob[te, , drop = FALSE]) %*% cf
      g0[te] <- cbind(1, S0[te, , drop = FALSE] / pmax(m0[te], 1)) %*% cf
      g1[te] <- cbind(1, S1[te, , drop = FALSE] / pmax(m1[te], 1)) %*% cf
    }
    v <- (g0 - g1)^2 / (1 / pmax(m0, 1) + 1 / pmax(m1, 1))
    v[!both] <- mean(v[both])
    e2 <- (y - g)^2
    row <- p$oob$tau == tau
    expect_equal(p$oob$oob_rmse_raw[row], sqrt(mean(e2)), tolerance = 1e-10)
    expect_equal(p$oob$oob_rmse[row], sqrt(mean(e2 - (1 / m - 1 / 10) * v)), tolerance = 1e-10)
  }
  co <- oob_lstsq(fit, d$x, d$y, names(p$beta))
  expect_equal(unname(c(p$intercept, p$beta)), unname(co), tolerance = 1e-10)
})

test_that("D's scores do not depend on R's fit, and D and CV choices are fitted as without a choice", {
  d <- oobr_data(seed = 5)
  fit <- oobr_fit(d, max_interaction = 2, splits = 50)
  # the D rows of a choice among both rules are those of a choice among D's thresholds alone
  pDR <- prune_components(fit, d$x, d$y, seed = 4)
  pD <- prune_components(fit, d$x, d$y, rule = "D", seed = 4)
  score_cols <- c("tau", "oob_rmse", "se_mse", "oob_rmse_raw", "se_mse_raw", "n_kept")
  expect_identical(pDR$oob[pDR$oob$rule == "D", score_cols], pD$oob[score_cols], ignore_attr = TRUE)
  # D chosen out of bag: the forest's intercept and kept components
  expect_identical(pD$rule, "D")
  gD <- prune_components(fit, d$x, d$y, rule = "D", tau = pD$tau)
  expect_identical(pD$beta, gD$beta)
  expect_identical(pD$intercept, gD$intercept)
  expect_identical(pD$components, gD$components)
  expect_identical(predict(pD, d$x), predict(gD, d$x))
  out <- capture.output(print(pD))
  expect_false(any(grepl("fitted out of bag", out)))
  expect_true(any(grepl("Rule D \\(drop\\)", out)))
  # R chosen by CV: in sample, as R given at the chosen tau with select_by = "cv"
  pC <- prune_components(fit, d$x, d$y, rule = "R", tau = c(0, 0.01), select_by = "cv", seed = 4)
  gC <- prune_components(fit, d$x, d$y, rule = "R", tau = pC$tau, select_by = "cv")
  expect_identical(pC$beta, gC$beta)
  expect_identical(pC$intercept, gC$intercept)
  expect_no_match(format(pC), "fitted out of bag")
})

test_that("a given R falls back to in-sample coefficients where the choice falls back to CV", {
  d <- oobr_data(n = 300, seed = 3)
  in_sample <- function(fit, p) {
    own <- predict_components(fit, d$x)$m
    stats::lm.fit(cbind(1, as.matrix(own[, names(p$beta), with = FALSE])), d$y)$coefficients
  }
  # fewer than 10 families: in sample, with a message; an explicit select_by = "oob" warns and fits out of bag
  f8 <- oobr_fit(d, max_interaction = 2, ntrees = 8, splits = 40)
  expect_message(p8 <- prune_components(f8, d$x, d$y, rule = "R", tau = 0.01), "fitted in sample")
  expect_null(p8$select_by)
  expect_equal(unname(c(p8$intercept, p8$beta)), unname(in_sample(f8, p8)), tolerance = 1e-10)
  expect_warning(o8 <- prune_components(f8, d$x, d$y, rule = "R", tau = 0.01, select_by = "oob"), "fewer than 10")
  expect_identical(o8$select_by, "oob")
  co8 <- oob_lstsq(f8, d$x, d$y, names(o8$beta))
  expect_equal(unname(c(o8$intercept, o8$beta)), unname(co8), tolerance = 1e-10)
  # a given D needs nothing out of bag: no message
  expect_silent(pD <- prune_components(f8, d$x, d$y, rule = "D", tau = 0.01))
  expect_null(pD$select_by)

  # a deterministic forest: in sample, silently; an explicit select_by = "oob" is an error
  set.seed(1)
  det <- rpf(d$x, d$y, max_interaction = 2, ntrees = 1, splits = 20, deterministic = TRUE)
  expect_silent(pdet <- prune_components(det, d$x, d$y, rule = "R", tau = 0.01))
  expect_null(pdet$select_by)
  expect_equal(unname(c(pdet$intercept, pdet$beta)), unname(in_sample(det, pdet)), tolerance = 1e-10)
  expect_error(prune_components(det, d$x, d$y, rule = "R", tau = 0.01, select_by = "oob"), "deterministic")

  # no in-bag counts (restored from a blob of an earlier version, with its data): in sample, with a message
  f20 <- oobr_fit(d, max_interaction = 2, splits = 40)
  blob <- rpf_marshal(f20, include_data = TRUE)
  blob$fit_state$inbag <- NULL
  expect_message(pn <- prune_components(rpf_unmarshal(blob), d$x, d$y, rule = "R", tau = 0.01), "no in-bag counts")
  expect_null(pn$select_by)
  pin <- prune_components(f20, d$x, d$y, rule = "R", tau = 0.01, select_by = "cv")
  expect_equal(pn$beta, pin$beta, tolerance = 1e-12)
  expect_error(
    prune_components(rpf_unmarshal(blob), d$x, d$y, rule = "R", tau = 0.01, select_by = "oob"),
    "no in-bag counts"
  )
})

test_that("the default 1se rule is sparser than min, on the same scores", {
  d <- oobr_data(seed = 2)
  fit <- oobr_fit(d)
  p1 <- prune_components(fit, d$x, d$y, seed = 3)
  pm <- prune_components(fit, d$x, d$y, seed = 3, select = "min")
  expect_identical(p1$oob, pm$oob)
  expect_identical(c(p1$rule, p1$tau), unlist(p1$oob[p1$oob$one_se, c("rule", "tau")], use.names = FALSE))
  expect_identical(c(pm$rule, pm$tau), unlist(pm$oob[pm$oob$min, c("rule", "tau")], use.names = FALSE))
  # 1se keeps the true components only; min also keeps small false ones here
  expect_setequal(names(p1$beta), c("x1", "x4", "x2:x3"))
  expect_gt(length(pm$beta), length(p1$beta))
  expect_true(all(names(p1$beta) %in% names(pm$beta)))
})
