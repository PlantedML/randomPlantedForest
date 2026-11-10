# prune_components(max_panels = ): the display. Every kept component of order >= 3 is coarsened to at most
# max_panels heat maps, then the same rule prunes again.

pd_data <- function(n = 400, p = 6, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  y <- x[, 1] + 2 * x[, 2] * x[, 3] * x[, 4] + 2 * x[, 1] * x[, 2] * x[, 5] * x[, 6] + rnorm(n, sd = 0.1)
  list(x = x, y = y)
}

pd_fit <- function(d, max_interaction = ncol(d$x), nthreads = 1L) {
  set.seed(2)
  rpf(
    d$x,
    d$y,
    max_interaction = max_interaction,
    ntrees = 12,
    splits = 150,
    split_try = 5,
    t_try = 0.5,
    nthreads = nthreads
  )
}

# brute force: f_S averaged over the training values of the coarse features in the row's intervals (product of the
# intervals' empirical marginals), evaluated with predict_components_matrix()
pd_brute <- function(fit, X, S, C, cuts, row) {
  grids <- lapply(seq_along(C), function(i) {
    iv <- findInterval(X[row, C[i]], cuts[[i]])
    X[findInterval(X[, C[i]], cuts[[i]]) == iv, C[i]]
  })
  Z <- X[rep(row, prod(lengths(grids))), , drop = FALSE]
  G <- as.matrix(expand.grid(grids))
  for (i in seq_along(C)) {
    Z[, C[i]] <- G[, i]
  }
  mean(fit$fit$predict_components_matrix(Z, list(S), 0L))
}

test_that("the panel form: C = {} is the component, the coarsened values are the cell means of f_S", {
  d <- pd_data()
  fit <- pd_fit(d)
  purify(fit)
  sets <- prune_stored_sets(fit)
  hi <- sets[lengths(sets) >= 3]
  expect_true(length(hi) > 0)
  X <- d$x
  M <- fit$fit$predict_components_matrix(X, hi, 0L)
  for (k in seq_along(hi)) {
    v <- fit$fit$coarsen_eval(list(S = hi[[k]], C = integer(0), cuts = list()), X)
    expect_lt(max(abs(v - M[, k])), 1e-12 * max(1, sd(M[, k])))
  }
  sel <- fit$fit$coarsen_select(hi[1:min(4, length(hi))], 3L, 1L)
  co <- Filter(function(s) s$form == "coarsened", sel)
  expect_true(length(co) > 0)
  for (s in co) {
    expect_lte(s$cells, 3)
    # stored form = the training values used for the choice (bitwise) = coarsen_eval
    expect_identical(rpf_coarse_predict(s$stored, X), s$g)
    expect_equal(fit$fit$coarsen_eval(list(S = s$S, C = s$C, cuts = s$cuts), X), s$g, tolerance = 1e-14)
    for (row in c(3, 77)) {
      expect_equal(s$g[row], pd_brute(fit, X, s$S, s$C, s$cuts, row), tolerance = 1e-10)
    }
    # the choice: the candidate keeping the most variance (the first of ties)
    expect_equal(s$var_kept, max(s$candidates$var_kept))
  }
})

test_that("the OOB sums: with C = {} those of purify_oob(); threads do not change them", {
  d <- pd_data(n = 300)
  fit <- pd_fit(d, 4)
  o <- fit$fit$purify_oob(0L, 1L)
  nm <- vapply(o$sets, paste, "", collapse = ",")
  hi <- o$sets[lengths(o$sets) >= 3]
  oo <- fit$fit$coarsen_oob(lapply(hi, function(S) list(S = S, C = integer(0), cuts = list())), 2L)
  j <- match(vapply(hi, paste, "", collapse = ","), nm)
  expect_lt(max(abs(oo$sum0 - o$sum0[, j])), 1e-12)
  expect_lt(max(abs(oo$sum1 - o$sum1[, j])), 1e-12)
  S <- hi[[1]]
  sp <- list(list(S = S, C = S[1], cuts = list(stats::quantile(d$x[, S[1]], c(0.3, 0.6), names = FALSE))))
  expect_identical(fit$fit$coarsen_oob(sp, 1L), fit$fit$coarsen_oob(sp, 3L))
})

test_that("the C++ checks its input", {
  d <- pd_data(n = 200)
  fit <- pd_fit(d, 3)
  expect_error(fit$fit$coarsen_select(list(1:3), 3L, 1L), "not purified")
  purify(fit)
  expect_error(fit$fit$coarsen_select(list(1:2), 3L, 1L), "at least 3 features")
  expect_error(fit$fit$coarsen_select(list(1:3), 0L, 1L), "at least 1")
  expect_error(fit$fit$coarsen_eval(list(S = 1:3, C = 7L, cuts = list(0.5)), d$x), "subset of S")
  expect_error(fit$fit$coarsen_eval(list(S = 1:3, C = 1L, cuts = list()), d$x), "one vector per feature")
  st <- list(
    S = 1:3,
    cfeat = 1L,
    ccuts = list(0.5),
    dfeat = 2:3,
    dcuts = list(0.5, 0.5),
    box = matrix(0L, 1, 3),
    w = matrix(1, 1, 2),
    n_trees = 1L
  )
  expect_error(rpf_coarse_build(st), "Corrupt")
})

test_that("max_panels: at most P heat maps, decomposition, marginal means, Inf = the first pruning", {
  d <- pd_data()
  fit <- pd_fit(d, 4)
  pd <- prune_components(fit, d$x, d$y, max_panels = 3, seed = 5, nthreads = 2L)
  expect_s3_class(pd, "rpf_pruned")
  disp <- pd$display
  expect_true(disp$second)
  expect_true(length(disp$coarse) > 0)
  for (z in disp$coarse) {
    expect_lte(z$cells, 3)
    expect_lte(prod(lengths(z$cuts) + 1), 3)
    expect_length(z$C, length(z$S) - 2)
  }
  expect_true(all(c("form", "heat_maps") %in% names(pd$components)))
  expect_true(all(pd$components$heat_maps[pd$components$kept] <= 3))
  # the first pruning is the model without the display, and that is Inf's result
  ref <- prune_components(fit, d$x, d$y, max_panels = Inf, seed = 5, nthreads = 2L)
  expect_null(ref$display)
  expect_identical(disp$first$beta, ref$beta)
  expect_identical(disp$first$oob, ref$oob)
  # the second pruning: the same candidates and folds
  expect_identical(pd$fold_id, ref$fold_id)
  expect_identical(pd$oob[, c("rule", "tau")], ref$oob[, c("rule", "tau")])
  expect_true(all(names(pd$beta) %in% names(ref$beta)))
  # predict = intercept + the sum of the components; predict_components agrees
  p <- predict(pd, d$x[1:50, ])$.pred
  pc <- predict_components(pd, d$x[1:50, ])
  expect_equal(pd$intercept + rowSums(as.matrix(pc$m)), p, tolerance = 1e-12)
  expect_identical(colnames(pc$m), names(pd$beta))
  # the coarsened components are what predict uses
  nm <- intersect(names(disp$coarse), names(pd$beta))
  expect_true(length(nm) > 0)
  expect_equal(pc$m[[nm[1]]], pd$beta[[nm[1]]] * rpf_coarse_predict(disp$coarse[[nm[1]]]$stored, d$x[1:50, ]))
  # the marginal means: averaging along j commutes with the cell means (j coarse: g against f_S coarsened on C - {j})
  X <- d$x
  z <- disp$coarse[[1]]
  S <- match(z$S, colnames(X))
  C <- match(z$C, colnames(X))
  gv <- rpf_coarse_predict(z$stored, X)
  gsd <- sqrt(mean((gv - mean(gv))^2))
  for (r in c(10, 200)) {
    for (j in S) {
      Z <- X[rep(r, nrow(X)), , drop = FALSE]
      Z[, j] <- X[, j]
      keep <- C != j
      a <- mean(rpf_coarse_predict(z$stored, Z))
      b <- mean(fit$fit$coarsen_eval(list(S = S, C = C[keep], cuts = unname(z$cuts)[keep]), Z))
      expect_lt(abs(a - b), 1e-12 * gsd)
    }
  }
  # nothing coarsened: the first pruning's model, bitwise
  pmax <- prune_components(fit, d$x, d$y, max_panels = .Machine$integer.max, seed = 5, nthreads = 1L)
  expect_false(pmax$display$second)
  expect_length(pmax$display$coarse, 0)
  expect_identical(pmax$beta, ref$beta)
  expect_identical(pmax$intercept, ref$intercept)
  expect_identical(predict(pmax, d$x)$.pred, predict(ref, d$x)$.pred)
  # threads do not change the result
  pd1 <- prune_components(fit, d$x, d$y, max_panels = 3, seed = 5, nthreads = 1L)
  expect_identical(pd1$beta, pd$beta)
  no_time <- function(co) lapply(co, function(z) z[names(z) != "seconds"])
  expect_identical(no_time(pd1$display$coarse), no_time(pd$display$coarse))
  expect_identical(predict(pd1, d$x)$.pred, predict(pd, d$x)$.pred)
  out <- capture.output(print(pd))
  expect_true(any(grepl("at most 3 heat maps", out)))
  expect_true(any(grepl("heat_maps", out)))
})

test_that("full order: components of order >= 6 coarsened with q vectors containing 1s", {
  skip_on_cran()
  d <- pd_data(n = 500, p = 7, seed = 3)
  fit <- pd_fit(d, nthreads = 2L)
  purify(fit)
  sets <- prune_stored_sets(fit)
  big <- sets[lengths(sets) >= 6]
  expect_true(length(big) > 0)
  sel <- fit$fit$coarsen_select(big, 10L, 2L)
  co <- Filter(function(s) s$form == "coarsened", sel)
  expect_true(length(co) > 0)
  X <- d$x
  # the maximal q vectors of n axes with product <= 10, in lexicographic order
  qv <- function(n) {
    g <- as.matrix(expand.grid(rep(list(1:10), n)))[, n:1, drop = FALSE]
    g <- g[do.call(order, as.data.frame(g)), , drop = FALSE]
    pr <- apply(g, 1, prod)
    ok <- pr <= 10 & apply(g, 1, function(q) all(prod(q) / q * (q + 1) > 10))
    apply(g[ok, , drop = FALSE], 1, paste, collapse = ",")
  }
  for (s in co) {
    k <- length(s$S)
    expect_lte(s$cells, 10)
    expect_length(s$C, k - 2)
    expect_equal(nrow(s$candidates), choose(k, 2) * length(qv(k - 2)))
    expect_identical(s$candidates$q[s$candidates$combo == 1], qv(k - 2))
    # with >= 4 coarse features some feature is averaged out
    expect_true(any(s$q == 1))
    expect_identical(rpf_coarse_predict(s$stored, X), s$g)
    # brute force where the product grid is small (a coarse feature with one interval spans all rows)
    sizes <- vapply(
      seq_along(s$C),
      function(i) sum(findInterval(X[, s$C[i]], s$cuts[[i]]) == findInterval(X[9, s$C[i]], s$cuts[[i]])),
      numeric(1)
    )
    if (prod(sizes) <= 2e5) {
      expect_equal(s$g[9], pd_brute(fit, X, s$S, s$C, s$cuts, 9), tolerance = 1e-10)
    }
  }
  # and through prune_components() on the full-order forest
  pd <- prune_components(fit, d$x, d$y, seed = 1, nthreads = 2L)
  for (z in pd$display$coarse) {
    expect_lte(z$cells, 10)
  }
  expect_equal(
    pd$intercept + rowSums(as.matrix(predict_components(pd, d$x)$m)),
    predict(pd, d$x)$.pred,
    tolerance = 1e-12
  )
})

test_that("the display after a given rule, a CV choice and a deterministic forest", {
  d <- pd_data(n = 300)
  fit <- pd_fit(d, 4)
  # a given R: out of bag again, on the coarsened OOB components
  gr <- prune_components(fit, d$x, d$y, max_panels = 3, rule = "R", tau = 0.001)
  expect_true(gr$display$second)
  expect_identical(gr$select_by, "oob")
  expect_identical(gr$rule, "R")
  expect_identical(gr$tau, 0.001)
  expect_null(gr$oob)
  # a given D: the shares of the coarsened components, the forest's intercept
  gd <- prune_components(fit, d$x, d$y, max_panels = 3, rule = "D", tau = 0.001)
  expect_true(gd$display$second)
  expect_true(all(gd$beta == 1))
  expect_identical(gd$intercept, gd$display$first$intercept)
  expect_true(all(gd$components$share[gd$components$kept] >= 0.001))
  # a CV choice: the chosen rule and tau applied again (in sample)
  cv <- prune_components(fit, d$x, d$y, max_panels = 3, select_by = "cv", folds = 3L, seed = 2)
  expect_identical(cv$rule, cv$display$first$rule)
  expect_identical(cv$tau, cv$display$first$tau)
  expect_identical(cv$cv, cv$display$first$cv)
  for (m in list(gr, gd, cv)) {
    expect_equal(
      m$intercept + rowSums(as.matrix(predict_components(m, d$x)$m)),
      predict(m, d$x)$.pred,
      tolerance = 1e-12
    )
  }
  # a deterministic forest: the default scores by CV, the display runs
  set.seed(3)
  fd <- rpf(d$x, d$y, max_interaction = 3, ntrees = 1, splits = 60, deterministic = TRUE)
  pdd <- prune_components(fd, d$x, d$y, max_panels = 2, folds = 3L, seed = 1)
  expect_identical(pdd$select_by, "cv")
  expect_false(is.null(pdd$display))
  for (z in pdd$display$coarse) {
    expect_lte(z$cells, 2)
  }
})

test_that("the panels of the components kept as they are: own cut points, forest_panels, coarse and asis cover", {
  d <- pd_data()
  fit <- pd_fit(d, 4)
  P <- 50
  pd <- prune_components(fit, d$x, d$y, max_panels = P, seed = 5, nthreads = 1L)
  disp <- pd$display
  fo <- disp$forms
  expect_identical(names(disp$asis), fo$component[fo$form == "asis"])
  expect_identical(names(disp$coarse), fo$component[fo$form == "coarsened"])
  expect_true(any(names(disp$asis) %in% names(pd$beta)))
  expect_true(any(names(disp$coarse) %in% names(pd$beta)))
  # coarse and asis together: every kept component of three or more predictors, once
  kept_of <- function(z) names(z)[vapply(z, function(e) e$kept, logical(1))]
  shown <- c(kept_of(disp$coarse), kept_of(disp$asis))
  expect_setequal(shown, names(pd$beta)[lengths(pd$sets) >= 3])
  expect_false(anyDuplicated(shown) > 0)
  # the cuts: f_S's own cut points (the leaf bounds of the trees containing S, inside the training range) on the
  # |S| - 2 predictors with the fewest (the first of ties); their cells are the forest_panels
  X <- d$x
  model <- fit$fit$get_model()
  if (!is_purified(fit)) {
    purify(fit)
  }
  for (nm in names(disp$asis)) {
    z <- disp$asis[[nm]]
    S <- match(z$S, colnames(X))
    own <- lapply(S, function(k) {
      b <- unlist(lapply(model, function(fam) {
        lapply(seq_along(fam$variables), function(t) {
          if (all(S %in% fam$variables[[t]])) vapply(fam$intervals[[t]], function(iv) iv[, k], numeric(2))
        })
      }))
      b <- sort(unique(b))
      b[b > min(X[, k]) & b <= max(X[, k])]
    })
    pick <- sort(order(lengths(own))[seq_len(length(S) - 2)])
    expect_identical(z$S, colnames(X)[sort(S)])
    expect_identical(z$C, colnames(X)[S[pick]])
    expect_named(z$cuts, z$C)
    expect_identical(unname(z$cuts), own[pick])
    expect_identical(z$forest_panels, prod(lengths(z$cuts) + 1))
    expect_lte(z$forest_panels, P)
    expect_identical(fo$heat_maps[fo$component == nm], z$forest_panels)
    expect_identical(z$kept, nm %in% names(pd$beta))
    # f_S is constant on each panel: averaged within the cells of these cuts, it is itself
    f <- fit$fit$coarsen_eval(list(S = S, C = integer(0), cuts = list()), X)
    g <- fit$fit$coarsen_eval(list(S = S, C = S[pick], cuts = unname(z$cuts)), X)
    expect_lt(max(abs(g - f)), 1e-12 * max(1, sd(f)))
  }
  # at max_panels = 3 every one is coarsened: asis is empty
  expect_length(prune_components(fit, d$x, d$y, max_panels = 3, seed = 5, nthreads = 1L)$display$asis, 0)
})

test_that("max_panels is checked", {
  d <- pd_data(n = 200)
  fit <- pd_fit(d, 3)
  expect_error(prune_components(fit, d$x, d$y, max_panels = 0), "max_panels")
  expect_error(prune_components(fit, d$x, d$y, max_panels = 2.5), "max_panels")
  expect_error(prune_components(fit, d$x, d$y, max_panels = NA), "max_panels")
})

test_that("predict(): the built indices are cached, rebuilt after serialisation; marshal keeps the display", {
  d <- pd_data()
  fit <- pd_fit(d, 4)
  pd <- prune_components(fit, d$x, d$y, max_panels = 3, seed = 5, nthreads = 1L)
  nm <- intersect(names(pd$display$coarse), names(pd$beta))[1]
  expect_false(is.na(nm))
  p1 <- predict(pd, d$x)$.pred
  expect_true(nm %in% ls(pd$display$cache))
  expect_identical(predict(pd, d$x)$.pred, p1)
  expect_identical(
    rpf_coarse_predict(pd$display$coarse[[nm]]$stored, d$x),
    rpf_coarse_eval(get(nm, envir = pd$display$cache), d$x)
  )
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(pd$display$cache[[nm]], f)
  xp <- readRDS(f)
  expect_false(rpf_coarse_valid(xp))
  expect_error(rpf_coarse_eval(xp, d$x), "gone")
  assign(nm, xp, envir = pd$display$cache)
  expect_identical(predict(pd, d$x)$.pred, p1)
  expect_true(rpf_coarse_valid(get(nm, envir = pd$display$cache)))
  # marshal: the display survives saveRDS() and predicts the same
  blob <- rpf_marshal(pd)
  expect_null(blob$display$cache)
  saveRDS(blob, f)
  back <- rpf_unmarshal(readRDS(f))
  expect_s3_class(back, "rpf_pruned")
  expect_identical(back$display$coarse, pd$display$coarse)
  expect_identical(predict(back, d$x)$.pred, p1)
  expect_identical(as.matrix(predict_components(back, d$x)$m), as.matrix(predict_components(pd, d$x)$m))
})
