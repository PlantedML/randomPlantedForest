# The out-of-bag (OOB) path along the splits: rpf(oob_path = ), and the prefix property behind it

structures <- c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")

path_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n))
  y <- 2 * x[, 1] + sin(4 * x[, 2]) + 3 * x[, 1] * x[, 2] + 0.3 * x[, 3] + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

path_fit <- function(d, splits, y = d$y, ...) {
  set.seed(5)
  rpf_050(d$x, y, max_interaction = 2, n_families = 6, n_splits = splits, ...)
}

# One family's prediction at the rows of X, directly from its leaves
family_pred <- function(fam, X) {
  out <- numeric(nrow(X))
  for (t in seq_along(fam$variables)) {
    for (l in seq_along(fam$values[[t]])) {
      iv <- fam$intervals[[t]][[l]]
      ok <- rep(TRUE, nrow(X))
      for (u in seq_len(ncol(X))) {
        ok <- ok & X[, u] >= iv[1, u] & X[, u] < iv[2, u]
      }
      out <- out + fam$values[[t]][[l]][1] * ok
    }
  }
  out
}

# The OOB sums of a fit, brute force: each family's prediction from its leaves at its OOB rows, summed over the
# families of even and of odd index
brute_oob_sums <- function(fit, X) {
  g <- fit$fit$get_model()
  ib <- fit$fit$inbag()
  s <- list(numeric(nrow(X)), numeric(nrow(X)))
  for (a in seq_along(g)) {
    h <- (a - 1) %% 2 + 1
    s[[h]] <- s[[h]] + ifelse(ib[, a] == 0, family_pred(g[[a]], X), 0)
  }
  s
}

# Every tree's leaves (values and intervals), keyed by family and the tree's variables
leaf_list <- function(fit) {
  g <- fit$fit$get_model()
  out <- list()
  for (a in seq_along(g)) {
    for (t in seq_along(g[[a]]$variables)) {
      key <- paste(a, paste(g[[a]]$variables[[t]], collapse = ","))
      out[[key]] <- lapply(seq_along(g[[a]]$values[[t]]), function(l) {
        c(g[[a]]$values[[t]][[l]], as.vector(g[[a]]$intervals[[t]][[l]]))
      })
    }
  }
  out
}

test_that("the prefix property: the k-split forest is the first k splits of a longer fit", {
  # the growth, before the level a classification fit adds after its last split
  # With delete_leaves = FALSE a split only appends leaves, so every tree of the k-split forest must be the first
  # leaves of that tree in the longer fit: every structure, regression and classification.
  d <- path_sim()
  yc <- factor(d$y > stats::median(d$y), labels = c("a", "b"))
  for (s in structures) {
    for (y in list(d$y, yc)) {
      big <- leaf_list(path_fit(d, 60, y = y, split_structure = s, delete_leaves = FALSE, level = FALSE))
      for (k in c(1, 7, 30)) {
        small <- leaf_list(path_fit(d, k, y = y, split_structure = s, delete_leaves = FALSE, level = FALSE))
        is_prefix <- vapply(
          names(small),
          function(t) {
            length(big[[t]]) >= length(small[[t]]) && identical(small[[t]], big[[t]][seq_along(small[[t]])])
          },
          logical(1)
        )
        expect_true(all(is_prefix), label = paste(s, if (is.factor(y)) "classification" else "regression", k))
      }
    }
  }
})

test_that("the OOB path at checkpoint k equals the OOB sums of a separate fit with k splits", {
  d <- path_sim()
  cp <- c(1, 4, 13, 40, 70)
  for (s in structures) {
    for (del in c(TRUE, FALSE)) {
      big <- path_fit(d, 70, split_structure = s, delete_leaves = del, oob_path = cp)
      P <- big$fit$oob_path()
      expect_identical(P$checkpoints, as.integer(cp))
      for (c in seq_along(cp)) {
        sums <- brute_oob_sums(path_fit(d, cp[c], split_structure = s, delete_leaves = del), d$x)
        expect_equal(P$sum0[, c], sums[[1]], tolerance = 1e-12)
        expect_equal(P$sum1[, c], sums[[2]], tolerance = 1e-12)
      }
      ib <- big$fit$inbag()
      expect_identical(P$count[, 1], as.integer(rowSums(ib[, c(TRUE, FALSE)] == 0)))
      expect_identical(P$count[, 2], as.integer(rowSums(ib[, c(FALSE, TRUE)] == 0)))
    }
  }
})

test_that("the rows are really out of bag: a row's OOB prediction does not depend on its own response", {
  d <- path_sim()
  P <- path_fit(d, 60, oob_path = c(10, 60))$fit$oob_path()
  for (i in c(1, 17, 150)) {
    y2 <- d$y
    y2[i] <- y2[i] + 100
    P2 <- path_fit(d, 60, y = y2, oob_path = c(10, 60))$fit$oob_path()
    expect_identical(P2$sum0[i, ], P$sum0[i, ])
    expect_identical(P2$sum1[i, ], P$sum1[i, ])
    # other rows do change (the families that drew row i fit differently)
    expect_false(identical(P2$sum0[-i, ], P$sum0[-i, ]))
  }
})

test_that("the OOB path does not depend on the number of threads", {
  d <- path_sim()
  a <- path_fit(d, 80, oob_path = TRUE, nthreads = 1)
  b <- path_fit(d, 80, oob_path = TRUE, nthreads = 2)
  expect_identical(a$fit$oob_path(), b$fit$oob_path())
  expect_identical(a$oob_path, b$oob_path)
})

test_that("a fit with the path is the fit without; the default fit has no path", {
  d <- path_sim()
  for (s in structures) {
    plain <- path_fit(d, 50, split_structure = s)
    with_path <- path_fit(d, 50, split_structure = s, oob_path = TRUE)
    expect_false("oob_path" %in% names(plain))
    expect_identical(predict(plain, d$x), predict(with_path, d$x))
    expect_identical(plain$fit$get_model(), with_path$fit$get_model())
    expect_identical(plain$fit$inbag(), with_path$fit$inbag())
  }
  expect_error(plain$fit$oob_path(), "no OOB path")
})

test_that("the path table: corrected and plain OOB errors from the sums", {
  d <- path_sim()
  fit <- path_fit(d, 50, oob_path = c(5, 50))
  P <- fit$fit$oob_path()
  tab <- fit$oob_path
  expect_s3_class(tab, "data.frame")
  expect_named(tab, c("n_splits", "oob_rmse", "oob_rmse_raw", "se_mse", "se_mse_raw"))
  expect_identical(tab$n_splits, c(5L, 50L))
  m0 <- P$count[, 1]
  m1 <- P$count[, 2]
  m <- m0 + m1
  ok <- m > 0
  both <- m0 > 0 & m1 > 0
  for (c in 1:2) {
    g <- (P$sum0[, c] + P$sum1[, c]) / m
    e2 <- (d$y - g)^2
    v <- ifelse(both, (P$sum0[, c] / m0 - P$sum1[, c] / m1)^2 / (1 / m0 + 1 / m1), 0)
    corr <- e2 - (1 / m - 1 / 6) * v
    expect_equal(tab$oob_rmse_raw[c], sqrt(mean(e2[ok])))
    expect_equal(tab$oob_rmse[c], sqrt(mean(corr[ok])))
    expect_equal(tab$se_mse[c], stats::sd(corr[ok]) / sqrt(sum(ok)))
    expect_equal(tab$se_mse_raw[c], stats::sd(e2[ok]) / sqrt(sum(ok)))
  }
  expect_equal(attr(tab, "n_rows"), sum(ok))
  # the forest's own OOB prediction at the end: the mean over each row's OOB families
  s <- brute_oob_sums(fit, d$x)
  expect_equal(tab$oob_rmse_raw[2], sqrt(mean((d$y - (s[[1]] + s[[2]]) / m)[ok]^2)))
})

test_that("oob_path works with every interface and survives marshaling", {
  d <- path_sim()
  ref <- path_fit(d, 40, oob_path = TRUE)
  df <- data.frame(d$x, y = d$y)
  set.seed(5)
  f1 <- rpf_050(y ~ ., data = df, max_interaction = 2, n_families = 6, n_splits = 40, oob_path = TRUE)
  set.seed(5)
  f2 <- rpf_050(df[, 1:4], df$y, max_interaction = 2, n_families = 6, n_splits = 40, oob_path = TRUE)
  expect_identical(f1$oob_path, ref$oob_path)
  expect_identical(f2$oob_path, ref$oob_path)
  restored <- rpf_unmarshal(rpf_marshal(ref))
  expect_identical(restored$oob_path, ref$oob_path)
})

test_that("oob_path arguments are checked", {
  d <- path_sim()
  expect_error(path_fit(d, 20, oob_path = c(5, 3)), "increasing")
  expect_error(path_fit(d, 20, oob_path = c(5, 30)), "increasing")
  expect_error(path_fit(d, 20, oob_path = 2.5), "increasing")
  expect_error(path_fit(d, 20, oob_path = NA), "increasing")
  expect_error(rpf_050(d$x, d$y, n_families = 1, n_splits = 20, oob_path = TRUE, deterministic = TRUE), "deterministic")
  # classification: the path of its criterion (test-classification-pipeline.R)
  yc <- factor(d$y > stats::median(d$y))
  expect_named(
    rpf_050(d$x, yc, n_splits = 20, oob_path = TRUE)$oob_path,
    c("n_splits", "oob_loss", "oob_loss_raw", "se_loss", "se_loss_raw")
  )
  # the core checks the checkpoints too (the params as stored: the old names)
  fit <- new(
    randomPlantedForest:::RandomPlantedForest,
    randomPlantedForest:::rpf_param_vector(
      list(
        max_interaction = 2,
        ntrees = 2,
        splits = 10,
        split_try = 10,
        t_try = 0.4,
        purify = FALSE,
        deterministic = FALSE,
        nthreads = 1,
        split_decay_rate = 0.1,
        max_candidates = 50,
        delete_leaves = TRUE,
        split_structure = "leaves"
      ),
      "regression"
    )
  )
  expect_error(fit$set_oob_path(c(3L, 11L)), "increasing")
})

test_that("oob_path_checkpoints(): four per doubling, ending at splits", {
  expect_identical(oob_path_checkpoints(10), 1:10)
  cp <- oob_path_checkpoints(3200)
  expect_true(all(c(12L, 25L, 50L, 100L, 200L, 400L, 800L, 1600L, 3200L) %in% cp))
  expect_identical(cp, sort(unique(cp)))
  expect_identical(head(cp, 3), c(6L, 7L, 9L))
  expect_identical(tail(oob_path_checkpoints(1000), 2), c(951L, 1000L))
  expect_error(oob_path_checkpoints(0))
})
