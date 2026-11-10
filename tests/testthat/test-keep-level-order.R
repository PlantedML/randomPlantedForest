# keep_level_order: factor predictors coded in the order of levels() instead of by the outcome

klo_data <- function(n = 120, outcome = c("regression", "binary", "multiclass")) {
  outcome <- match.arg(outcome)
  set.seed(3)
  d <- data.frame(
    g = factor(sample(letters[1:5], n, replace = TRUE)),
    o = factor(sample(c("lo", "mid", "hi"), n, replace = TRUE), levels = c("lo", "mid", "hi"), ordered = TRUE),
    x = runif(n)
  )
  # effects that are not monotone in the level order, so that ordering by the outcome changes it
  eta <- c(a = 0, b = 3, c = 0.5, d = 2.5, e = 1)[as.character(d$g)] +
    c(lo = 0, mid = 2, hi = 1)[as.character(d$o)] +
    d$x
  d$y <- switch(
    outcome,
    regression = eta + rnorm(n, sd = 0.3),
    binary = factor(ifelse(eta + rnorm(n, sd = 0.5) > 2.5, "yes", "no")),
    multiclass = factor(cut(eta + rnorm(n, sd = 0.5), c(-Inf, 2, 3.5, Inf), labels = c("A", "B", "C")))
  )
  d
}

# the levels each box of family 1 covers on `variable`, as positions in levels(): a box covers the codes c with
# lower <= c < upper, and code k stands for x$factor_levels[[variable]][k]
klo_box_positions <- function(fit, variable, lev) {
  boxes <- rpf_boxes(fit)
  lower <- boxes[[paste0(variable, "_lower")]]
  upper <- boxes[[paste0(variable, "_upper")]]
  used <- !is.na(lower)
  codes <- seq_along(fit$factor_levels[[variable]])
  lapply(which(used), function(i) {
    covered <- fit$factor_levels[[variable]][codes[codes >= lower[i] & codes < upper[i]]]
    sort(match(covered, lev))
  })
}

is_run <- function(pos) length(pos) == 0 || all(diff(pos) == 1)

test_that("keep_level_order = FALSE is the default and orders every factor by the outcome", {
  d <- klo_data()
  set.seed(1)
  f_default <- rpf_050(y ~ ., data = d, n_families = 2, n_splits = 20)
  set.seed(1)
  f_false <- rpf_050(y ~ ., data = d, n_families = 2, n_splits = 20, keep_level_order = FALSE)
  expect_identical(f_default$factor_levels, f_false$factor_levels)
  expect_identical(predict(f_default, d), predict(f_false, d))
  expect_identical(f_default$params$keep_level_order, character(0))
  # ordered by the mean outcome, ordered factors too
  expect_identical(f_default$factor_levels$g, names(sort(tapply(d$y, d$g, mean))))
  expect_identical(f_default$factor_levels$o, names(sort(tapply(d$y, d$o, mean))))
  expect_false(identical(f_default$factor_levels$g, levels(d$g)))
  expect_false(identical(f_default$factor_levels$o, levels(d$o)))
})

test_that("keep_level_order keeps the order of levels() for all or the named factors", {
  d <- klo_data()
  f_all <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 10, keep_level_order = TRUE)
  expect_identical(f_all$factor_levels, list(g = levels(d$g), o = levels(d$o)))
  expect_identical(f_all$params$keep_level_order, c("g", "o"))

  f_o <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 10, keep_level_order = "o")
  expect_identical(f_o$factor_levels$o, levels(d$o))
  expect_identical(f_o$factor_levels$g, names(sort(tapply(d$y, d$g, mean))))
  expect_identical(f_o$params$keep_level_order, "o")

  # names are stored in the order of the predictors, duplicates once
  f_both <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 10, keep_level_order = c("o", "g", "o"))
  expect_identical(f_both$params$keep_level_order, c("g", "o"))

  # a character predictor becomes a factor with alphabetical levels
  d_chr <- d
  d_chr$g <- as.character(d_chr$g)
  f_chr <- rpf_050(x = d_chr[c("g", "x")], y = d_chr$y, n_families = 1, n_splits = 10, keep_level_order = "g")
  expect_identical(f_chr$factor_levels$g, sort(unique(d_chr$g)))
})

test_that("with keep_level_order, splits only separate runs of neighbouring levels", {
  structures <- c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")
  for (outcome in c("regression", "binary", "multiclass")) {
    d <- klo_data(outcome = outcome)
    for (s in structures) {
      set.seed(2)
      f <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 25, split_structure = s, keep_level_order = TRUE)
      for (v in c("g", "o")) {
        pos <- klo_box_positions(f, v, levels(d[[v]]))
        expect_true(length(pos) > 0, label = paste(outcome, s, v, "is split on"))
        expect_true(all(vapply(pos, is_run, logical(1))), label = paste(outcome, s, v, "boxes are runs"))
      }
    }
  }

  # the check can fail: ordered by the outcome, some box of g covers levels that are not neighbours in levels()
  d <- klo_data()
  set.seed(2)
  f <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 25)
  pos <- klo_box_positions(f, "g", levels(d$g))
  expect_false(all(vapply(pos, is_run, logical(1))))
})

test_that("predictions use the stored coding, whatever the levels of new_data", {
  for (outcome in c("regression", "binary", "multiclass")) {
    d <- klo_data(outcome = outcome)
    for (keep in list(TRUE, "o")) {
      set.seed(4)
      f <- rpf_050(y ~ ., data = d, n_families = 2, n_splits = 20, max_interaction = 2, keep_level_order = keep)
      type <- if (outcome == "regression") "numeric" else "prob"
      ref <- predict(f, d, type = type)

      # levels in another order (an ordered factor must have the fit's levels: hardhat's forge() checks its type)
      d_rev <- d
      d_rev$g <- factor(as.character(d$g), levels = rev(levels(d$g)))
      expect_equal(predict(f, d_rev, type = type), ref)

      # a subset of the levels
      rows <- which(d$g %in% c("b", "d"))
      d_sub <- d[rows, ]
      d_sub$g <- droplevels(d_sub$g)
      expect_equal(predict(f, d_sub, type = type), ref[rows, , drop = FALSE], ignore_attr = TRUE)

      # characters instead of factors
      d_chr <- d
      d_chr$g <- as.character(d$g)
      expect_equal(predict(f, d_chr, type = type), ref)

      # predict_components() codes the same way: intercept + components = the prediction
      comp <- predict_components(f, d_rev)
      comp_ref <- predict_components(f, d)
      expect_equal(comp$m, comp_ref$m)
      if (outcome == "regression") {
        expect_equal(comp$intercept + rowSums(comp$m), predict(f, d)$.pred)
      }

      # purified predictions are the same
      purify(f)
      expect_equal(predict(f, d_rev, type = type), ref)

      # marshal round trip keeps the coding
      restored <- rpf_unmarshal(rpf_marshal(f))
      expect_identical(restored$factor_levels, f$factor_levels)
      expect_identical(restored$params$keep_level_order, f$params$keep_level_order)
      expect_equal(predict(restored, d_rev, type = type), ref)
    }
  }
})

test_that("keep_level_order is validated", {
  d <- klo_data()
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = "x"), "Not a factor")
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = c("g", "zz")), "Not a predictor")
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = NA), "must be")
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = c(TRUE, FALSE)), "must be")
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = 1), "must be")
  expect_error(rpf_050(y ~ ., data = d, n_families = 1, keep_level_order = NA_character_), "must be")
  # TRUE or character(0) without any factor predictor: nothing to keep
  f <- rpf_050(y ~ x, data = d, n_families = 1, n_splits = 5, keep_level_order = TRUE)
  expect_identical(f$params$keep_level_order, character(0))
  f <- rpf_050(y ~ ., data = d, n_families = 1, n_splits = 5, keep_level_order = character(0))
  expect_identical(f$params$keep_level_order, character(0))
})

test_that("print() names the factors that keep their level order", {
  d <- klo_data()
  out_all <- capture.output(print(rpf_050(y ~ ., data = d, n_families = 1, n_splits = 5, keep_level_order = TRUE)))
  expect_true(any(grepl("Factors `g` and `o` keep their level order.", out_all, fixed = TRUE)))
  out_o <- capture.output(print(rpf_050(y ~ ., data = d, n_families = 1, n_splits = 5, keep_level_order = "o")))
  expect_true(any(grepl("Factor `o` keeps its level order; `g` is ordered by the outcome.", out_o, fixed = TRUE)))
  out_none <- capture.output(print(rpf_050(y ~ ., data = d, n_families = 1, n_splits = 5)))
  expect_false(any(grepl("level order", out_none)))
})
