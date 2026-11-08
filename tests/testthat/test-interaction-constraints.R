# interaction_constraints: a tree on the features T is grown only if T is a subset of a given set.

set.seed(149)
n_ic <- 300
X_ic <- matrix(runif(n_ic * 5), n_ic, 5, dimnames = list(NULL, paste0("x", 1:5)))
y_ic <- X_ic[, 1] * X_ic[, 2] + X_ic[, 3] * X_ic[, 4] * X_ic[, 5] + X_ic[, 5] + rnorm(n_ic, sd = 0.1)
yc_ic <- factor(y_ic > median(y_ic))
structures_ic <- c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees")

# the split dimensions (1-based features) of every tree with leaves, over all families
tree_dims <- function(fit) {
  out <- list()
  for (fam in fit$fit$get_model()) {
    for (v in fam$variables) {
      v <- sort(v[v > 0])
      if (length(v) > 0) out[[paste(v, collapse = ",")]] <- v
    }
  }
  unname(out)
}

allowed <- function(dims, sets) any(vapply(sets, function(s) all(dims %in% s), logical(1)))

leaves_of <- function(fit) {
  g <- fit$fit$get_model()
  unlist(lapply(g, function(fam) c(unlist(fam$values), unlist(fam$intervals))))
}

test_that("no tree lies outside the allowed sets, for every split structure, regression and classification", {
  sets <- list(c("x1", "x2"), c("x3", "x4"), "x5")
  idx <- list(1:2, 3:4, 5)
  for (ss in structures_ic) {
    for (task in c("regression", "L2", "logit", "exponential")) {
      set.seed(3)
      fit <- if (task == "regression") {
        rpf(
          X_ic,
          y_ic,
          max_interaction = 0,
          ntrees = 5,
          splits = 60,
          split_structure = ss,
          interaction_constraints = sets
        )
      } else {
        rpf(
          X_ic,
          yc_ic,
          max_interaction = 0,
          ntrees = 5,
          splits = 60,
          split_structure = ss,
          loss = task,
          interaction_constraints = sets
        )
      }
      dims <- tree_dims(fit)
      info <- paste(ss, task)
      expect_true(all(vapply(dims, allowed, logical(1), sets = idx)), info = info)
      # the allowed pairs are grown (the data have the interaction {x1, x2})
      expect_true(any(vapply(dims, function(d) identical(d, 1:2), logical(1))), info = info)
    }
  }
})

test_that("interaction_constraints composes with max_interaction", {
  for (ss in structures_ic) {
    set.seed(4)
    fit <- suppressMessages(rpf(
      X_ic,
      y_ic,
      max_interaction = 2,
      ntrees = 5,
      splits = 80,
      split_structure = ss,
      interaction_constraints = list(1:2, 3:5)
    ))
    dims <- tree_dims(fit)
    expect_true(all(lengths(dims) <= 2), info = ss)
    expect_true(all(vapply(dims, allowed, logical(1), sets = list(1:2, 3:5))), info = ss)
  }
  # without max_interaction's cap, the triple {x3, x4, x5} can be grown
  set.seed(4)
  fit <- rpf(X_ic, y_ic, max_interaction = 0, ntrees = 5, splits = 200, interaction_constraints = list(1:2, 3:5))
  expect_true(max(lengths(tree_dims(fit))) == 3)
})

test_that("a feature in no set is not used", {
  set.seed(5)
  fit <- rpf(X_ic, y_ic, max_interaction = 0, ntrees = 5, splits = 50, interaction_constraints = list(c(1, 2)))
  expect_true(all(unlist(tree_dims(fit)) %in% 1:2))
})

test_that("names, indices and the formula interface give the same fit", {
  set.seed(6)
  a <- rpf(
    X_ic,
    y_ic,
    max_interaction = 0,
    ntrees = 3,
    splits = 40,
    interaction_constraints = list(c("x1", "x2"), "x3")
  )
  set.seed(6)
  b <- rpf(
    X_ic,
    y_ic,
    max_interaction = 0,
    ntrees = 3,
    splits = 40,
    interaction_constraints = list(c(2, 1), 3, c(3, 3))
  )
  dat <- data.frame(y = y_ic, X_ic)
  set.seed(6)
  d <- rpf(
    y ~ .,
    data = dat,
    max_interaction = 0,
    ntrees = 3,
    splits = 40,
    interaction_constraints = list(c("x1", "x2"), "x3")
  )
  expect_identical(leaves_of(a), leaves_of(b))
  expect_identical(leaves_of(a), leaves_of(d))
  expect_identical(a$params$interaction_constraints, list(c("x1", "x2"), "x3"))
  expect_output(print(a), "Interaction constraints.*x1, x2")
})

test_that("without interaction_constraints the fit is unchanged; all sets allowed equals no constraint", {
  for (ss in structures_ic) {
    set.seed(7)
    a <- rpf(X_ic, y_ic, max_interaction = 3, ntrees = 3, splits = 40, split_structure = ss)
    set.seed(7)
    b <- rpf(
      X_ic,
      y_ic,
      max_interaction = 3,
      ntrees = 3,
      splits = 40,
      split_structure = ss,
      interaction_constraints = NULL
    )
    set.seed(7)
    d <- suppressMessages(rpf(
      X_ic,
      y_ic,
      max_interaction = 3,
      ntrees = 3,
      splits = 40,
      split_structure = ss,
      interaction_constraints = list(1:5)
    ))
    expect_identical(leaves_of(a), leaves_of(b), info = ss)
    expect_identical(leaves_of(a), leaves_of(d), info = ss)
    expect_null(a$params$interaction_constraints)
  }
  # classification too
  set.seed(8)
  a <- rpf(X_ic, yc_ic, max_interaction = 2, ntrees = 3, splits = 40, loss = "logit")
  set.seed(8)
  d <- suppressMessages(rpf(
    X_ic,
    yc_ic,
    max_interaction = 2,
    ntrees = 3,
    splits = 40,
    loss = "logit",
    interaction_constraints = list(1:5)
  ))
  expect_identical(leaves_of(a), leaves_of(d))
})

test_that("bad interaction_constraints give clear errors", {
  fit <- function(ic) rpf(X_ic, y_ic, ntrees = 2, splits = 5, interaction_constraints = ic)
  expect_error(fit(c("x1", "x2")), "must be .*NULL.* or a non-empty list")
  expect_error(fit(list()), "must be .*NULL.* or a non-empty list")
  expect_error(fit(list(c("x1", "nope"))), "unknown predictor")
  expect_error(fit(list(1:2, 7)), "between 1 and 5")
  expect_error(fit(list(0)), "between 1 and 5")
  expect_error(fit(list(1.5)), "between 1 and 5")
  expect_error(fit(list(character(0))), "non-empty vector without")
  expect_error(fit(list(c(1, NA))), "non-empty vector without")
  expect_error(fit(list(TRUE)), "must be predictor names or column indices")
  expect_message(
    rpf(X_ic, y_ic, ntrees = 2, splits = 5, max_interaction = 1, interaction_constraints = list(1:2)),
    "Trees use at most 1 predictor"
  )
})

test_that("with a formula, indices refer to the formula's predictors in the user's order", {
  set.seed(2)
  n <- 200
  df <- data.frame(x1 = runif(n), g = factor(sample(letters[1:3], n, TRUE)), x2 = runif(n))
  df$y <- (df$g == "b") * 2 + df$x1 * df$x2 + rnorm(n, sd = 0.2)
  fit <- rpf(
    y ~ x1 + g + x2,
    data = df,
    interaction_constraints = list(c(1, 2)),
    max_interaction = 0,
    ntrees = 2,
    splits = 20
  )
  # the model matrix puts the factor last (x1, x2, g): the indices must not refer to it
  expect_setequal(fit$params$interaction_constraints[[1]], c("x1", "g"))
  blueprint <- hardhat::default_formula_blueprint(intercept = FALSE, indicators = "none")
  processed <- hardhat::mold(y ~ x1 + g + x2, df, blueprint = blueprint)
  cols <- colnames(preprocess_predictors_fit(processed)$predictors_matrix)
  expect_identical(cols, c("x1", "x2", "g"))
  expect_false(match("x2", cols) %in% unlist(tree_dims(fit)))
  # the same as giving the names, and `y ~ .` follows the columns of `data`
  set.seed(1)
  a <- rpf(y ~ x1 + g + x2, data = df, interaction_constraints = list(c(3, 2)), max_interaction = 0, ntrees = 2)
  set.seed(1)
  b <- rpf(y ~ x1 + g + x2, data = df, interaction_constraints = list(c("x2", "g")), max_interaction = 0, ntrees = 2)
  expect_identical(predict(a, df), predict(b, df))
  fit_dot <- rpf(y ~ ., data = df, interaction_constraints = list(c(2, 3)), max_interaction = 0, ntrees = 1, splits = 5)
  expect_setequal(fit_dot$params$interaction_constraints[[1]], c("g", "x2"))
  expect_error(
    rpf(y ~ x1 + g + x2, data = df, interaction_constraints = list(c(1, 4))),
    "formula's predictors"
  )
  # a term that is not one model column can only be given by column names
  expect_error(
    rpf(y ~ poly(x1, 2) + x2, data = df, interaction_constraints = list(c(1, 2))),
    "not a single model column"
  )
})

test_that("with a recipe, interaction_constraints takes names only", {
  skip_if_not_installed("recipes")
  dat <- data.frame(y = y_ic, X_ic)
  rec <- recipes::recipe(y ~ ., data = dat)
  set.seed(9)
  fit <- rpf(rec, dat, max_interaction = 0, ntrees = 2, splits = 30, interaction_constraints = list(c("x1", "x2")))
  expect_true(all(unlist(tree_dims(fit)) %in% 1:2))
  expect_error(rpf(rec, dat, interaction_constraints = list(1:2)), "by predictor names")
})
