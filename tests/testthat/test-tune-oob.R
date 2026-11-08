# rpf(tune = ): the greediness and min_leaf_rows tuned out of bag in one call. Small data and a short ladder keep
# these fast.

to_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1))
  list(x = x, y = x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3))
}
to_ladder <- data.frame(split_try = c(10L, 5L, 3L, 1L), t_try = c(0.4, 0.2, 0.1, 0.025))
to_tune <- function(d, search = 2:3, min_leaf_rows = c(1L, 20L), neighbours = TRUE, ladder = to_ladder, seed = 3, ...) {
  set.seed(seed)
  rpf(
    d$x,
    d$y,
    max_interaction = 2,
    ntrees = 10,
    tune = tune_control(ladder = ladder, search = search, min_leaf_rows = min_leaf_rows, neighbours = neighbours),
    ...
  )
}
# the model without its search table, the C++ forest (compared through its model), the display's cache of built
# indices and the timings
strip_tune <- function(p) {
  p$tune <- NULL
  p$fit$fit <- NULL
  p$display$cache <- NULL
  p$display$seconds <- NULL
  for (nm in names(p$display$coarse)) {
    p$display$coarse[[nm]]$seconds <- NULL
  }
  p
}

test_that("rpf(tune = ) returns the pruned model of the chosen settings", {
  d <- to_sim()
  r <- to_tune(d)
  expect_s3_class(r, "rpf_pruned")
  p <- predict(r, d$x)
  expect_equal(nrow(p), nrow(d$x))
  expect_true(all(is.finite(p$.pred)))
  tab <- r$tune$search
  fin <- tab[tab$final, ]
  expect_equal(nrow(fin), 1L)
  s <- r$tune$settings
  expect_identical(s, as.list(fin[c("split_try", "t_try", "min_leaf_rows", "splits")]))
  # the final forest is the chosen setting's
  pf <- r$fit$params
  expect_identical(pf$split_try, s$split_try)
  expect_identical(pf$t_try, s$t_try)
  expect_identical(pf$min_leaf_rows, s$min_leaf_rows)
  expect_identical(pf$splits, s$splits)
  expect_identical(r$fit$splits_oob$splits, s$splits)
  expect_true(is.finite(fin$oob_rmse) && fin$oob_rmse > 0)
  expect_s3_class(r$tune$control, "rpf_tune_control")
  # the display ran on the chosen model (max_panels = 10 by default)
  expect_false(is.null(r$display))
})

test_that("a one-point search equals set.seed(s); rpf(splits = NULL, prune = TRUE) at that point", {
  d <- to_sim()
  one <- data.frame(split_try = 5L, t_try = 0.2)
  for (m in c(1L, 20L)) {
    r <- to_tune(d, ladder = one, search = 1, min_leaf_rows = m, seed = 11)
    after_tune <- .Random.seed
    expect_equal(nrow(r$tune$search), 1L)
    set.seed(11)
    p <- rpf(
      d$x,
      d$y,
      max_interaction = 2,
      ntrees = 10,
      splits = NULL,
      split_try = 5L,
      t_try = 0.2,
      min_leaf_rows = m,
      prune = TRUE
    )
    # and leaves the random number stream where that call leaves it
    expect_identical(after_tune, .Random.seed)
    expect_identical(strip_tune(r), strip_tune(p))
    expect_identical(r$fit$fit$get_model(), p$fit$fit$get_model())
    expect_identical(r$fit$params, p$fit$params)
    expect_identical(r$fit$oob_path, p$fit$oob_path)
    expect_identical(predict(r, d$x), predict(p, d$x))
    # the score is the pruned model's corrected OOB RMSE without the display, as prune_components() chose it
    set.seed(11)
    p0 <- rpf(
      d$x,
      d$y,
      max_interaction = 2,
      ntrees = 10,
      splits = NULL,
      split_try = 5L,
      t_try = 0.2,
      min_leaf_rows = m,
      prune = TRUE,
      max_panels = Inf
    )
    i <- which(p0$oob$rule == p0$rule & p0$oob$tau == p0$tau)
    expect_identical(r$tune$search$oob_rmse, p0$oob$oob_rmse[i])
    expect_identical(
      r$tune$search[c("rule", "tau", "n_kept")],
      data.frame(rule = p0$rule, tau = p0$tau, n_kept = length(p0$beta))
    )
    path <- p0$fit$oob_path
    expect_identical(r$tune$search$oob_rmse_forest, path$oob_rmse[match(p0$fit$params$splits, path$splits)])
  }
})

test_that("a one-point search passes split_structure to the forest", {
  d <- to_sim()
  r <- to_tune(
    d,
    ladder = data.frame(split_try = 5L, t_try = 0.2),
    search = 1,
    min_leaf_rows = 1,
    seed = 11,
    split_structure = "cur_trees_1"
  )
  set.seed(11)
  p <- rpf(
    d$x,
    d$y,
    max_interaction = 2,
    ntrees = 10,
    splits = NULL,
    split_try = 5L,
    t_try = 0.2,
    split_structure = "cur_trees_1",
    prune = TRUE
  )
  expect_identical(r$fit$params$split_structure, "cur_trees_1")
  expect_identical(strip_tune(r), strip_tune(p))
  expect_identical(r$fit$fit$get_model(), p$fit$fit$get_model())
})

test_that("the search table is complete: every forest of the three steps, each step's choice, one final", {
  d <- to_sim()
  r <- to_tune(d)
  tab <- r$tune$search
  key <- paste(tab$split_try, tab$t_try, tab$min_leaf_rows)
  expect_false(anyDuplicated(key) > 0)
  # step 1: the searched points at the first m
  expect_length(r$tune$steps$step1$candidates, 2L)
  expect_true(all(tab$min_leaf_rows[tab$step == 1] == 1L))
  expect_setequal(tab$t_try[tab$step == 1], to_ladder$t_try[2:3])
  expect_equal(sum(tab$chosen_step1, na.rm = TRUE), 1L)
  i1 <- which(tab$chosen_step1)
  # step 2: every m at the chosen point
  expect_length(r$tune$steps$step2$candidates, 2L)
  at_i1 <- tab$split_try == tab$split_try[i1] & tab$t_try == tab$t_try[i1]
  expect_setequal(tab$min_leaf_rows[at_i1], c(1L, 20L))
  expect_equal(sum(tab$chosen_step2, na.rm = TRUE), 1L)
  m2 <- tab$min_leaf_rows[which(tab$chosen_step2)]
  # step 3: the chosen point's ladder neighbours at the chosen m are there
  li <- match(paste(tab$split_try[i1], tab$t_try[i1]), paste(to_ladder$split_try, to_ladder$t_try))
  nb <- intersect(c(li - 1L, li + 1L), seq_len(nrow(to_ladder)))
  for (j in nb) {
    expect_true(any(
      tab$split_try == to_ladder$split_try[j] & tab$t_try == to_ladder$t_try[j] & tab$min_leaf_rows == m2
    ))
  }
  expect_equal(sum(tab$chosen_step3, na.rm = TRUE), 1L)
  expect_true(tab$final[which(tab$chosen_step3)])
  expect_equal(sum(tab$final), 1L)
  expect_true(all(is.finite(tab$oob_rmse) & tab$splits >= 1 & tab$stopped_at >= tab$splits))
  # each step chose the lowest score among its candidates (ties: greedier, larger m)
  keys <- paste0(
    match(paste(tab$split_try, tab$t_try), paste(to_ladder$split_try, to_ladder$t_try)),
    "_",
    tab$min_leaf_rows
  )
  for (s in names(r$tune$steps)) {
    st <- r$tune$steps[[s]]
    mse <- tab$oob_rmse[match(st$candidates, keys)]^2
    expect_identical(st$chosen, st$candidates[which(mse == min(mse))[1]])
  }
})

test_that("the choice: the lowest MSE, ties to the simpler", {
  ch <- randomPlantedForest:::tune_oob_choose
  expect_identical(ch(c(2, 1, 3), simpler = 1:3), 2L)
  expect_identical(ch(c(1, 1, 3), simpler = c(2, 1, 3)), 2L)
  expect_identical(ch(c(1, 1 + 1e-14, 3), simpler = c(2, 1, 3)), 2L)
  expect_identical(ch(c(1, 1 + 1e-6, 3), simpler = c(2, 1, 3)), 1L)
})

test_that("the same random number state gives the same search and model", {
  skip_on_cran()
  d <- to_sim()
  a <- to_tune(d, seed = 4)
  b <- to_tune(d, seed = 4)
  cols <- setdiff(names(a$tune$search), "seconds")
  expect_identical(a$tune$search[cols], b$tune$search[cols])
  expect_identical(predict(a, d$x), predict(b, d$x))
  expect_identical(strip_tune(a), strip_tune(b))
})

test_that("TRUE and \"quick\" are the documented controls; the quick preset grows two forests", {
  resolve <- randomPlantedForest:::tune_resolve
  expect_null(resolve(FALSE, NULL))
  expect_identical(resolve(TRUE, NULL), tune_control())
  q <- resolve("quick", NULL)
  expect_identical(q, tune_control(search = c(5, 7), min_leaf_rows = 1, neighbours = FALSE))
  expect_identical(q$ladder$split_try[q$search], c(3L, 1L))
  expect_identical(q$ladder$t_try[q$search], c(0.1, 0.025))
  expect_identical(tune_control()$search, 3:8)
  expect_identical(tune_control()$ladder[3, "t_try"], 0.4)
  expect_identical(tune_control()$ladder, greediness_ladder())

  skip_on_cran()
  d <- to_sim()
  set.seed(2)
  r <- rpf(d$x, d$y, max_interaction = 2, ntrees = 10, tune = "quick")
  tab <- r$tune$search
  expect_identical(tab$t_try, c(0.1, 0.025))
  expect_identical(tab$min_leaf_rows, c(1L, 1L))
  expect_identical(names(r$tune$steps), "step1")
  expect_identical(which(tab$final), which.min(tab$oob_rmse))
})

test_that("the formula interface tunes the same as x/y; print() and marshal show and keep the search", {
  skip_on_cran()
  d <- to_sim()
  df <- data.frame(d$x, y = d$y)
  ctl <- tune_control(ladder = to_ladder, search = 2:3, min_leaf_rows = 1L, neighbours = FALSE)
  set.seed(5)
  a <- rpf(y ~ ., data = df, max_interaction = 2, ntrees = 10, tune = ctl)
  set.seed(5)
  b <- rpf(d$x, d$y, max_interaction = 2, ntrees = 10, tune = ctl)
  cols <- setdiff(names(a$tune$search), "seconds")
  expect_identical(a$tune$search[cols], b$tune$search[cols])
  expect_identical(predict(a, df)$.pred, predict(b, d$x)$.pred)
  out <- capture.output(print(a))
  expect_true(any(grepl("Tuned out of bag", out, fixed = TRUE)))
  restored <- rpf_unmarshal(rpf_marshal(a))
  expect_identical(restored$tune, a$tune)
  expect_identical(predict(restored, df), predict(a, df))
})

test_that("rpf(tune = ) and tune_control() check their arguments", {
  d <- to_sim(n = 100)
  f <- function(...) rpf(d$x, d$y, ntrees = 10, ...)
  expect_error(f(tune = "full"), "tune")
  expect_error(f(tune = NA), "tune")
  expect_error(f(tune = TRUE, split_try = 3), "split_try")
  expect_error(f(tune = TRUE, t_try = 0.1), "t_try")
  expect_error(f(tune = TRUE, min_leaf_rows = 5), "min_leaf_rows")
  expect_error(f(tune = TRUE, splits = 30), "splits")
  expect_error(f(tune = TRUE, prune = FALSE), "prune")
  expect_error(f(tune = TRUE, deterministic = TRUE), "deterministic")
  expect_error(rpf(d$x, d$y, ntrees = 5, tune = TRUE), "at least 10")
  expect_error(rpf(Species ~ ., data = iris, ntrees = 10, tune = TRUE), "regression only")
  expect_error(f(tune = TRUE, split_structure = "hist"), "min_leaf_rows")
  expect_error(f(tune = "quick", oob_path = 1:3), "oob_path")
  expect_error(tune_control(ladder = data.frame(a = 1)), "split_try")
  expect_error(tune_control(ladder = data.frame(split_try = 1, t_try = 0)), "positive")
  expect_error(tune_control(search = 11), "search")
  expect_error(tune_control(search = c(3, 3)), "search")
  expect_error(tune_control(min_leaf_rows = 0), "min_leaf_rows")
  expect_error(tune_control(neighbours = NA), "neighbours")
})
