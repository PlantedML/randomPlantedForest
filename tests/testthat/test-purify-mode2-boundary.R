# purify(mode = 2) must count the training rows that lie exactly on a leaf's lower
# bound. Split points are data values, so with integer-valued features many rows
# sit on a bound. Such a row used to be lost when it held the largest value of its
# KD-tree node (up to 32 rows), so the data have many ties and several nodes.

boundary_data <- function(n = 400) {
  set.seed(4302)
  x <- data.frame(
    x1 = sample(0:9, n, replace = TRUE),
    x2 = sample(0:9, n, replace = TRUE),
    x3 = sample(0:4, n, replace = TRUE),
    x4 = sample(0:9, n, replace = TRUE)
  )
  x$y <- x$x1 + 0.5 * x$x2 * x$x3 - 0.3 * x$x1 * x$x4 + 0.2 * x$x2 * x$x3 * x$x4 + rnorm(n)
  x
}

fit_boundary <- function(dat, max_interaction) {
  set.seed(43)
  rpf(y ~ ., data = dat, max_interaction = max_interaction, ntrees = 10, splits = 30)
}

for (order in 2:3) {
  test_that(paste0("purify mode 2 equals mode 1 with tied features, order ", order), {
    dat <- boundary_data()
    fit1 <- fit_boundary(dat, order)
    fit2 <- fit_boundary(dat, order)
    fit3 <- fit_boundary(dat, order)

    purify(fit1, mode = 1L)
    purify(fit2, mode = 2L)
    m1 <- predict_components(fit1, dat)
    m2 <- predict_components(fit2, dat)

    expect_equal(colnames(m1$m), colnames(m2$m))
    expect_equal(as.matrix(m2$m), as.matrix(m1$m), tolerance = 1e-10)
    expect_equal(m2$intercept, m1$intercept, tolerance = 1e-10)

    # predict_components() purifies an unpurified fit itself, with mode 2
    expect_false(is_purified(fit3))
    m3 <- predict_components(fit3, dat)
    expect_equal(as.matrix(m3$m), as.matrix(m1$m), tolerance = 1e-10)
    expect_equal(m3$intercept, m1$intercept, tolerance = 1e-10)
  })
}
