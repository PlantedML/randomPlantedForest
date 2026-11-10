set.seed(275)
link_dat <- data.frame(
  y2 = factor(sample(c("a", "b"), 120, replace = TRUE)),
  y3 = factor(sample(c("hi", "mid", "lo"), 120, replace = TRUE)),
  x1 = rnorm(120),
  x2 = rnorm(120),
  x3 = rnorm(120)
)

test_that("type = \"link\" equals type = \"numeric\" for every classification loss", {
  for (outcome in c("y2", "y3")) {
    for (loss in c("L1", "L2", "logit", "exponential")) {
      fit <- rpf_050(
        stats::reformulate(c("x1", "x2", "x3"), outcome),
        data = link_dat,
        n_families = 5,
        max_interaction = 2,
        n_splits = 10,
        loss = loss,
        nthreads = 1
      )
      link <- predict(fit, link_dat, type = "link")
      num <- predict(fit, link_dat, type = "numeric")
      expect_identical(link, num)
      expect_equal(nrow(link), nrow(link_dat))
    }
  }
})
