# The probabilities of multiclass logit from its K - 1 scores W against the first level: p_1 = 1 / (1 + sum(exp(W)))
# and p_k = exp(W_k) p_1. Where a score overflows exp() (above about 709), they are computed with the row's largest
# score taken out, so they stay finite (a pruned model whose rule R refit diverged can have such scores).

softmax_fit <- function() {
  set.seed(3)
  n <- 150
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- factor(sample(c("a", "b", "c"), n, TRUE), levels = c("a", "b", "c"))
  set.seed(4)
  list(x = x, y = y, fit = rpf(x, y, loss = "logit", n_families = 10, n_splits = 5, prune = FALSE))
}

# the reference: the softmax of (0, W), each row's largest entry taken out
softmax_ref <- function(W) {
  t(apply(cbind(0, W), 1, function(s) {
    e <- exp(s - max(s))
    e / sum(e)
  }))
}

test_that("multiclass logit probabilities stay finite where the scores overflow exp()", {
  d <- softmax_fit()
  W <- rbind(c(10, 0), c(700, 0), c(710, 0), c(1000, -5), c(900, 905), c(-1000, -1000), c(2000, 2000))
  P <- as.matrix(randomPlantedForest:::predict_rpf_prob_raw(d$fit, W))
  expect_true(all(is.finite(P)))
  expect_equal(unname(rowSums(P)), rep(1, nrow(W)), tolerance = 1e-14)
  expect_equal(unname(P), softmax_ref(W), tolerance = 1e-12)
  # the two largest equal scores split the probability
  expect_equal(unname(P[7, ]), c(0, 0.5, 0.5))
  expect_equal(unname(P[5, ]), unname(softmax_ref(W[5, , drop = FALSE])[1, ]), tolerance = 1e-12)
})

test_that("rows without overflow keep the plain formula bit for bit", {
  d <- softmax_fit()
  set.seed(1)
  W <- matrix(rnorm(40, sd = 30), 20, 2)
  W[3, ] <- c(800, 1)
  P <- as.matrix(randomPlantedForest:::predict_rpf_prob_raw(d$fit, W))
  ok <- -3
  e <- exp(W[ok, ])
  den <- 1 + rowSums(e)
  expect_identical(unname(P[ok, ]), unname(cbind(1 / den, e / den)))
  expect_true(all(is.finite(P[3, ])))
})

test_that("a pruned model with scores beyond exp()'s range predicts finite probabilities", {
  d <- softmax_fit()
  m <- prune_components(d$fit, d$x, d$y, rule = "D", tau = 0, max_panels = Inf)
  # as a diverged refit of rule R can leave it: an intercept far beyond exp()'s range
  m$intercept <- m$intercept + c(1000, 0)
  P <- as.matrix(predict(m, d$x, type = "prob"))
  expect_true(all(is.finite(P)))
  expect_equal(unname(rowSums(P)), rep(1, nrow(d$x)), tolerance = 1e-14)
  expect_true(all(predict(m, d$x, type = "class")$.pred_class == "b"))
  expect_true(all(P[, 2] > 1 - 1e-12))
})
