# How classif_before_n18.rds was made (test-driftfree-losses.R): three classification fits of the version before the
# drift-free losses (with delta and epsilon in their params), marshalled, with their predictions. Not run by the tests.
#   Rscript make-classif-before-n18.R <library holding that version> <out.rds>
args <- commandArgs(TRUE)
library(randomPlantedForest, lib.loc = args[1])
out <- args[2]
set.seed(1)
n <- 120
d <- data.frame(x1 = round(runif(n), 3), x2 = round(runif(n), 3), x3 = sample(0:2, n, TRUE))
eta <- 2 * d$x1 - 1 + d$x2 * d$x3
d$yb <- factor(ifelse(eta + rnorm(n, sd = 0.3) > 0.3, "hi", "lo"))
d$y3 <- factor(c("a", "b", "c")[1 + (eta > 0) + (eta > 0.8)])
nd <- d[1:12, c("x1", "x2", "x3")]
set.seed(2)
f1 <- rpf(
  yb ~ x1 + x2 + x3,
  data = d,
  loss = "exponential",
  n_families = 3,
  n_splits = 10,
  max_interaction = 2,
  prune = FALSE,
  delta = 0.01,
  epsilon = 0.2
)
set.seed(3)
f2 <- rpf(
  y3 ~ x1 + x2 + x3,
  data = d,
  loss = "logit",
  n_families = 3,
  n_splits = 10,
  max_interaction = 2,
  prune = FALSE
)
set.seed(4)
f3 <- rpf(y3 ~ x1 + x2 + x3, data = d, loss = "exponential", n_families = 10, max_interaction = 2)
stopifnot(identical(f1$params$delta, 0.01), !is.null(f2$params$epsilon), inherits(f3, "rpf_pruned"))
fixture <- list(
  version = as.character(utils::packageVersion("randomPlantedForest")),
  data = d,
  new_data = nd,
  exponential2 = rpf_marshal(f1),
  logit3 = rpf_marshal(f2, include_data = TRUE),
  pruned3 = rpf_marshal(f3),
  pred = list(
    exponential2 = predict(f1, nd, type = "prob"),
    logit3 = predict(f2, nd, type = "prob"),
    logit3_link = predict(f2, nd, type = "link"),
    pruned3 = predict(f3, nd, type = "prob")
  )
)
saveRDS(fixture, out, version = 2)
cat("written", out, file.size(out), "bytes\n")
