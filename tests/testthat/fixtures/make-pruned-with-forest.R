# How pruned_with_forest.rds was made (test-keep-forest.R): pruned models of the version before keep_forest (D2,
# c0c34ea8), which kept their forest, marshalled, with their predictions, components and print then. Not run by the
# tests.
#   Rscript make-pruned-with-forest.R <library holding that version> <out.rds>
args <- commandArgs(TRUE)
library(randomPlantedForest, lib.loc = args[1])
out <- args[2]
set.seed(1)
n <- 300
d <- data.frame(
  x1 = round(runif(n), 2),
  x2 = round(runif(n), 2),
  x3 = round(runif(n), 2),
  x4 = round(runif(n), 2),
  x5 = sample(0:2, n, TRUE)
)
eta <- 2 * d$x1 - 1 + 4 * d$x2 * d$x3 * d$x4
d$y <- eta + rnorm(n, sd = 0.1)
d$yb <- factor(ifelse(eta + rnorm(n, sd = 0.3) > 0.3, "hi", "lo"))
nd <- d[1:12, c("x1", "x2", "x3", "x4", "x5")]
# regression, order 3, with the display at 3 heat maps (a coarsened three-way component)
set.seed(2)
reg <- rpf(y ~ x1 + x2 + x3 + x4 + x5, data = d, n_families = 12, max_panels = 3)
# binary classification, the default loss (logit)
set.seed(3)
cls <- rpf(yb ~ x1 + x2 + x3 + x4 + x5, data = d, n_families = 10, max_interaction = 2)
stopifnot(inherits(reg, "rpf_pruned"), inherits(cls, "rpf_pruned"), length(reg$display$coarse) > 0)
fixture <- list(
  version = as.character(utils::packageVersion("randomPlantedForest")),
  data = d,
  new_data = nd,
  reg = rpf_marshal(reg),
  cls = rpf_marshal(cls),
  pred = list(
    reg = predict(reg, nd),
    cls = predict(cls, nd, type = "prob"),
    cls_link = predict(cls, nd, type = "link"),
    cls_class = predict(cls, nd, type = "class")
  ),
  components = list(reg = predict_components(reg, nd)$m, cls = predict_components(cls, nd)$m),
  print = list(reg = format(reg), cls = format(cls))
)
saveRDS(fixture, out, version = 2)
cat("written", out, file.size(out), "bytes\n")
