# Sanity checks in the spirit of mlr3's autotest(): any configuration must learn
# an outcome fully determined by a single discrete predictor. The signal sits on
# the lowest, middle or top value, which need the lowest threshold, two splits
# on the same variable, or the highest threshold, respectively.

sanity_predictors <- function(n = 60) {
  binary <- rep(c(FALSE, TRUE), length.out = n)
  three <- rep(0:2, length.out = n)
  list(
    numeric = list(x = as.numeric(binary), signal = binary),
    factor = list(x = factor(binary), signal = binary),
    logical = list(x = binary, signal = binary),
    three_low = list(x = three, signal = three == 0),
    three_mid = list(x = three, signal = three == 1),
    three_high = list(x = three, signal = three == 2)
  )
}

# Fits y ~ x for every predictor and split structure, returns failing labels
sanity_failures <- function(make_y, check, ...) {
  failures <- character()
  for (name in names(sanity_predictors())) {
    predictor <- sanity_predictors()[[name]]
    dat <- data.frame(x = predictor$x, y = make_y(predictor$signal))
    for (structure in c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees")) {
      set.seed(1)
      fit <- rpf(y ~ x, data = dat, ntrees = 5, splits = 10, split_structure = structure, ...)
      if (!check(fit, dat, predictor$signal)) {
        failures <- c(failures, paste(name, structure))
      }
    }
  }
  failures
}

test_that("regression learns from a single discrete predictor", {
  failures <- sanity_failures(
    make_y = \(signal) signal + 0.5,
    check = \(fit, dat, signal) {
      pred <- predict(fit, dat)$.pred
      min(pred[signal]) > max(pred[!signal])
    }
  )
  expect_identical(failures, character())
})

test_that("classification learns from a single discrete predictor", {
  classes <- list(
    binary = \(signal) factor(ifelse(signal, "a", "b")),
    multiclass = \(signal) factor(ifelse(signal, "a", rep(c("b", "c"), length.out = length(signal))))
  )
  for (task in names(classes)) {
    for (loss in c("L1", "L2", "logit", "exponential")) {
      failures <- sanity_failures(
        make_y = classes[[task]],
        check = \(fit, dat, signal) {
          pred <- predict(fit, dat, type = "class")$.pred_class
          all((pred == "a") == signal)
        },
        loss = loss
      )
      expect_identical(failures, character(), label = paste(task, loss))
    }
  }
})
