xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100)
)

test_that("Probability prediction: y = factor / L2 loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  
  bin_pred <- predict(bin_fit_L2, new_data = xdat, type = "prob")
})

test_that("Class prediction: y = factor / L2 loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")
  
})

test_that("Class prediction: y = 0,1 / L2 loss", {
  bin_fit <- suppressWarnings(rpf(y01 ~ x1 + x2, data = xdat, loss = "L2"))
  
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")
  
})
