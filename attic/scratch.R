library(randomPlantedForest)
train <-  tibble::tibble(
  x1 = rnorm(100, 30, 2),
  x2 = rnorm(100, -50, 2),
  x3 = factor(rbinom(100, 2, 1/2), labels = LETTERS[1:3]),
  y = 1.5 * x1 + -3 * x2 + 2 * as.numeric(x3) + rnorm(100)
)


purrr::map_if(train, is.numeric, range)

set.seed(23)
rpfit <- rpf(y ~ ., data = train, max_interaction = 1, ntrees = 2, splits = 5)


mod <- rpfit$fit$get_model()

tree1 <- mod[[1]]
tree1$variables[[2]]
tree1$values[[2]]
tree1$intervals[[2]]

lapply(tree1, function(x) x[[2]] |> length())

names(rpfit$blueprint$ptypes$predictors)
levels(rpfit$blueprint$ptypes$predictors$x3)
rpfit$factor_levels[["x3"]]

lm(y ~ ., data = train) |> summary()


# new_purify crash ------------------------------------------------------------------------------------------------

library(randomPlantedForest)
train <-  mtcars[1:20, ]
test <-  mtcars[21:32, ]

set.seed(23)
rpfit <- rpf(mpg ~., data = train, max_interaction = 3)
pred <- predict(rpfit, test)

# This probably does something
rpfit$fit$new_purify()

# This reliably crashes RStudio
rpfit$fit$purify()

# As well as this
rpfit <- rpf(mpg ~., data = train, purify = TRUE)


# Test purify -----------------------------------------------------------------------------------------------------

library(randomPlantedForest)
train <-  mtcars[1:20, ]
test <-  mtcars[21:32, ]

set.seed(23)
rpfit <- rpf(mpg ~., data = train, max_interaction = 3)
pred <- predict(rpfit, test)

mod1 <- rpfit$fit$get_model()

# This probably does something
rpfit$fit$new_purify()

mod2 <- rpfit$fit$get_model()

identical(mod1, mod2)


# MSE -------------------------------------------------------------------------------------------------------------

library(randomPlantedForest)
train <-  mtcars[1:20, ]
test <-  mtcars[21:32, ]

ntrees <- 35
set.seed(23)
rpfit <- rpf(mpg ~., data = train, max_interaction = 3, ntrees = ntrees)
pred <- predict(rpfit, test)

mse <- rpfit$fit$MSE(as.matrix(pred$.pred), as.matrix(test$mpg))
mse


# Component prediction --------------------------------------------------------------------------------------------
library(randomPlantedForest)
train <-  mtcars[1:20, ]
test <-  mtcars[21:32, ]

ntrees <- 35
set.seed(23)
rpfit <- rpf(mpg ~., data = train, max_interaction = 3, ntrees = ntrees)

# Default behavior, "normal" prediction
pred0 <- predict(rpfit, test)
pred0

# 10 predictors, so components must be vector of length 10
predict(rpfit, test, type = "numeric", components = 1:10)


# Factor levels ---------------------------------------------------------------------------------------------------

rpfit <- rpf(species ~ ., data = na.omit(palmerpenguins::penguins))
levels(palmerpenguins::penguins$sex)
levels(rpfit$blueprint$ptypes$predictors$sex)
rpfit$factor_levels$sex

levels(palmerpenguins::penguins$island)
levels(rpfit$blueprint$ptypes$predictors$island)
rpfit$factor_levels$island


# Weird types -----------------------------------------------------------------------------------------------------

library(randomPlantedForest)
train <-  tibble::tibble(
  x1 = rnorm(100, 30, 2),
  x2 = rnorm(100, -50, 2),
  x3 = factor(rbinom(100, 2, 1/2), labels = LETTERS[1:3]),
  x4 = rep(Sys.time(), 100),
  y = 1.5 * x1 + -3 * x2 + 2 * as.numeric(x3) + rnorm(100)
)


purrr::map_if(train, is.numeric, range)

set.seed(23)
rpfit <- rpf(y ~ ., data = train, max_interaction = 1, ntrees = 2, splits = 5)

pred0 <- predict(rpfit, train)



# Multiclass probability sums --------------------------------------------------------------------------------------


library(randomPlantedForest)
xdf <- na.omit(palmerpenguins::penguins)
train_ids <- sample(333, 222)
train <- xdf[train_ids, ]
test <- xdf[setdiff(1:333, train_ids), ]

rpfit <- randomPlantedForest::rpf(species ~ ., data = train, loss = "L1")
pred <- predict(rpfit, test, type = "prob")
pred$sum <- rowSums(pred)

# Implausible cases are returned as-is:
pred[pred$sum > 1 + 10 * .Machine$double.eps, ]



train_mat <- randomPlantedForest:::preprocess_predictors_predict(rpfit, test[ , -1])
pred_raw <- rpfit$fit$predict_matrix(train_mat, 0)

# What we did so far
pred_trunc <- apply(pred_raw, 2, function(col) pmax(0, pmin(1, col)))
# what we do for exp and logit
pred_softmx <- randomPlantedForest:::softmax(pred_raw)

pred_raw[1:10, ]
pred_trunc[1:10, ]
pred_softmx[1:10, ]

all(pred_raw[, 1] > 1)
all(pred_raw[, 2] < 0)
all(pred_raw[, 3] < 0)

pred_raw - randomPlantedForest:::softmax(pred_raw)



as.data.frame(pred_raw) |>
  tidyr::pivot_longer(cols = 1:3) |>
  # dplyr::filter(name == "V1") |>
  ggplot2::ggplot(ggplot2::aes(x = value)) +
  ggplot2::facet_wrap(~name) +
  ggplot2::geom_histogram()


# reprex
library(randomPlantedForest)
# real data, train/test split
xdf <- na.omit(palmerpenguins::penguins)
train_ids <- sample(333, 222)
train <- xdf[train_ids, ]
test <- xdf[setdiff(1:333, train_ids), ]

# fit
rpfit <- randomPlantedForest::rpf(species ~ ., data = train, loss = "L2")

# needs preprocessing due do factor variables etc.
test_mat <- randomPlantedForest:::preprocess_predictors_predict(rpfit, test[ , -1])

# Raw predictions from C++
pred_raw <- rpfit$fit$predict_matrix(test_mat, 0)

pred_raw[1:10, ]

# All predicted values are either greater 1 or smaller 0
all(pred_raw[, 1] > 1)
all(pred_raw[, 2] < 0)
all(pred_raw[, 3] < 0)

#
predict(rpfit, test, type = "numeric", components = 1)
rpfit$fit$predict_matrix()
