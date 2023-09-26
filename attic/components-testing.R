
library(randomPlantedForest)
set.seed(4)
max_interaction <- 2
object <- rpf(mpg ~ ., data = mtcars, max_interaction = max_interaction)
new_data <- mtcars[, -1]
predictors <- names(new_data)

tictoc::tic()
all_components_orig <- lapply(seq_len(max_interaction), function(i) {
  combinations <- utils::combn(predictors, i, simplify = FALSE)
  components <- lapply(combinations, function(x) .predict_single_component(object, new_data, x))
  do.call(cbind, args = components)
})
all_components_orig <- do.call(cbind, all_components_orig)
tictoc::toc()

tictoc::tic()
all_components_rcpp <- lapply(seq_len(max_interaction), function(i) {
  combinations <- RcppAlgos::comboGeneral(predictors, m = i, nThreads = 1)
  components <- apply(combinations, 1, function(j) .predict_single_component(object, new_data, j))
  do.call(cbind, args = components)
})
all_components_rcpp <- do.call(cbind, all_components_rcpp)
tictoc::toc()

testthat::expect_equal(all_components_orig, all_components_rcpp)



predict_components(rp, mtcars, predictors = c("cyl", "am"))
predict_components_rcppAlgos(rp, mtcars, predictors = c("am", "cyl"))

.predict_single_component(rp, mtcars, predictors = c("cyl", "am"))



rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 10, purify = TRUE)
predict_components(rpfit, mtcars, predictors = c("cyl", "am"))

bench::mark(
  combn = utils::combn(names(mtcars), 5, simplify = FALSE),
  RcppAlgos = RcppAlgos::comboGeneral(names(mtcars), m = 5, nThreads = 1),
  check = FALSE
)


## repr
library(randomPlantedForest)

rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 9, purify = TRUE)
prof <- profvis::profvis(predict_components(rpfit, mtcars))

# "wrong" order in new_data
rp$fit$predict_matrix(as.matrix(mtcars[, c("am", "cyl")]), c(1, 2))
rp$fit$predict_matrix(as.matrix(mtcars[, c("am", "cyl")]), c(2, 1))

# correct order in new_data
rp$fit$predict_matrix(as.matrix(mtcars[, c("cyl", "am")]), c(1, 2))
rp$fit$predict_matrix(as.matrix(mtcars[, c("cyl", "am")]), c(2, 1))

#####

rpfit <- rpf(mpg ~ cyl + am + gear + hp + wt, data = mtcars, max_interaction = 5, purify = TRUE)

predict_components(rpfit, mtcars, predictors = c("cyl", "am"))

bench::mark(
  old = predict_components_old(rpfit, mtcars),
  new = predict_components(rpfit, mtcars),
  rcpp = predict_components_Rcpp(rpfit, mtcars),
  check = TRUE
)


RcppAlgos::comboGeneral(c("cyl", "am", "gear"), 2, FUN = identity, FUN.VALUE = NULL)
max_interaction = 3
expand.grid(
  max_interaction = seq_len(max_interaction),
  combs = c("cyl", "am", "gear")
)
