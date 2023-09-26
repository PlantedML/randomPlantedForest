set.seed(23)
rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 10, ntrees = 30, purify = TRUE)

# Extract all components, including main effects and interaction terms up to `max_interaction`

components_orig <- predict_components(rpfit, mtcars)
components_new <- predict_components_int(rpfit, mtcars)

setequal(names(components_orig$m), names(components_new$m))

names(components_orig$m)
names(components_new$m)

testthat::expect_equal(components_orig$m, components_new$m, list_as_map = TRUE)
testthat::expect_equal(components_orig$intercept, components_new$intercept)
testthat::expect_equal(components_orig$x, components_new$x)

doParallel::registerDoParallel()

bm <- bench::mark(
  orig = predict_components(rpfit, mtcars),
  new = predict_components_int(rpfit, mtcars),
  check = function(x1, x2) {
    stopifnot(setequal(names(x1$m), names(x2$m)))
    testthat::expect_equal(x1$m, x2$m, list_as_map = TRUE)
    testthat::expect_equal(x1$intercept, x2$intercept)
    testthat::expect_equal(x1$x, x1$x)
    return(TRUE)
  },
  min_iterations = 3,
  memory = FALSE
)
bm
plot(bm)

profvis::profvis({predict_components(rpfit, mtcars)})
profvis::profvis({predict_components_int(rpfit, mtcars)})

# ----
library(randomPlantedForest)
library(data.table)
data(Bikeshare, package = "ISLR2")
bike <- data.table(Bikeshare)
bike[, hr := as.numeric(as.character(hr))]
bike[, workingday := factor(workingday, levels = c(0, 1), labels = c("No Workingday", "Workingday"))]
bike[, season := factor(season, levels = 1:4, labels = c("Winter", "Spring", "Summer", "Fall"))]
bike <- bike[weathersit != "heavy rain/snow", ]

rp <- rpf(
  bikers ~ day + hr + temp + windspeed + workingday + hum + weathersit + season,
  data = bike,
  max_interaction = 5, ntrees = 30, purify = TRUE, parallel = TRUE
)

bm <- bench::mark(
  orig = predict_components(rp, bike),
  new = predict_components_int(rp, bike),
  check = function(x1, x2) {
    stopifnot(setequal(names(x1$m), names(x2$m)))
    testthat::expect_equal(x1$m, x2$m, list_as_map = TRUE)
    testthat::expect_equal(x1$intercept, x2$intercept)
    testthat::expect_equal(x1$x, x1$x)
    return(TRUE)
  },
  min_iterations = 3
)
bm
plot(bm)

profvis::profvis({predict_components_int(rp, bike)})



if (foreach::getDoParRegistered()) {
  m_all <- foreach(j = idx, .combine = "+") %dopar% tree_fun(j)
} else {
  m_all <- foreach(j = idx, .combine = "+") %do% tree_fun(j)
}

# mlbench? --------------------------------------------------------------------------------------------------------

xdat <- as.data.frame(mlbench::mlbench.peak(100, d = 20))
rp <- rpf(y ~ ., data = xdat, purify = TRUE, max_interaction = 10, parallel = TRUE)

comp <- predict_components_int(rp, xdat)
vi <- glex::glex_vi(comp)
ggplot2::autoplot(vi, max_interaction = 4, threshold = .01, scale = "relative")
ggplot2::autoplot(vi, by_degree = TRUE, scale = "relative")

RcppAlgos::comboGeneral(5, 3, FUN = c)
sum(sapply(1:10, \(i) RcppAlgos::comboCount(20, i)))

sum(choose(20, seq_len(3)))
