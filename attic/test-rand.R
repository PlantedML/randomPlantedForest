library(randomPlantedForest)


#' Simple vector shuffling using R's rng
#' @param input A numeric vector.
#' @return `input`, but shuffled randomly.
#' @export
shuffle_R_vector <- function(input) {
  .Call(`_randomPlantedForest_shuffle_R_vector`, input)
}


for (i in sample(1e5, 10)) {
  vec <- sample(i)
  randomPlantedForest:::shuffle_R_vector(vec) |>
    unique() |>
    length() |>
    testthat::expect_equal(i)
}



set.seed(12)
vec <- 1:10

(shuffle1 <- randomPlantedForest:::shuffle_R_vector(vec))
(shuffle2 <- randomPlantedForest:::shuffle_R_vector(vec))

testthat::expect_failure(testthat::expect_equal(shuffle1, shuffle2))

set.seed(12)
(shuffle_reset1 <- randomPlantedForest:::shuffle_R_vector(vec))
(shuffle_reset2 <- randomPlantedForest:::shuffle_R_vector(vec))

testthat::expect_failure(testthat::expect_equal(shuffle_reset1, shuffle_reset2))

testthat::expect_equal(shuffle1, shuffle_reset1)
testthat::expect_equal(shuffle2, shuffle_reset2)


# Debug -----------------------------------------------------------------------------------------------------------

library(Rcpp)

set.seed(10)
sourceCpp("src/shuffle_debug.cpp")
sourceCpp("src/shuffle_debug.cpp")


set.seed(10)
sourceCpp("src/shuffle_debug.cpp")
