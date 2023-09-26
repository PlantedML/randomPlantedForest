xdf <- na.omit(palmerpenguins::penguins)

rpfitlogit <- randomPlantedForest::rpf(species ~ ., data = xdf, loss = "logit", ntrees = 50, max_interaction = 3)
predict(rpfitlogit, xdf, type = "prob") |>
  cbind(xdf[, "species"]) |>
  head(20)

# Without logistic transformation = raw predictions
predict(rpfitlogit, xdf, type = "numeric") |>
        cbind(xdf[, "species"]) |>
        head(20)

# Exponential
rpfitexp <- randomPlantedForest::rpf(species ~ ., data = xdf, loss = "exponential", ntrees = 50, max_interaction = 3)
predict(rpfitexp, xdf, type = "prob") |>
  cbind(xdf[, "species"]) |>
  head(20)

# Without logistic transformation = raw predictions
predict(rpfitexp, xdf, type = "numeric") |>
  cbind(xdf[, "species"]) |>
  head(20)

# L1
rpfitl1 <- randomPlantedForest::rpf(species ~ ., data = xdf, loss = "L1", ntrees = 50, max_interaction = 3)
predict(rpfitl1, xdf, type = "prob") |>
  cbind(xdf[, "species"]) |>
  head(20)

# Without logistic transformation = raw predictions
predict(rpfitl1, xdf, type = "numeric") |>
  cbind(xdf[, "species"]) |>
  head(20)
