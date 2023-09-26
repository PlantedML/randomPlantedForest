library(randomPlantedForest)
# Regression task, select target + 3 predictors
train <-  mtcars[1:20, 1:4]
test <-  mtcars[21:32, 1:4]

set.seed(23)
rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)

# Not purified
rpfit$fit$is_purified()

(components <- predict_components(rpfit, test))

# Auto-purifies
rpfit$fit$is_purified()

cbind(m_sum = rowSums(components), predict(rpfit, test), truth = test$mpg)


grid <- expand.grid(
  p = 1:30,
  max_interaction = 1:30
) |>
  subset(max_interaction <= p)

grid$combs <- mapply(function(p, m) {
  sum(choose(p, seq_len(m)))
}, p = grid$p, m = grid$max_interaction)


library(ggplot2)

ggplot(grid, aes(x = p, y = max_interaction, fill = log10(combs))) +
  geom_rect(aes(xmin = p - 0.5, xmax = p + 0.5, ymin = max_interaction - 0.5, ymax = max_interaction + 0.5)) +
  scale_fill_viridis_b(breaks = seq(0, 10)) +
  theme_minimal()

ggplot(grid, aes(x = p, y = combs, fill = max_interaction)) +
  geom_path(aes(group = max_interaction)) +
  geom_point(shape = 22, size = 4) +
  #scale_y_log10() +
  scale_fill_viridis_b(breaks = seq(0, 30, 5)) +
  theme_minimal()


# binary ----------------------------------------------------------------------------------------------------------

set.seed(124)
xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2),
  x5 = factor(sample(letters[1:4], 100, replace = TRUE))
)

bin_fit <- rpf(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)
bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")

#.predict_single_component(bin_fit, xdat, "x2")
m <- predict_components(bin_fit, xdat)
cbind(msum = rowSums(m), bin_pred)

# multiclass ------------------------------------------------------------------------------------------------------

library(randomPlantedForest)
library(palmerpenguins)

penguins <- na.omit(penguins)

set.seed(23)
rpfit <- rpf(species ~ island + bill_length_mm + bill_depth_mm,
             data = penguins, max_interaction = 2, ntrees = 30)
purify(rpfit)

components <- predict_components(rpfit, penguins)
str(components)

# numeric predictions: raw, untransformed from C++
preds <- predict(rpfit, penguins, type = "numeric")

cbind(
  m_sum = rowSums(components),
  preds,
  predsum = rowSums(preds),
  truth = penguins$species
) |> head()

# Transform such that we have 1 row per level
components_multi_long <- do.call(rbind, lapply(levels(penguins$species), function(x) {
  # Extract components for given y level
  dat <- components[, endsWith(names(components), x), drop = FALSE]
  names(dat) <- gsub(paste0("_", x), "", names(dat))
  # Append y level (needs check for pre-existing column with that name!)
  dat[["y_level"]] <- x
  # Reorder such that y_level comes first
  dat <- dat[, c(ncol(dat), seq_len(ncol(dat) - 1))]
  dat
}))

head(components_multi_long)


cbind(
  y_level = components_multi_long$y_level,
  msum = rowSums(components_multi_long[2:ncol(components_multi_long)]),
  preds,
  predsum = rowSums(preds)
)



components |>
  dplyr::mutate(id = seq_len(nrow(components))) |>
  tidyr::pivot_longer(cols = seq_len(ncol(components))) |>
  tidyr::separate(col = "name", into = c("m", "y_level"), sep = "(_)(?!.*_)") |>
  dplyr::group_by(m, id) |>
  dplyr::summarise(msum = sum(value)) |>
  dplyr::arrange(id, m)
