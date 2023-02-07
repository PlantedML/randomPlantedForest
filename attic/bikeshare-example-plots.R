# From https://github.com/PlantedML/shap_decomposition/blob/master/bike_example/bike_example.R
# Adapted to plot rpf m's rather than xgb
library(ISLR2)
library(xgboost)
library(randomPlantedForest)
library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(glex)

set.seed(2022)

# Prepare data
data(Bikeshare)
bike <- data.table(Bikeshare)
bike[, hr := as.numeric(as.character(hr))]
# x <- as.matrix(bike[, .(day, hr, temp, windspeed, workingday, hum)])
# y <- bike$bikers

# xgboost
# xg <- xgboost(data = x, label = y, params = list(max_depth = 4, eta = .1), nrounds = 20)
rp <- rpf(bikers ~ day + hr + temp + windspeed + workingday + hum, data = bike, max_interaction = 3, ntrees = 200, parallel = TRUE)
# SHAP decomposition
# res <- glex::glex(xg, x)

# Plot SHAP, main effects, 2-way and 3-way interactions together ----------
vars <- c("hr", "temp", "workingday")

# Main effects

ms <- extract_components(rp, bike, predictors = vars)
ms <- as.data.table(ms)

ms2 <- rbindlist(lapply(vars, function(colname) {
  # have to append cols as vectors ([[1]]) otherwise weird stuff happens
  data.table(variable = colname, value = bike[, ..colname][[1]], m = ms[, ..colname][[1]])
}))


p2 <- ggplot(ms2, aes(x = value, y = m)) +
  facet_wrap(~variable, scales = "free_x") +
  geom_abline(intercept = 0, slope = 0, col = "red") +
  geom_line() +
  theme_bw() +
  xlab("Feature value") +
  ylab(expression(Main~effects~italic(m[j])))

p2

# 2-way interactions
int2way <- rbind(data.table(variable = "hr", value = bike[, "hr"][[1]],
                            workingday = (bike[, "workingday"][[1]]),
                            m = ms[, "hr:workingday"][[1]]),
                 data.table(variable = "temp", value = bike[, "temp"][[1]],
                            workingday = (bike[, "workingday"][[1]]),
                            m = ms[, "temp:workingday"][[1]]))
int2way[, workingday := factor(workingday, levels = 0:1,
                               labels = c("No working day", "Working day"))]

p3_1 <- ggplot(int2way, aes(x = value, y = m, col = workingday)) +
  facet_wrap(~variable, scales = "free_x") +
  geom_abline(intercept = 0, slope = 0, col = "red") +
  geom_line() +
  scale_color_manual(values = viridisLite::viridis(4)[c(2, 3)], name = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Feature value") +
  ylab(expression(2-way~interactions~italic(m[jk])))
p3_2 <- ggplot(data.frame(hr = x[, "hr"], temp = x[, "temp"], m = ms[, "hr:temp"][[1]]), aes(x = hr, y = temp, col = m)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(name = expression(italic(m[jk])))
p3 <- plot_grid(p3_1, p3_2, rel_widths = c(.55, .45))

p3

# 3-way interaction
p4 <- ggplot(data.frame(hr = bike[, "hr"][[1]], temp = bike[, "temp"][[1]],
                        workingday = factor(bike[, "workingday"][[1]], levels = c("0", "1"),
                                            labels = c("No working day", "Working day")),
                        m = ms[, "hr:temp:workingday"][[1]]),
             aes(x = hr, y = temp, col = m)) +
  facet_wrap(~workingday) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(name = expression(italic(m[jkl]))) +
  ylab(expression(atop(3-way~interactions~italic(m[jkl]), "temp")))

p4

plot_grid(p2, p3, p4, ncol = 1)
#ggsave("bike_example.pdf", width = 10, height = 10.5)
#ggsave("bike_example.png", width = 10, height = 10.5)
