library(mvtnorm)
source("rpf.R")
source("predict_rpf.R")
library(randomPlantedForest)
library(microbenchmark)
library(data.table)
library(ggplot2)

# Regression --------------------------------------------------------------
n <- 1000
p <- 4
beta <- c(.5, 1, 0, -1.5)
beta0 <- 0
cov_base <- 0
sigma <- toeplitz(cov_base^(0:(p-1)))
x <- matrix(rmvnorm(n = n, sigma = sigma), ncol = p,
            dimnames = list(NULL, paste0('x', seq_len(p))))
lp <- x %*% beta + beta0
y <- lp[, 1] + rnorm(n)
#dat <- data.frame(y = y, x)

ntrees <- 10
max_interaction <- 4
splits <- 10
split_try <- 10
t_try <- 0.4
loss <- "L2"
delta = 0.1
epsilon = 0.1

# R version
rpf_r <- rpf(Y = y, X = x, ntrees = ntrees, max_interaction = max_interaction, 
                  splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_r <- predict_rpf(x, rpf_r)

misclass_R_reg=mean((lp[,1]-pred_r)^2)

# C++ version
rpf_cpp <- randomPlantedForest::rpf(y = y, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                         splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_cpp <- predict(rpf_cpp, x)

misclass_C_reg=mean((lp[,1]-pred_cpp$.pred)^2)

diff_reg=mean((pred_r-pred_cpp$.pred)^2)

pred_regr <- data.table(Task = "regr", R = pred_r, Cpp = pred_cpp$.pred)

# # Plot predictions
# plot(pred_r, pred_cpp$.pred)
# abline(0, 1, col = "red")
# 
# # Compare runtime
# microbenchmark(r = {rpf::rpf(Y = y, X = x, ntrees = ntrees, max_interaction = max_interaction, 
#                              splits = splits, split_try = split_try, t_try = t_try, loss = loss)}, 
#                cpp = {randomPlantedForest::rpf(y = y, x = x, ntrees = ntrees, max_interaction = max_interaction, 
#                                                splits = splits, split_try = split_try, t_try = t_try, loss = loss)}, 
#                times = 1)

# Classification --------------------------------------------------------------
y <- rbinom(n, 1, plogis(lp[, 1]))
yfac <- factor(y)

lop<- (sign(lp[,1])+1)/2
#dat <- data.frame(y = y, x)

# L1 loss --------------------------------------------------------------
loss <- "L1" # "logit" # "L1" # "L2" # "exponential"

# R version
rpf_r <- rpf(Y = y, X = x, ntrees = ntrees, max_interaction = max_interaction, 
                  splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_r <- predict_rpf(x, rpf_r)
if (loss %in% c("logit", "exponential")) {
  pred_r <- 1 / (1 + exp(-pred_r))
}
pred_r <- pmin(1, pred_r)
pred_r <- pmax(0, pred_r)

misclass_R_L1=mean((lop-pred_r)^2)

# C++ version
rpf_cpp <- randomPlantedForest::rpf(y = yfac, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                                    splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_cpp <- predict(rpf_cpp, x, type = "prob")

misclass_C_L1=mean((lop-pred_cpp$.pred_1)^2)

diff_L1=mean((pred_r-pred_cpp$.pred_1)^2)

pred_classif_L1 <- data.table(Task = "classif_L1", R = pred_r, Cpp = pred_cpp$.pred_1)

# L2 loss --------------------------------------------------------------
loss <- "L2" # "logit" # "L1" # "L2" # "exponential"

# R version
rpf_r <- rpf(Y = y, X = x, ntrees = ntrees, max_interaction = max_interaction, 
                  splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_r <- predict_rpf(x, rpf_r)
if (loss %in% c("logit", "exponential")) {
  pred_r <- 1 / (1 + exp(-pred_r))
}
pred_r <- pmin(1, pred_r)
pred_r <- pmax(0, pred_r)

misclass_R_L2=mean((lop-pred_r)^2)

# C++ version
rpf_cpp <- randomPlantedForest::rpf(y = yfac, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                                    splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_cpp <- predict(rpf_cpp, x, type = "prob")

misclass_C_L2=mean((lop-pred_cpp$.pred_1)^2)

diff_L2=mean((pred_r-pred_cpp$.pred_1)^2)

pred_classif_L2 <- data.table(Task = "classif_L2", R = pred_r, Cpp = pred_cpp$.pred_1)

# Logit loss --------------------------------------------------------------
loss <- "logit" # "logit" # "L1" # "L2" # "exponential"

# R version
rpf_r <- rpf(Y = y, X = x, ntrees = ntrees, max_interaction = max_interaction, 
                  splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_r <- predict_rpf(x, rpf_r)
if (loss %in% c("logit", "exponential")) {
  pred_r <- 1 / (1 + exp(-pred_r))
}
pred_r <- pmin(1, pred_r)
pred_r <- pmax(0, pred_r)

misclass_R_logit=mean((lop-pred_r)^2)

# C++ version

# loss <- "logit_2" # "logit" # "L1" # "L2" # "exponential"

rpf_cpp <- randomPlantedForest::rpf(y = yfac, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                                    splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_cpp <- predict(rpf_cpp, x, type = "prob")

pred_classif_logit <- data.table(Task = "classif_logit", R = pred_r, Cpp = pred_cpp$.pred_1)

diff_logit=mean((pred_r-pred_cpp$.pred_1)^2)

misclass_C_logit=mean((lop-pred_cpp$.pred_1)^2)

# Exponential loss --------------------------------------------------------------
loss <- "exponential" # "logit" # "L1" # "L2" # "exponential"

yexp <- y
yexp[yexp == 0] <- -1

# R version
rpf_r <- rpf(Y = yexp, X = x, ntrees = ntrees, max_interaction = max_interaction, 
                  splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_r <- predict_rpf(x, rpf_r)
if (loss %in% c("logit", "exponential")) {
  pred_r <- 1 / (1 + exp(-pred_r))
}
pred_r <- pmin(1, pred_r)
pred_r <- pmax(0, pred_r)

misclass_R_expo=mean((lop-pred_r)^2)

# C++ version

# loss <- "exponential_2" # "logit" # "L1" # "L2" # "exponential"

rpf_cpp <- randomPlantedForest::rpf(y = yfac, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                                    splits = splits, split_try = split_try, t_try = t_try, loss = loss, delta = delta, epsilon = epsilon)
pred_cpp <- predict(rpf_cpp, x, type = "prob")

misclass_C_expo=mean((lop-pred_cpp$.pred_1)^2)

diff_expo=mean((pred_r-pred_cpp$.pred_1)^2)

pred_classif_exp <- data.table(Task = "classif_exp", R = pred_r, Cpp = pred_cpp$.pred_1)


# Plot --------------------------------------------------------------------
df <- rbind(pred_regr, pred_classif_L1, pred_classif_L2, 
            pred_classif_logit, pred_classif_exp)

ggplot(df, aes(x = R, y = Cpp)) + 
  facet_wrap(~ Task, scales = "free") + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw() #+ 
  #xlim(0,1) + ylim(0, 1)

# Plot distributions
df1 <- df[, .(Task, R)]
df1[, Version := "R"]
df2 <- df[, .(Task, Cpp)]
df2[, Version := "Cpp"]
colnames(df1)[2] <- "Pred"
colnames(df2)[2] <- "Pred"
df_both <- rbind(df1, df2)

ggplot(df_both, aes(x = Pred, col = Version)) + 
  facet_wrap(~Task, scales = "free")+ 
  stat_ecdf()

misclass_R_reg
misclass_C_reg
misclass_R_L2
misclass_C_L2
misclass_R_L1
misclass_C_L1
misclass_R_logit
misclass_C_logit
misclass_R_expo
misclass_C_expo

diff_reg
diff_L2
diff_L1
diff_logit
diff_expo