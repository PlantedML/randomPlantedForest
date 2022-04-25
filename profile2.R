
library(randomPlantedForest)
library(mvtnorm)

n <- 1000
p <- 4
beta <- c(1, 2, 0, -3)
beta0 <- 0
cov_base <- 0
sigma <- toeplitz(cov_base^(0:(p-1)))
x <- matrix(rmvnorm(n = n, sigma = sigma), ncol = p,
            dimnames = list(NULL, paste0('x', seq_len(p))))
lp <- x %*% beta + beta0
y <- rbinom(n, 1, plogis(lp[, 1]))
yfac <- factor(y)

ntrees <- 10
max_interaction <- 4
splits <- 30
split_try <- 10
t_try <- 0.4
loss <- "L2" # "logit" # "L1" # "L2" # "exponential"

randomPlantedForest:::start_profiler("/tmp/profile.out")
rpf_cpp <- randomPlantedForest::rpf(y = yfac, x = x, ntrees = ntrees, max_interaction = max_interaction, 
                                    splits = splits, split_try = split_try, t_try = t_try, loss = loss)
randomPlantedForest:::stop_profiler()

system2("/home/wright/go/bin/pprof", 
        "-pdf src/randomPlantedForest.so /tmp/profile.out", 
        stdout = "profile.pdf")
