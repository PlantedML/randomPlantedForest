# rpf() with the defaults of randomPlantedForest 0.5.0, for the tests written for them: a forest (no pruning) of
# order 2 with 50 tree families, 30 splits and 10 / 0.4. Each is added only where the call does not give the argument,
# so a call is evaluated exactly as if it spelled them out (a test that gives `splits = NULL` keeps it). Calls with
# `tune` get nothing added. The draw of candidates and the classification loss are not set here: a test that needs the
# exact draw or `loss = "exponential"` says so. Tests of the current defaults call rpf() itself.
defaults_050 <- list(max_interaction = 2, ntrees = 50, splits = 30, split_try = 10, t_try = 0.4, prune = FALSE)

rpf_050 <- function(...) {
  given <- ...names()
  call <- sys.call()
  call[[1]] <- quote(rpf)
  if (!"tune" %in% given) {
    for (a in setdiff(names(defaults_050), given)) {
      call[[a]] <- defaults_050[[a]]
    }
  }
  eval(call, parent.frame())
}
