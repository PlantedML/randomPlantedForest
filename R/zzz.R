Rcpp::loadModule("mod_rpf", TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("randomPlantedForest", libpath) # nocov
}
