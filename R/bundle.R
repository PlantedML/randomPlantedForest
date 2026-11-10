#' Bundle an rpf model
#'
#' Method for [bundle::bundle()] wrapping [rpf_marshal()]/[rpf_unmarshal()],
#' so rpf models work with the standard tidymodels serialization workflow.
#' Training data is not included; see [rpf_marshal()] for the implications.
#' Pruned forests from [prune_components()] or [rpf()] are bundled the same way.
#'
#' @param x `[rpf | rpf_pruned]`: A fitted [`rpf`] model, or a pruned one from
#'   [prune_components()].
#' @param ... Unused.
#' @return An object of class `bundled_rpf` (`bundled_rpf_pruned`) for `bundle()`,
#'   or the restored [rpf] (`rpf_pruned`) object for `unbundle()`.
#' @examplesIf requireNamespace("bundle", quietly = TRUE)
#' fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10)
#' b <- bundle::bundle(fit)
#' tmp <- tempfile(fileext = ".rds")
#' saveRDS(b, tmp)
#' restored <- bundle::unbundle(readRDS(tmp))
#' predict(restored, mtcars)
#' @exportS3Method bundle::bundle
bundle.rpf <- function(x, ...) {
  bundle::bundle_constr(
    object = rpf_marshal(x),
    situate = bundle::situate_constr(function(object) randomPlantedForest::rpf_unmarshal(object)),
    desc_class = "rpf"
  )
}

#' @rdname bundle.rpf
#' @exportS3Method bundle::bundle
bundle.rpf_pruned <- function(x, ...) {
  bundle::bundle_constr(
    object = rpf_marshal(x),
    situate = bundle::situate_constr(function(object) randomPlantedForest::rpf_unmarshal(object)),
    desc_class = "rpf_pruned"
  )
}
