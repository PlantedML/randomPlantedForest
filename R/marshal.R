#' Serialize and restore Random Planted Forests
#'
#' `rpf_marshal()` converts an [rpf] object into a plain R list safe for
#' [saveRDS()]; `rpf_unmarshal()` reverses it. The C++ forest behind an `rpf`
#' object does not survive R serialization, so this explicit round-trip is
#' required to store or transfer fitted models.
#'
#' Training data is only included with `include_data = TRUE`; without it, a
#' restored forest can predict (including purified prediction if the forest
#' was purified before marshaling) but can never be purified afterwards.
#'
#' @param x `[rpf]`: A fitted [`rpf`] model.
#' @param include_data `[logical(1): FALSE]`: Store training data in the blob, enabling
#'   [purify()] after restoring.
#' @param blob `[rpf_marshaled]`: Object created by `rpf_marshal()`.
#' @return `rpf_marshal()` returns a list of class `rpf_marshaled`;
#'   `rpf_unmarshal()` returns an object of class [rpf].
#' @export
#' @examples
#' fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10)
#' blob <- rpf_marshal(fit)
#' tmp <- tempfile(fileext = ".rds")
#' saveRDS(blob, tmp)
#' restored <- rpf_unmarshal(readRDS(tmp))
#' all.equal(predict(fit, mtcars), predict(restored, mtcars))
rpf_marshal <- function(x, include_data = FALSE) {
  checkmate::assert_class(x, "rpf")
  checkmate::assert_flag(include_data)
  check_rpf_alive(x)

  purified <- x$fit$is_purified()
  fit_state <- list(
    version = 1L,
    pkg_version = utils::packageVersion("randomPlantedForest"),
    model = x$fit$get_model(),
    bounds = x$fit$get_bounds(),
    dims = x$fit$get_shape(),
    purified = purified,
    grid = if (purified) x$fit$get_grid_leaves(),
    data = if (include_data) x$fit$get_data(),
    # $forest (export_forest = TRUE) duplicates $model; rebuilt on unmarshal
    had_forest = !is.null(x$forest)
  )

  out <- unclass(x)
  out$fit <- NULL
  out$forest <- NULL
  out$fit_state <- fit_state
  structure(out, class = "rpf_marshaled")
}

#' @rdname rpf_marshal
#' @export
rpf_unmarshal <- function(blob) {
  checkmate::assert_class(blob, "rpf_marshaled")
  state <- blob$fit_state
  if (!identical(state$version, 1L)) {
    cli::cli_abort("Unsupported {.cls rpf_marshaled} version {.val {state$version}}.")
  }
  installed <- utils::packageVersion("randomPlantedForest")
  if (!is.null(state$pkg_version) && state$pkg_version > installed) {
    cli::cli_warn(c(
      "This model was marshaled with {.pkg randomPlantedForest} {state$pkg_version}, but {installed} is installed.",
      "!" = "Restoring may not be reliable."
    ))
  }

  pars <- rpf_param_vector(blob$params, blob$mode)
  fit <- if (blob$mode == "classification") {
    methods::new(ClassificationRPF, blob$params$loss, pars)
  } else {
    methods::new(RandomPlantedForest, pars)
  }
  fit$set_shape(
    state$dims$feature_size,
    state$dims$value_size,
    state$dims$sample_size,
    state$bounds$lower,
    state$bounds$upper
  )
  if (!is.null(state$data)) {
    fit$set_training_data(as.matrix(state$data$Y), as.matrix(state$data$X))
  }
  fit$set_model(state$model)
  if (state$purified) {
    fit$set_grid_leaves(state$grid)
  }

  forest <- NULL
  if (isTRUE(state$had_forest)) {
    forest <- fit$get_model()
    class(forest) <- "rpf_forest"
  }

  out <- blob
  out$fit_state <- NULL
  out <- unclass(out)
  fields <- out[setdiff(names(out), c("fit", "blueprint", "forest"))]
  do.call(new_rpf, c(list(fit = fit, blueprint = blob$blueprint, forest = forest), fields))
}

#' Check whether an rpf object's C++ forest is still alive
#'
#' The C++ forest does not survive [saveRDS()]; an `rpf` object restored via
#' [readRDS()] without [rpf_marshal()]/[rpf_unmarshal()] is unusable.
#' @param x `[rpf]`: A fitted [`rpf`] model, possibly restored from disk.
#' @return `TRUE` if the underlying model can be used, `FALSE` otherwise.
#' @export
rpf_is_valid <- function(x) {
  checkmate::assert_class(x, "rpf")
  !inherits(try(x$fit$is_purified(), silent = TRUE), "try-error")
}

check_rpf_alive <- function(x) {
  if (!rpf_is_valid(x)) {
    cli::cli_abort(
      c(
        "The C++ forest behind this {.cls rpf} object is gone.",
        "i" = "Most likely it was saved with {.fn saveRDS} and restored with {.fn readRDS}.",
        "i" = "Use {.code blob <- rpf_marshal(x)} before saving and {.code rpf_unmarshal(blob)} after loading, see {.help rpf_marshal}."
      ),
      call = rlang::caller_env()
    )
  }
  invisible(x)
}
