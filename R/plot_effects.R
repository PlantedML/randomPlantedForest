#' Plot Prediction Components
#'
#' @rdname plot_components
#' @param components,object Predicted components of an `rpf()` including the original data the model was fit on, as
#'   returned by `extract_components()`
#' @param predictor,predictors `[character]` vector of predictor names, e.g. `"x1"` to plot main effect of `x1`, and
#'   `c("x1", "x2")` to plot the interaction term `x1:x2`.
#' @param ... Unused
#'
#' @return A `ggplot2` object.
#' @import ggplot2
#' @importFrom scales alpha
#' @export
#'
#' @examples
#'
#' # introduce factor variables to show categorical feature handling
#' mtcars$cyl <- factor(mtcars$cyl)
#' mtcars$vs <- factor(mtcars$vs)
#'
#' # Fit forest, extract components
#' set.seed(12)
#' rpfit <- rpf(mpg ~ cyl + wt + hp + drat + vs, data = mtcars, ntrees = 25, max_interaction = 3)
#' components <- extract_components(rpfit, mtcars)
#'
#' # Main effects ----
#' plot_main_effect(components, "wt")
#' plot_main_effect(components, "drat")
#' plot_main_effect(components, "cyl")
#'
#' # 2-degree interaction effects ----
#' # 2d continuous, scatterplot of arbitrary orientation
#' plot_twoway_effects(components, c("wt", "drat"))
#' plot_twoway_effects(components, c("drat", "wt"))
#'
#' # continuous + categorical (forces continuous on x axis, colors by categorical)
#' plot_twoway_effects(components, c("wt", "cyl"))
#' # identical: plot_twoway_effects(components, c("cyl", "wt"))
#'
#' # 2d categorical, heatmap of arbitrary orientation
#' plot_twoway_effects(components, c("vs", "cyl"))
#' plot_twoway_effects(components, c("cyl", "vs"))
plot_main_effect <- function(components, predictor, ...) {
  checkmate::assert_class(components, "rpf_components")
  checkmate::assert_string(predictor) # Must be a single predictor
  checkmate::assert_subset(predictor, names(components$x), empty.ok = FALSE)

  xdf <- assemble_components(components, predictor)
  x_type <- get_x_types(components, predictor)

  if (x_type == "continuous") {
    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[predictor]], y = .data[["m"]])
    )
    p <- p + ggplot2::geom_line(size = 1.2, key_glyph = "rect")

    # p <- p + ggplot2::geom_point()
  }

  if (x_type == "categorical") {
    p <- ggplot2::ggplot(unique(xdf), ggplot2::aes(x = .data[[predictor]], y = .data[["m"]]))
    p <- p + ggplot2::geom_col(alpha = .8, width = 2/3)
    p <- p + theme(panel.grid.major.x = element_blank())
  }

  p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(y = label_m(predictor))
}

#' @rdname plot_components
#' @export
plot_twoway_effects <- function(components, predictors, ...) {
  checkmate::assert_class(components, "rpf_components")
  checkmate::assert_character(predictors, len = 2, unique = TRUE)
  checkmate::assert_subset(predictors, names(components$x))

  x_types <- get_x_types(components, predictors)
  xdf <- assemble_components(components, predictors)
  x_cat <- names(xdf)[which(x_types == "categorical")]
  x_cont <- names(xdf)[which(x_types == "continuous")]

  fillaes <- ggplot2::aes(
    color = .data[["m"]],
    fill = ggplot2::after_scale(.data[["colour"]])
  )

  # 2x continuous ----
  if (setequal(x_types, "continuous")) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cont[[1]]]],
      y = .data[[x_cont[[2]]]],
      !!!fillaes
    ))
    p <- p + ggplot2::geom_point(size = 2.5, shape = 21, stroke = 1)
    p <- p + diverging_palette(limits = get_m_limits(xdf))
    p <- p + ggplot2::labs(color = label_m(predictors))
  }

  # 2x categorical ----
  if (setequal(x_types, "categorical")) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cat[[1]]]], y = .data[[x_cat[[2]]]],
      !!!fillaes
      ))
    p <- p + ggplot2::geom_tile()
    p <- p + diverging_palette(limits = get_m_limits(xdf))
    p <- p + ggplot2::labs(color = label_m(predictors))
  }

  # 1x categorical 1x continuous ----
  if (setequal(x_types, c("categorical", "continuous"))) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
        x = .data[[x_cont]], y = .data[["m"]], color = .data[[x_cat]]
      )) +
      ggplot2::geom_line(size = 1.2, key_glyph = "rect") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::labs(y = label_m(predictors))
  }

  # Final cleanup ----
  p + ggplot2::theme(legend.position = "bottom")
}

#' @rdname plot_components
#' @export
#' @examples
#' # plot_threeway_effects(components, c("hr", "temp", "workingday"))
plot_threeway_effects <- function(components, predictors, ...) {
  checkmate::assert_class(components, "rpf_components")
  checkmate::assert_character(predictors, len = 3, unique = TRUE)
  checkmate::assert_subset(predictors, names(components$x), empty.ok = FALSE)

  # Create look-up table for predictors and their types
  x_types <- get_x_types(components, predictors)
  xdf <- assemble_components(components, predictors)
  x_cat <- names(xdf)[which(x_types == "categorical")]
  x_cont <- names(xdf)[which(x_types == "continuous")]

  # 3x continuous ----
  if (all(x_types == "continuous")) {
    stop("Can't visualize 3 continuous predictor effects (yet?), feel free to make a suggestion!")
  }

  # 1x categorical 2x continuous ----
  if (all(sort(x_types) == c("categorical", "continuous", "continuous"))) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cont[[1]]]],
      y = .data[[x_cont[[2]]]],
      colour = .data[["m"]],
      fill = ggplot2::after_scale(.data[["colour"]])
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[[x_cat]])) +
      ggplot2::geom_point(size = 2.5, shape = 21, stroke = 1) +
      diverging_palette(limits = get_m_limits(xdf)) +
      ggplot2::labs(color = label_m(predictors))
  }

  # 2x categorical 1x continuous ----
  if (all(sort(x_types) == c("categorical", "categorical", "continuous"))) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cont]],
      y = .data[["m"]],
      colour = .data[[x_cat[[1]]]],
      fill = ggplot2::after_scale(.data[["colour"]])
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[[x_cat[[2]]]]), scales = "free_y") +
      ggplot2::geom_line(size = 1.2, key_glyph = "rect") +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::labs(y = label_m(predictors))
  }

  # 3x categorical ----
  if (all(x_types == "categorical")) {

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cat[[1]]]],
      y = .data[[x_cat[[2]]]],
      colour = .data[["m"]],
      fill = ggplot2::after_scale(.data[["colour"]])
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[[x_cat[[3]]]])) +
      ggplot2::geom_tile() +
      diverging_palette(limits = get_m_limits(xdf)) +
      ggplot2::labs(
        color = label_m(predictors)
      )
  }

  # Final cleanup ----
  p + ggplot2::theme(legend.position = "bottom")
}

#' @rdname plot_components
#' @export
autoplot.rpf_components <- function(object, predictors, ...) {
  np <- length(predictors)

  if (np == 1) p <- plot_main_effect(object, predictors, ...)
  if (np == 2) p <- plot_twoway_effects(object, predictors, ...)
  if (np == 3) p <- plot_threeway_effects(object, predictors, ...)

  p
}

#' Assemble original data and corresponding component
#' @keywords internal
#' @noRd
#' @inheritParams plot_twoway_effects
assemble_components <- function(components, predictors) {
  xtemp <- data.table::copy(components)
  xdf <- xtemp$x[, predictors, with = FALSE]
  xdf$m <- xtemp$m[[paste(sort(predictors), collapse = ":")]]
  xdf[]
}

#' Create component label from predictors
#'
#' @keywords internal
#' @noRd
label_m <- function(predictors, mathy = TRUE) {

  preds <- paste0(sort(predictors), collapse = ", ")

  if (mathy) {
    as.expression(bquote(hat(m)[plain(.(preds))]))
  } else {
    sprintf("m(%s)", preds)
  }
}

#' Utility to get symmetric range of component
#' @keywords internal
#' @noRd
get_m_limits <- function(xdf) {
  c(-1,1) * max(abs(xdf[["m"]]))
}

#' Utility to get predictor types from training data
#' @keywords internal
#' @noRd
get_x_types <- function(components, predictors) {
  # Create look-up table for predictors and their types
  tp <- c(numeric = "continuous", integer = "continuous", character = "categorical", factor = "categorical")
  x_types <- vapply(predictors, function(p) {
    cl <- class(components[["x"]][[p]])[1]
    tp[cl]
  }, "")
  checkmate::assert_subset(x_types, c("categorical", "continuous"), empty.ok = FALSE)
  x_types
}

diverging_palette <- function(...) {

  # guide_colorbar <- ggplot2::guide_colorbar(
  #   barwidth = ggplot2::unit(15, "char"),
  #   barheight = ggplot2::unit(1, "char"),
  #   title.position = "right",
  #   title.hjust = .5
  # )
  #
  guide_colorbar <- ggplot2::guide_colorbar(
    barwidth = ggplot2::unit(10.2, "lines"),
    barheight = ggplot2::unit(1, "char"),
    title.position = "bottom",
    title.hjust = .5, title.vjust = 1
  )

  if (!requireNamespace("scico")) {
    ggplot2::scale_color_distiller(
      palette = "PRGn", type = "div",
      guide = guide_colorbar,
      ...
    )
  } else {
    scico::scale_color_scico(
      palette = "vikO",
      guide = guide_colorbar,
      ...
    )
  }

}
