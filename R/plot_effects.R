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

  xdf <- data.table::data.table(
    x = components$x[[predictor]],
    m = components$m[[predictor]],
    variable = predictor
  )

  x_type <- switch(class(xdf[["x"]]),
                   numeric = "continuous",
                   integer = "continuous",
                   factor = "categorical",
                   character = "categorical",
                   stop("Can't guess type of predictor: ", class(xdf[["x"]]))
  )

  p <- ggplot2::ggplot(xdf, ggplot2::aes(x = .data[["x"]], y = .data[["m"]]))

  if (x_type == "continuous") {
    p <- p + ggplot2::geom_line()
    p <- p + ggplot2::geom_point()
  } else if (x_type == "categorical") {
    p <- p + ggplot2::geom_col(position = "dodge", alpha = 1/3, width = 2/3)
  }

  p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      x = predictor, y = sprintf("m(%s)", predictor)
    ) +
    ggplot2::theme_minimal()
}


#' @rdname plot_components
#' @export
plot_twoway_effects <- function(components, predictors, ...) {
  checkmate::assert_class(components, "rpf_components")
  checkmate::assert_character(predictors, len = 2, unique = TRUE)
  checkmate::assert_subset(predictors, names(components$x))

  # Create look-up table for predictors and their types
  tp <- c(numeric = "continuous", integer = "continuous", character = "categorical", factor = "categorical")
  x_types <- vapply(predictors, function(p) {
    cl <- class(components[["x"]][[p]])[1]
    tp[cl]
  }, "")
  checkmate::assert_subset(x_types, c("categorical", "continuous"), empty.ok = FALSE)

  xdf <- data.table::data.table(
    x = components[["x"]][[predictors[[1]]]],
    y = components[["x"]][[predictors[[2]]]],
    m = components[["m"]][[paste0(sort(predictors), collapse = ":")]]
  )

  # 2x continuous ----
  if (setequal(x_types, "continuous")) {
    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[["x"]], y = .data[["y"]],
      color = .data[["m"]],
      fill = ggplot2::after_scale(scales::alpha(.data[["colour"]], 0.9))
      ))
    p <- p + ggplot2::geom_point(size = 2, shape = 21, stroke = 2)
    p <- p + ggplot2::scale_color_distiller(
      palette = "PRGn", type = "div",
      limits = c(-1,1) * max(abs(xdf$m)),
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(15, "char"),
        title.position = "bottom", title.hjust = .5
      )
    )
    p <- p +
      ggplot2::labs(
        x = predictors[[1]], y = predictors[[2]],
        color = sprintf("m(%s, %s)", predictors[[1]], predictors[[2]])
      )
  }

  # 2x categorical ----
  if (setequal(x_types, "categorical")) {
    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[["x"]], y = .data[["y"]],
      color = .data[["m"]],
      fill = ggplot2::after_scale(scales::alpha(.data[["colour"]], 0.9))
      ))
    p <- p + ggplot2::geom_tile()
    p <- p + ggplot2::scale_color_distiller(
      palette = "RdBu",
      limits = c(-1,1) * max(abs(xdf$m)),
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(15, "char"),
        title.position = "bottom", title.hjust = .5
      )
    )
    p <- p +
      ggplot2::labs(
        x = predictors[[1]], y = predictors[[2]],
        color = sprintf("m(%s, %s)", predictors[[1]], predictors[[2]])
      )
  }

  # 1x categorical 1x continuous ----
  if (setequal(x_types, c("categorical", "continuous"))) {

    data.table::setnames(xdf, names(xdf), c(x_types, "m"))

    p <- ggplot2::ggplot(xdf, aes(x = .data[["continuous"]], y = .data[["m"]], color = .data[["categorical"]])) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::labs(
        x = predictors[match("continuous", x_types)],
        y = sprintf("m(%s, %s)", predictors[[1]], predictors[[2]]),
        color = predictors[match("categorical", x_types)]
      )

  }

  p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

#' @rdname plot_components
#' @export
#' @examples
#' # plot_threeway_effects(components, c("hr", "temp", "workingday"))
plot_threeway_effects <- function(components, predictors, ...) {
  checkmate::assert_class(components, "rpf_components")
  checkmate::assert_character(predictors, len = 3, unique = TRUE)
  checkmate::assert_subset(predictors, names(components$x))

  # Create look-up table for predictors and their types
  tp <- c(numeric = "continuous", integer = "continuous", character = "categorical", factor = "categorical")
  x_types <- vapply(predictors, function(p) {
    cl <- class(components[["x"]][[p]])[1]
    tp[cl]
  }, "")
  checkmate::assert_subset(x_types, c("categorical", "continuous"), empty.ok = FALSE)

  xdf <- data.table::data.table(
    x1 = components[["x"]][[predictors[[1]]]],
    x2 = components[["x"]][[predictors[[2]]]],
    x3 = components[["x"]][[predictors[[3]]]],
    m = components[["m"]][[paste0(sort(predictors), collapse = ":")]]
  )

  # 3x continuous ----
  if (all(x_types == "continuous")) {
    stop("Can't visualize 3 continuous predictor effects (yet?), feel free to make a suggestion!")
  }

  # 1x categorical 2x continuous ----
  if (all(sort(x_types) == c("categorical", "continuous", "continuous"))) {
    x_cat <- names(xdf)[which(x_types == "categorical")]
    x_cont <- names(xdf)[which(x_types == "continuous")]

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cont[[1]]]],
      y = .data[[x_cont[[2]]]],
      colour = .data[["m"]],
      fill = ggplot2::after_scale(scales::alpha(.data[["colour"]], 0.9))
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[[x_cat]])) +
      ggplot2::geom_point(size = 2, shape = 21, stroke = 2) +
      ggplot2::scale_color_distiller(
        palette = "PRGn", type = "div",
        limits = c(-1,1) * max(abs(xdf$m)),
        guide = ggplot2::guide_colorbar(
          barwidth = ggplot2::unit(15, "char"),
          title.position = "bottom", title.hjust = .5
        )
      ) +
      ggplot2::labs(
        x = predictors[[1]], y = predictors[[2]],
        color = sprintf("m(%s)", paste0(predictors, collapse = ", "))
      )
  }

  # 2x categorical 1x continuous ----
  if (all(sort(x_types) == c("categorical", "categorical", "continuous"))) {
    x_cat <- names(xdf)[which(x_types == "categorical")]
    x_cont <- names(xdf)[which(x_types == "continuous")]

    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[[x_cont]],
      y = .data[["m"]],
      colour = .data[[x_cat[[1]]]],
      fill = ggplot2::after_scale(scales::alpha(.data[["colour"]], 0.9))
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[[x_cat[[2]]]]), scales = "free_y") +
      ggplot2::geom_line() +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::labs(
        x = predictors[which(x_types == "continuous")],
        y = sprintf("m(%s)", paste0(predictors, collapse = ", ")),
        color = predictors[which(x_types == "categorical")][[2]]
      )
  }
  # 3x categorical ----
  if (all(x_types == "categorical")) {
    p <- ggplot2::ggplot(xdf, ggplot2::aes(
      x = .data[["x1"]],
      y = .data[["x2"]],
      colour = .data[["m"]],
      fill = ggplot2::after_scale(scales::alpha(.data[["colour"]], 0.9))
    )) +
      ggplot2::facet_wrap(ggplot2::vars(.data[["x3"]])) +
      ggplot2::geom_tile() +
      ggplot2::scale_color_distiller(
        palette = "PRGn", type = "div",
        limits = c(-1,1) * max(abs(xdf$m)),
        guide = ggplot2::guide_colorbar(
          barwidth = ggplot2::unit(15, "char"),
          title.position = "bottom", title.hjust = .5
        )
      ) +
      ggplot2::labs(
        x = predictors[[1]], y = predictors[[2]],
        color = sprintf("m(%s)", paste0(predictors, collapse = ", "))
      )
  }

  p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom"
    )
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
