#' Extract the boxes of a tree family
#'
#' Every tree in a Random Planted Forest is a set of boxes (leaves) with a
#' value each, and the prediction is the sum of the values of all boxes
#' containing the point, averaged over tree families.
#' `rpf_boxes()` returns the boxes of one tree family as a data frame.
#'
#' Before purification, these are the leaves as grown: boxes in one tree may
#' overlap, and the trees are not an identified decomposition.
#' After [purify()], they are the cells of the purified grid, one
#' non-overlapping set per component of the functional decomposition.
#' As [purify()] modifies the forest in place, extract the unpurified boxes
#' before purifying.
#'
#' Factor predictors are encoded as integer codes in the order of
#' `x$factor_levels`, with a box covering the codes `c` with
#' `lower <= c < upper`.
#'
#' @param x `[rpf]`: A fitted [`rpf`] model.
#' @param family `[integer(1): 1]`: Which tree family to extract, between 1 and `ntrees`.
#' @param ... Reserved for future expansion.
#'
#' @return A [base::data.frame] with one row per box and columns
#'   * `tree`: The tree, named by its sorted variables as in [predict_components()],
#'     e.g. `"x1:x2"`, or `"(Intercept)"`.
#'   * `order`: The number of variables of the tree, `0` for the intercept.
#'   * `box`: The index of the box within its tree.
#'   * `value`: The value of the box. For multiclass classification, one
#'     column `value_<class>` per outcome level instead.
#'   * `<variable>_lower`, `<variable>_upper` for every predictor: The bounds of
#'     the box, `NA` for variables the tree does not use.
#'
#' @seealso [plot.rpf()] to plot the boxes.
#' @export
#' @examples
#' rpfit <- rpf(mpg ~ wt + hp, data = mtcars, ntrees = 5, splits = 10)
#' head(rpf_boxes(rpfit))
rpf_boxes <- function(x, family = 1L, ...) {
  checkmate::assert_class(x, "rpf")
  rlang::check_dots_empty()
  check_rpf_alive(x)
  checkmate::assert_int(family, lower = 1, upper = x$params$ntrees)

  # mold() moves factors last, so the model's column order can differ from the ptype's
  predictors <- names(hardhat::forge(x$blueprint$ptypes$predictors, x$blueprint)$predictors)
  boxes <- if (is_purified(x)) {
    grid_boxes(x$fit$get_grid_leaves()[[family]])
  } else {
    leaf_boxes(x$fit$get_model()[[family]])
  }

  tree <- vapply(
    boxes$variables,
    \(vars) if (all(vars == 0)) "(Intercept)" else paste(sort(predictors[vars]), collapse = ":"),
    character(1)
  )
  out <- data.frame(
    tree = tree,
    order = vapply(boxes$variables, \(vars) sum(vars > 0), integer(1)),
    box = boxes$box
  )

  values <- boxes$values
  if (ncol(values) == 1) {
    out$value <- values[, 1]
  } else {
    colnames(values) <- paste0("value_", levels(x$blueprint$ptypes$outcomes[[1]]))
    out <- cbind(out, values)
  }

  for (j in seq_along(predictors)) {
    used <- vapply(boxes$variables, \(vars) j %in% vars, logical(1))
    out[[paste0(predictors[j], "_lower")]] <- ifelse(used, boxes$lower[, j], NA_real_)
    out[[paste0(predictors[j], "_upper")]] <- ifelse(used, boxes$upper[, j], NA_real_)
  }
  out
}

# Raw leaves from get_model(): one 2 x p bounds matrix per leaf
leaf_boxes <- function(family) {
  n_leaves <- lengths(family$values)
  leaves <- unlist(family$intervals, recursive = FALSE)
  list(
    variables = rep(family$variables, n_leaves),
    box = unlist(lapply(n_leaves, seq_len)),
    values = do.call(rbind, unlist(family$values, recursive = FALSE)),
    lower = do.call(rbind, lapply(leaves, \(bounds) bounds[1, ])),
    upper = do.call(rbind, lapply(leaves, \(bounds) bounds[2, ]))
  )
}

# Purified grid from get_grid_leaves(): cells between consecutive limits of
# each tree variable, values in column-major order over the tree's dims.
# The last index per dimension is padding and never used for prediction.
grid_boxes <- function(family) {
  p <- length(family$lim_list)
  trees <- lapply(family$trees, \(tree) {
    vars <- tree$variables
    if (all(vars == 0)) {
      cells <- matrix(1L, 1, 0)
      values <- tree$values[1, , drop = FALSE]
    } else {
      used <- lapply(family$lim_list[vars], \(lim) seq_len(length(lim) - 1))
      cells <- as.matrix(expand.grid(used))
      offsets <- cumprod(c(1, tree$dims[-length(tree$dims)]))
      values <- tree$values[drop(1 + (cells - 1) %*% offsets), , drop = FALSE]
    }
    lower <- upper <- matrix(NA_real_, nrow(values), p)
    for (k in seq_along(vars[vars > 0])) {
      lim <- family$lim_list[[vars[k]]]
      lower[, vars[k]] <- lim[cells[, k]]
      upper[, vars[k]] <- lim[cells[, k] + 1]
    }
    list(lower = lower, upper = upper, values = values)
  })
  n_boxes <- vapply(trees, \(tree) nrow(tree$values), integer(1))
  list(
    variables = rep(lapply(family$trees, `[[`, "variables"), n_boxes),
    box = unlist(lapply(n_boxes, seq_len)),
    values = do.call(rbind, lapply(trees, `[[`, "values")),
    lower = do.call(rbind, lapply(trees, `[[`, "lower")),
    upper = do.call(rbind, lapply(trees, `[[`, "upper"))
  )
}

#' Plot a tree family
#'
#' Illustrates how a Random Planted Forest predicts by drawing the trees of one
#' tree family. Each tree belongs to a set of variables, and its leaves are
#' boxes in those variables with a value each. The prediction for a point is
#' the root value plus the values of all leaves containing it, averaged over
#' tree families.
#'
#' * `type = "tree"` draws a split diagram in the style of Figure 1b of
#'   Hiabu et al. (2023), with the trees in rows by number of variables, each
#'   leaf filled by its value. The forest stores only its final leaves, not
#'   the order of splits or which leaf a tree grew from, so the splits within
#'   a tree are reconstructed: the diagram shows a split structure consistent
#'   with the leaves, and dashed lines link trees to the trees they could
#'   have grown from. Leaves that overlap other leaves of the same tree hang
#'   from dotted lines. Readable for small fits, e.g. `ntrees = 1` and
#'   `splits = 10`.
#' * `type = "boxes"` draws one panel per tree in variable space: leaves of
#'   trees on one variable as bars from zero to their value, leaves of trees
#'   on two variables as rectangles filled by their value. Overlapping leaves
#'   add up, which the translucent drawing shows.
#'
#' Trees on more than `max_interaction` variables, and trees not in `trees`,
#' are not drawn but summarized as a remainder, similar to the remainder term
#' of [predict_components()].
#'
#' After [purify()], the trees are the components of the functional
#' decomposition and only `type = "boxes"` applies. To plot components
#' averaged over all tree families, use [predict_components()], e.g. with the
#' plotting functions in `glex`.
#'
#' @inheritParams rpf_boxes
#' @param type `[character(1): "tree"]`: `"tree"` for a split diagram,
#'   `"boxes"` for the leaves in variable space.
#' @param max_interaction `[integer(1): 2]`: Draw trees on up to this many
#'   variables, at most 2 for `type = "boxes"`. The remaining trees are
#'   summarized.
#' @param trees `[character | NULL: NULL]`: Names of the trees to draw, as in
#'   the `tree` column of [rpf_boxes()], e.g. `c("x1", "x1:x2")`. `NULL` draws
#'   all trees up to `max_interaction`.
#' @param class `[character(1) | NULL: NULL]`: For multiclass classification,
#'   the outcome level whose values to plot. `NULL` uses the first level.
#' @param factor_levels `[logical(1): TRUE]`: Label factors with their
#'   levels instead of integer codes. Levels are ordered by their association
#'   with the outcome, as used in the fit.
#'
#' @return A `ggplot` for `type = "tree"`, a `patchwork` for `type = "boxes"`.
#' @references Hiabu, M., Mammen, E., Meyer, J. T. (2023). Random Planted
#'   Forest: a directly interpretable tree ensemble. arXiv:2012.14563.
#' @importFrom rlang .data
#' @export
#' @examplesIf rlang::is_installed(c("ggplot2", "patchwork"))
#' rpfit <- rpf(mpg ~ wt + hp + cyl, data = mtcars, ntrees = 1, splits = 10)
#' plot(rpfit)
#' plot(rpfit, type = "boxes")
#' purify(rpfit)
#' plot(rpfit, type = "boxes")
plot.rpf <- function(
  x,
  family = 1L,
  ...,
  type = c("tree", "boxes"),
  max_interaction = 2L,
  trees = NULL,
  class = NULL,
  factor_levels = TRUE
) {
  rlang::check_dots_empty()
  type <- rlang::arg_match(type)
  rlang::check_installed(
    c("ggplot2", if (type == "boxes") "patchwork"),
    reason = "to plot a tree family."
  )
  checkmate::assert_int(max_interaction, lower = 1, upper = if (type == "boxes") 2 else Inf)
  checkmate::assert_flag(factor_levels)
  if (type == "tree" && is_purified(x)) {
    cli::cli_abort(c(
      "The forest is purified, so its trees no longer have the leaves they were grown with.",
      "i" = "Use {.code type = \"boxes\"}, or plot before {.fn purify}."
    ))
  }
  boxes <- rpf_boxes(x, family = family)
  column <- value_column(x, boxes, class)
  boxes$.value <- boxes[[column]]
  levels <- if (factor_levels) x$factor_levels else list()
  value_label <- if (column == "value") "leaf value" else sprintf("leaf value (%s)", sub("^value_", "", column))
  if (!is.null(trees)) {
    checkmate::assert_subset(trees, unique(boxes$tree[boxes$order > 0]))
  }
  drawn <- boxes$order == 0 | (boxes$order <= max_interaction & (is.null(trees) | boxes$tree %in% trees))
  remainder <- summarize_remainder(boxes[!drawn, , drop = FALSE])
  boxes <- boxes[drawn, , drop = FALSE]

  switch(
    type,
    tree = plot_family_tree(boxes, family, levels, value_label, remainder),
    boxes = plot_family_boxes(boxes, family, levels, value_label, remainder)
  )
}

# Trees not drawn, summarized like the remainder term of predict_components();
# NULL if there are none
summarize_remainder <- function(high) {
  if (nrow(high) == 0) {
    return(NULL)
  }
  orders <- range(high$order)
  variables <- if (orders[1] == orders[2]) orders[1] else paste(orders, collapse = "-")
  sprintf(
    "Remainder: %d tree%s on %s variables with %d leaves, not drawn",
    length(unique(high$tree)),
    if (length(unique(high$tree)) == 1) "" else "s",
    variables,
    nrow(high)
  )
}

plot_family_boxes <- function(boxes, family, levels, value_label, remainder) {
  mains <- unique(boxes$tree[boxes$order == 1])
  main_limits <- range(0, boxes$.value[boxes$order == 1])
  pairs <- unique(boxes$tree[boxes$order == 2])
  pair_limits <- c(-1, 1) * max(abs(boxes$.value[boxes$order == 2]), 0)

  panels <- c(
    lapply(mains, \(tree) plot_main_tree(boxes[boxes$tree == tree, ], tree, main_limits, levels, value_label)),
    lapply(pairs, \(tree) plot_pair_tree(boxes[boxes$tree == tree, ], tree, pair_limits, levels, value_label))
  )
  intercept <- boxes$.value[boxes$order == 0]
  patchwork::wrap_plots(panels, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Tree family %d: leaves of each tree", family),
      subtitle = sprintf(
        "Prediction = root value (%s) + values of all leaves containing the point, summed over panels",
        format(sum(intercept), digits = 3)
      ),
      caption = remainder
    )
}

value_column <- function(x, boxes, class) {
  if (!"value" %in% names(boxes)) {
    outcome_levels <- levels(x$blueprint$ptypes$outcomes[[1]])
    class <- class %||% outcome_levels[1]
    checkmate::assert_choice(class, outcome_levels)
    return(paste0("value_", class))
  }
  if (!is.null(class)) {
    cli::cli_abort("{.arg class} only applies to multiclass classification.")
  }
  "value"
}

# Bounds of a box along one variable, as plot coordinates. A factor box
# [lower, upper) covers the codes from ceiling(lower) to ceiling(upper) - 1,
# drawn as unit-wide cells centred on the codes.
axis_bounds <- function(boxes, variable, levels) {
  lower <- boxes[[paste0(variable, "_lower")]]
  upper <- boxes[[paste0(variable, "_upper")]]
  if (variable %in% names(levels)) {
    list(min = ceiling(lower) - 0.5, max = ceiling(upper) - 0.5)
  } else {
    list(min = lower, max = upper)
  }
}

axis_scale <- function(aesthetic, variable, levels) {
  scale <- if (aesthetic == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
  if (variable %in% names(levels)) {
    scale(variable, breaks = seq_along(levels[[variable]]), labels = levels[[variable]])
  } else {
    scale(variable)
  }
}

plot_main_tree <- function(boxes, variable, limits, levels, value_label) {
  xs <- axis_bounds(boxes, variable, levels)
  rects <- data.frame(xmin = xs$min, xmax = xs$max, ymin = 0, ymax = boxes$.value)
  ggplot2::ggplot(rects) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      alpha = 0.5,
      fill = "grey30"
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    axis_scale("x", variable, levels) +
    ggplot2::scale_y_continuous(value_label, limits = limits) +
    ggplot2::labs(title = variable) +
    ggplot2::theme_minimal()
}

plot_pair_tree <- function(boxes, tree, limits, levels, value_label) {
  bounds <- grep("_lower$", names(boxes), value = TRUE)
  variables <- sort(sub("_lower$", "", bounds[!is.na(unlist(boxes[1, bounds]))]))
  xs <- axis_bounds(boxes, variables[1], levels)
  ys <- axis_bounds(boxes, variables[2], levels)
  rects <- data.frame(xmin = xs$min, xmax = xs$max, ymin = ys$min, ymax = ys$max, value = boxes$.value)
  ggplot2::ggplot(rects) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$value
      ),
      alpha = 0.7
    ) +
    axis_scale("x", variables[1], levels) +
    axis_scale("y", variables[2], levels) +
    ggplot2::scale_fill_gradient2(value_label, limits = limits) +
    ggplot2::labs(title = tree) +
    ggplot2::theme_minimal()
}
