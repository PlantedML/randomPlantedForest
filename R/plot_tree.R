# Split diagram of a tree family, in the style of Figure 1b of the rpf paper.
# The fit stores only final leaves, so within each tree a split structure is
# rebuilt from the leaves: any structure consistent with them, not the history.

# Split structure ------------------------------------------------------------

# A tree's leaves are one or more partitions of a region: a leaf split on a
# new variable stays and can be split again, adding another partition.
# Repeatedly take a subset of leaves that exactly tiles the region as a split
# tree; leaves left over hang from the region individually.
split_leaves <- function(leaves, variables, region) {
  groups <- list()
  remaining <- leaves
  budget <- new.env()
  budget$calls <- 0
  while (nrow(remaining) > 0) {
    tiling <- tile_region(remaining, variables, region, budget)
    if (is.null(tiling)) {
      break
    }
    groups[[length(groups) + 1]] <- tiling
    remaining <- remaining[!remaining$.id %in% tiling_ids(tiling), , drop = FALSE]
  }
  loose <- lapply(seq_len(nrow(remaining)), \(i) {
    list(kind = "leaf", leaf = remaining[i, ], region = leaf_region(remaining[i, ], variables))
  })
  children <- c(groups, loose)
  if (length(children) == 1) {
    return(children[[1]])
  }
  list(kind = "overlap", children = children, region = region)
}

# A split tree whose leaves exactly tile `region`, or NULL. Prefers a single
# leaf, then the most balanced cut. `budget` caps the search on large trees.
tile_region <- function(leaves, variables, region, budget) {
  budget$calls <- budget$calls + 1
  if (budget$calls > 5000) {
    return(NULL)
  }
  inside <- rep(TRUE, nrow(leaves))
  exact <- rep(TRUE, nrow(leaves))
  for (variable in variables) {
    lower <- leaves[[paste0(variable, "_lower")]]
    upper <- leaves[[paste0(variable, "_upper")]]
    inside <- inside & lower >= region[[variable]][1] & upper <= region[[variable]][2]
    exact <- exact & lower == region[[variable]][1] & upper == region[[variable]][2]
  }
  if (any(exact)) {
    leaf <- leaves[which(exact)[1], ]
    return(list(kind = "leaf", leaf = leaf, region = region))
  }
  leaves <- leaves[inside, , drop = FALSE]
  if (nrow(leaves) < 2) {
    return(NULL)
  }
  for (cut in candidate_cuts(leaves, variables, region)) {
    region_low <- region_high <- region
    region_low[[cut$variable]][2] <- cut$value
    region_high[[cut$variable]][1] <- cut$value
    low <- tile_region(leaves, variables, region_low, budget)
    if (is.null(low)) {
      next
    }
    high <- tile_region(leaves, variables, region_high, budget)
    if (!is.null(high)) {
      return(list(
        kind = "split",
        variable = cut$variable,
        value = cut$value,
        region = region,
        children = list(low, high)
      ))
    }
  }
  NULL
}

# Leaf bounds strictly inside the region, most balanced first
candidate_cuts <- function(leaves, variables, region) {
  cuts <- list()
  for (variable in variables) {
    lower <- leaves[[paste0(variable, "_lower")]]
    upper <- leaves[[paste0(variable, "_upper")]]
    values <- unique(c(lower, upper))
    values <- values[values > region[[variable]][1] & values < region[[variable]][2]]
    for (value in values) {
      balance <- abs(sum(upper <= value) - sum(lower >= value))
      cuts[[length(cuts) + 1]] <- list(variable = variable, value = value, balance = balance)
    }
  }
  cuts[order(vapply(cuts, `[[`, numeric(1), "balance"))]
}

tiling_ids <- function(node) {
  if (node$kind == "leaf") {
    return(node$leaf$.id)
  }
  unlist(lapply(node$children, tiling_ids))
}

leaf_region <- function(leaf, variables) {
  stats::setNames(
    lapply(variables, \(v) c(leaf[[paste0(v, "_lower")]], leaf[[paste0(v, "_upper")]])),
    variables
  )
}

# Labels ---------------------------------------------------------------------

# Describe where `region` is narrower than `parent`, e.g. "x1 < 0.5" or "g in {a, c}"
region_label <- function(region, parent, levels) {
  parts <- character()
  for (variable in names(region)) {
    bounds <- region[[variable]]
    outer <- parent[[variable]]
    if (variable %in% names(levels)) {
      codes <- seq_along(levels[[variable]])
      covered <- codes[codes >= bounds[1] & codes < bounds[2]]
      if (length(covered) < sum(codes >= outer[1] & codes < outer[2])) {
        parts <- c(parts, sprintf("%s \u2208 {%s}", variable, paste(levels[[variable]][covered], collapse = ", ")))
      }
      next
    }
    low <- bounds[1] > outer[1]
    high <- bounds[2] < outer[2]
    lower <- format(bounds[1], digits = 3)
    upper <- format(bounds[2], digits = 3)
    if (low && high) {
      parts <- c(parts, sprintf("%s \u2264 %s < %s", lower, variable, upper))
    } else if (low) {
      parts <- c(parts, sprintf("%s \u2265 %s", variable, lower))
    } else if (high) {
      parts <- c(parts, sprintf("%s < %s", variable, upper))
    }
  }
  paste(parts, collapse = "\n")
}

# Layout ---------------------------------------------------------------------

box_half_width <- 0.42
box_half_height <- 0.22

# Place leaves on consecutive slots and inner nodes above the mean of their
# children; returns nodes and elbow edges in tree-local coordinates.
layout_split_tree <- function(node, levels) {
  nodes <- list()
  edges <- list()
  next_slot <- 0

  place <- function(node, depth, parent_region) {
    label <- region_label(node$region, parent_region, levels)
    if (node$kind == "leaf") {
      x <- next_slot
      next_slot <<- next_slot + 1
      nodes[[length(nodes) + 1]] <<- data.frame(
        x = x,
        y = -depth,
        kind = "leaf",
        value = node$leaf$.value,
        label = label
      )
      return(x)
    }
    xs <- vapply(node$children, \(child) place(child, depth + 1, node$region), numeric(1))
    x <- mean(xs)
    # inner nodes narrowing their parent's region are boxes; the top split and
    # repeated splits of the same region hang directly from the line above
    boxed <- depth > 0 && nzchar(label)
    if (boxed) {
      nodes[[length(nodes) + 1]] <<- data.frame(x = x, y = -depth, kind = node$kind, value = NA_real_, label = label)
    }
    # elbow: down from the node to a bar along the tops of the children
    top <- if (boxed) {
      -depth - box_half_height
    } else if (depth > 0) {
      -depth + box_half_height
    } else {
      -depth
    }
    bar <- -depth - 1 + box_half_height
    edges[[length(edges) + 1]] <<- data.frame(
      x = c(x, min(xs)),
      xend = c(x, max(xs)),
      y = c(top, bar),
      yend = c(bar, bar),
      overlap = node$kind == "overlap"
    )
    x
  }

  top_x <- place(node, 0, node$region)
  list(nodes = do.call(rbind, nodes), edges = do.call(rbind, edges), width = next_slot, top_x = top_x)
}

shift_layout <- function(df, dx, dy) {
  df$x <- df$x + dx
  df$y <- df$y + dy
  if (!is.null(df$xend)) {
    df$xend <- df$xend + dx
    df$yend <- df$yend + dy
  }
  df
}

# Arrange the trees in rows by interaction order below the root, each tree
# with a header naming its variables and, where narrower than the data, its region.
layout_family <- function(boxes, levels) {
  predictors <- sub("_lower$", "", grep("_lower$", names(boxes), value = TRUE))
  predictors <- predictors[vapply(predictors, \(v) any(!is.na(boxes[[paste0(v, "_lower")]])), logical(1))]
  data_range <- stats::setNames(
    lapply(predictors, \(v) range(boxes[[paste0(v, "_lower")]], boxes[[paste0(v, "_upper")]], na.rm = TRUE)),
    predictors
  )

  trees <- unique(boxes$tree[boxes$order > 0])
  tree_order <- boxes$order[match(trees, boxes$tree)]
  trees <- trees[order(tree_order)]
  tree_order <- sort(tree_order)

  layouts <- lapply(trees, \(tree) {
    leaves <- boxes[boxes$tree == tree, , drop = FALSE]
    leaves$.id <- seq_len(nrow(leaves))
    variables <- strsplit(tree, ":", fixed = TRUE)[[1]]
    region <- stats::setNames(
      lapply(variables, \(v) range(leaves[[paste0(v, "_lower")]], leaves[[paste0(v, "_upper")]])),
      variables
    )
    l <- layout_split_tree(split_leaves(leaves, variables, region), levels)
    l$region_label <- region_label(region, data_range[variables], levels)
    l
  })
  depths <- vapply(layouts, \(l) -min(l$nodes$y), numeric(1))

  gap <- 1
  row_top <- 0
  nodes <- edges <- headers <- list()
  for (o in unique(tree_order)) {
    in_row <- which(tree_order == o)
    widths <- vapply(layouts[in_row], `[[`, numeric(1), "width")
    row_top <- row_top - 1.3
    # centre the row on the root at x = 0
    x_offset <- -(sum(widths) + gap * (length(in_row) - 1) - 1) / 2
    for (i in in_row) {
      l <- layouts[[i]]
      nodes[[i]] <- shift_layout(l$nodes, x_offset, row_top)
      if (!is.null(l$edges)) {
        edges[[i]] <- shift_layout(l$edges, x_offset, row_top)
      }
      headers[[i]] <- data.frame(
        tree = trees[i],
        region = l$region_label,
        order = o,
        x = x_offset + l$top_x,
        y = row_top
      )
      x_offset <- x_offset + l$width + gap
    }
    row_top <- row_top - max(depths[in_row])
  }
  nodes <- do.call(rbind, nodes)
  edges <- do.call(rbind, edges)
  headers <- do.call(rbind, headers)
  list(nodes = nodes, edges = edges, headers = headers, trees = trees)
}

# Plot -----------------------------------------------------------------------

plot_family_tree <- function(boxes, family, levels, value_label, remainder) {
  if (!any(boxes$order > 0)) {
    cli::cli_abort("Tree family {family} has no trees to draw.")
  }
  n_leaves <- sum(boxes$order > 0)
  if (n_leaves > 60) {
    cli::cli_inform(c(
      "Drawing {n_leaves} leaves, which may be hard to read.",
      "i" = "Fit with fewer {.arg splits} for an illustration, or use {.code type = \"boxes\"}."
    ))
  }
  layout <- layout_family(boxes, levels)
  # widen slots so that split labels of long variable names fit
  label_chars <- max(nchar(unlist(strsplit(layout$nodes$label, "\n", fixed = TRUE))), 0)
  slot <- max(1, label_chars / 12)
  layout$nodes$x <- layout$nodes$x * slot
  layout$edges[c("x", "xend")] <- layout$edges[c("x", "xend")] * slot
  layout$headers$x <- layout$headers$x * slot
  nodes <- layout$nodes
  headers <- layout$headers
  headers$label <- sprintf("{%s}", gsub(":", ", ", headers$tree, fixed = TRUE))
  headers$label <- ifelse(
    nzchar(headers$region),
    paste0(headers$label, "\n", gsub("\n", ", ", headers$region, fixed = TRUE)),
    headers$label
  )
  root <- data.frame(x = 0, y = 0)

  # which leaf a higher-order tree grew from is not stored, so only trees on
  # one variable, which grow from the root, are linked
  links <- list(data.frame(
    x = root$x,
    y = root$y - 0.3,
    xend = headers$x[headers$order == 1],
    yend = headers$y[headers$order == 1] + 0.35
  ))
  links <- do.call(rbind, links)
  intercept <- boxes$.value[boxes$order == 0]
  row_labels <- unique(headers[c("order", "y")])
  row_labels$label <- sprintf("trees on\n%d variable%s", row_labels$order, ifelse(row_labels$order == 1, "", "s"))

  half_w <- box_half_width * slot
  half_h <- box_half_height
  leaves <- nodes[nodes$kind == "leaf", ]
  inner <- nodes[nodes$kind != "leaf", ]
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      data = links,
      colour = "grey70",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linetype = .data$overlap),
      data = layout$edges,
      colour = "grey30",
      show.legend = FALSE
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$x - half_w, xmax = .data$x + half_w, ymin = .data$y - half_h, ymax = .data$y + half_h),
      data = inner,
      fill = "white",
      colour = "grey30"
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$x - half_w,
        xmax = .data$x + half_w,
        ymin = .data$y - half_h,
        ymax = .data$y + half_h,
        fill = .data$value
      ),
      data = leaves,
      colour = "grey30"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(.data$x, .data$y, label = format(.data$value, digits = 2)),
      data = leaves,
      size = 2.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(.data$x, .data$y + half_h, label = .data$label),
      data = nodes,
      size = 2.3,
      vjust = -0.3,
      lineheight = 0.9
    ) +
    ggplot2::geom_label(
      ggplot2::aes(.data$x, .data$y + 0.05, label = .data$label),
      data = headers,
      size = 3,
      vjust = 0,
      lineheight = 0.9,
      label.size = 0
    ) +
    ggplot2::geom_label(
      ggplot2::aes(.data$x, .data$y, label = sprintf("root\n%s", format(intercept, digits = 3))),
      data = root,
      size = 3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(min(nodes$x) - 1, .data$y, label = .data$label),
      data = row_labels,
      hjust = 1,
      colour = "grey40",
      size = 3
    ) +
    remainder_layers(remainder, min(nodes$y) - half_h - 1) +
    ggplot2::scale_fill_gradient2(value_label) +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dotted")) +
    ggplot2::labs(
      title = sprintf("Tree family %d", family),
      subtitle = "Prediction = root value + values of all leaves containing the point"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 60))
}

# A label below the diagram standing in for the trees not drawn
remainder_layers <- function(remainder, y) {
  if (is.null(remainder)) {
    return(NULL)
  }
  ggplot2::annotate(
    "label",
    x = 0,
    y = y,
    label = remainder,
    size = 3,
    colour = "grey30",
    fill = "grey95"
  )
}
