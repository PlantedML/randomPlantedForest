# Plots of measured rows: one box per ref and cell, and scaling over n.
# Built on the raw values (times, peak_deltas_mb). Rows written before those
# columns existed fall back to their summaries.

# One row per measured value: seconds per call, or the absolute peak memory of
# a round with the level after setup (floor) alongside. raw = FALSE marks
# summaries standing in for missing raw values (median, min and max of time,
# the median peak). over names the column the setting leaves out, so that
# scaling draws one line per setting.
bench_long <- function(rows, what = "time", over = "n") {
  bench_check_choice(what, c("time", "memory"), "measure")
  raw_col <- if (what == "time") "times" else "peaks_mb"
  raw <- if (raw_col %in% names(rows)) rows[[raw_col]] else rep(NA_character_, nrow(rows))
  parts <- lapply(seq_len(nrow(rows)), function(i) {
    if (!is.na(raw[i]) && nzchar(raw[i])) {
      return(data.frame(value = as.numeric(strsplit(raw[i], ";", fixed = TRUE)[[1]]), raw = TRUE))
    }
    v <- if (what == "time") {
      unlist(rows[i, intersect(c("time_median", "time_min", "time_max"), names(rows))])
    } else {
      rows$peak_mb[i]
    }
    v <- unique(stats::na.omit(as.numeric(v)))
    data.frame(value = v, raw = rep(FALSE, length(v)))
  })
  keep <- rep(seq_len(nrow(rows)), vapply(parts, nrow, integer(1)))
  # the lowest level after setup over the rounds: every round's peak is above
  # its own level, so none falls below this line, as it could below the median
  floors <- if ("floors_mb" %in% names(rows)) rows$floors_mb else rep(NA_character_, nrow(rows))
  floor <- vapply(
    seq_len(nrow(rows)),
    function(i) {
      if (is.na(floors[i]) || !nzchar(floors[i])) {
        return(as.numeric(rows$floor_mb[i]))
      }
      min(as.numeric(strsplit(floors[i], ";", fixed = TRUE)[[1]]))
    },
    numeric(1)
  )
  without <- rows
  without[[over]] <- "*"
  out <- data.frame(
    ref = rows$ref[keep],
    op = rows$op[keep],
    task = rows$task[keep],
    n = rows$n[keep],
    n_new = rows$n_new[keep],
    floor = floor[keep],
    cell = .bench_cell_label(rows)[keep],
    setting = .bench_cell_label(without)[keep],
    do.call(rbind, parts),
    stringsAsFactors = FALSE
  )
  out$ref <- factor(out$ref, levels = unique(rows$ref))
  out
}

# Stands in for a plot without values, saying why it is empty.
.bench_empty_plot <- function(reason) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, size = 6, label = reason) +
    ggplot2::theme_void()
}

.bench_no_values <- function(what) {
  if (what == "memory") "Memory was not measured in this run." else "No times in these results."
}

.bench_measure_label <- function(what) {
  if (what == "time") "seconds per call" else "peak memory of the R process (MB)"
}

# Okabe-Ito without black: distinct for colour-blind readers, one per ref
.bench_ref_colours <- function(refs) {
  stats::setNames(rep_len(unname(grDevices::palette.colors(9, "Okabe-Ito")[-1]), length(refs)), refs)
}

.bench_theme <- function() {
  ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(face = "bold"))
}

# One panel per cell, one box per ref. Points on top, so a box built from
# three summary values does not pass for a distribution. For memory, a dashed
# line per ref marks the level after setup: the part of the peak the setup
# (data, and the fitted forest for prediction ops) already held.
bench_plot_cells <- function(rows, what = "time") {
  long <- bench_long(rows, what)
  if (!nrow(long)) {
    return(.bench_empty_plot(.bench_no_values(what)))
  }
  long$cell <- sub("^\\S+ ", "", long$cell) # the op is the tab or the title
  p <- ggplot2::ggplot(long, ggplot2::aes(ref, value, colour = ref, fill = ref)) +
    ggplot2::geom_boxplot(data = long[long$raw, ], outlier.shape = NA, alpha = 0.25, linewidth = 0.9, width = 0.6)
  if (what == "memory") {
    floors <- unique(long[!is.na(long$floor), c("ref", "cell", "floor")])
    p <- p +
      ggplot2::geom_errorbar(
        data = floors,
        ggplot2::aes(ref, ymin = floor, ymax = floor, linetype = "level after setup (lowest round)"),
        inherit.aes = FALSE,
        width = 0.7,
        linewidth = 0.9
      ) +
      ggplot2::scale_linetype_manual(values = c(`level after setup (lowest round)` = "dashed"))
  }
  p +
    ggplot2::geom_point(
      ggplot2::aes(shape = raw),
      position = ggplot2::position_jitter(width = 0.12, height = 0),
      size = 3.5,
      stroke = 1.6
    ) +
    ggplot2::facet_wrap(~cell, scales = "free_y", ncol = 3, labeller = ggplot2::label_wrap_gen(35)) +
    ggplot2::scale_colour_manual(values = .bench_ref_colours(levels(long$ref)), aesthetics = c("colour", "fill")) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 4),
      labels = c(`TRUE` = "measured", `FALSE` = "summary only")
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = .bench_measure_label(what), shape = NULL, colour = NULL, fill = NULL, linetype = NULL) +
    .bench_theme()
}

# Median per ref over training rows (n) or predicted rows (n_new), log axes,
# one panel per op and task, one line per remaining setting. Settings measured
# at a single size are left out.
bench_plot_scaling <- function(rows, what = "time", over = "n") {
  bench_check_choice(over, c("n", "n_new"), "scaling axis")
  long <- bench_long(rows, what, over = over)
  if (!nrow(long)) {
    return(.bench_empty_plot(.bench_no_values(what)))
  }
  long$x <- long[[over]]
  med <- stats::aggregate(value ~ ref + op + task + x + setting, data = long, FUN = stats::median)
  sizes <- tapply(med$x, med$setting, function(x) length(unique(x)))
  med <- med[med$setting %in% names(sizes)[sizes > 1L], ]
  if (!nrow(med)) {
    return(.bench_empty_plot(sprintf(
      "No cell in these results was measured at more than one number of %s.",
      if (over == "n") "training rows" else "predicted rows"
    )))
  }
  med[[over]] <- med$x
  ggplot2::ggplot(med, ggplot2::aes(x, value, colour = ref, group = interaction(ref, setting))) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(values = .bench_ref_colours(levels(long$ref))) +
    ggplot2::scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ op + task, scales = "free_y") +
    ggplot2::labs(
      x = if (over == "n") "training rows (n)" else "predicted rows",
      y = .bench_measure_label(what),
      colour = NULL
    ) +
    .bench_theme()
}
