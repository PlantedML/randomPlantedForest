# Plots of measured rows: one box per ref and cell, and scaling over n.
# Built on the raw values (times, peak_deltas_mb). Rows written before those
# columns existed fall back to their summaries.

# One row per measured value. raw = FALSE marks summaries standing in for
# missing raw values (median, min and max of time, the median of memory).
bench_long <- function(rows, what = "time") {
  bench_check_choice(what, c("time", "memory"), "measure")
  raw_col <- if (what == "time") "times" else "peak_deltas_mb"
  raw <- if (raw_col %in% names(rows)) rows[[raw_col]] else rep(NA_character_, nrow(rows))
  parts <- lapply(seq_len(nrow(rows)), function(i) {
    if (!is.na(raw[i]) && nzchar(raw[i])) {
      return(data.frame(value = as.numeric(strsplit(raw[i], ";", fixed = TRUE)[[1]]), raw = TRUE))
    }
    v <- if (what == "time") {
      unlist(rows[i, intersect(c("time_median", "time_min", "time_max"), names(rows))])
    } else {
      rows$peak_delta_mb[i]
    }
    v <- unique(stats::na.omit(as.numeric(v)))
    data.frame(value = v, raw = rep(FALSE, length(v)))
  })
  n <- vapply(parts, nrow, integer(1))
  keep <- rep(seq_len(nrow(rows)), n)
  out <- data.frame(
    ref = rows$ref[keep],
    op = rows$op[keep],
    task = rows$task[keep],
    n = rows$n[keep],
    cell = .bench_cell_label(rows)[keep],
    # the cell without n, to draw one line per setting over n
    setting = .bench_cell_label(transform(rows, n = "*"))[keep],
    do.call(rbind, parts),
    stringsAsFactors = FALSE
  )
  out$ref <- factor(out$ref, levels = unique(rows$ref))
  out
}

.bench_measure_label <- function(what) {
  if (what == "time") "seconds per call" else "peak memory above setup (MB)"
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
# three summary values does not pass for a distribution.
bench_plot_cells <- function(rows, what = "time") {
  long <- bench_long(rows, what)
  long$cell <- sub("^\\S+ ", "", long$cell) # the op is the tab or the title
  ggplot2::ggplot(long, ggplot2::aes(ref, value, colour = ref, fill = ref)) +
    ggplot2::geom_boxplot(data = long[long$raw, ], outlier.shape = NA, alpha = 0.25, linewidth = 0.9, width = 0.6) +
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
    ggplot2::labs(x = NULL, y = .bench_measure_label(what), shape = NULL, colour = NULL, fill = NULL) +
    .bench_theme()
}

# Median per ref over n, log axes, one panel per op and task, one line per
# remaining setting. Settings measured at a single size are left out.
bench_plot_scaling <- function(rows, what = "time") {
  long <- bench_long(rows, what)
  med <- stats::aggregate(value ~ ref + op + task + n + setting, data = long, FUN = stats::median)
  sizes <- tapply(med$n, med$setting, function(n) length(unique(n)))
  med <- med[med$setting %in% names(sizes)[sizes > 1L], ]
  ggplot2::ggplot(med, ggplot2::aes(n, value, colour = ref, group = interaction(ref, setting))) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(values = .bench_ref_colours(levels(long$ref))) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ op + task, scales = "free_y") +
    ggplot2::labs(x = "rows (n)", y = .bench_measure_label(what), colour = NULL) +
    .bench_theme()
}
