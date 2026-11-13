# From measured rows to deltas against a baseline ref, a Markdown report, and
# the committed history of anchor refs.

BENCH_TIME_THRESHOLD <- 10 # percent, below this a time change is noise
BENCH_MEM_THRESHOLD <- 3 # percent

.bench_cell_key <- function(rows) {
  do.call(paste, c(lapply(rows[BENCH_CELL_COLS], as.character), sep = "|"))
}

.bench_pct <- function(new, old) {
  ifelse(is.na(new) | is.na(old) | old == 0, NA_real_, 100 * (new - old) / old)
}

.bench_verdict <- function(pct, threshold, words) {
  ifelse(is.na(pct), NA_character_, ifelse(pct >= threshold, words[1], ifelse(pct <= -threshold, words[2], words[3])))
}

.bench_baseline <- function(rows, baseline) {
  if (baseline %in% rows$ref) baseline else rows$ref[1]
}

bench_deltas <- function(rows, baseline = "dev") {
  baseline <- .bench_baseline(rows, baseline)
  rows$key <- .bench_cell_key(rows)
  base <- rows[rows$ref == baseline, ]
  other <- rows[rows$ref != baseline, ]
  b <- base[match(other$key, base$key), ]
  out <- other[c("ref", BENCH_CELL_COLS)]
  out$baseline <- rep(baseline, nrow(out))
  out$time_pct <- .bench_pct(other$time_median, b$time_median)
  out$mem_pct <- .bench_pct(other$peak_delta_mb, b$peak_delta_mb)
  # peaks within the allocator's reuse of setup memory or the sampling interval
  out$mem_below <- !is.na(other$peak_delta_mb) & !is.na(b$peak_delta_mb) & other$peak_delta_mb < 1 & b$peak_delta_mb < 1
  out$time_verdict <- .bench_verdict(out$time_pct, BENCH_TIME_THRESHOLD, c("slower", "faster", "inconclusive"))
  out$mem_verdict <- .bench_verdict(out$mem_pct, BENCH_MEM_THRESHOLD, c("more", "less", "same"))
  # NA when either side has no output: unknown, not unchanged
  out$output_differs <- ifelse(is.na(other$digest) | is.na(b$digest), NA, other$digest != b$digest)
  out$skipped <- other$skipped
  rownames(out) <- NULL
  out
}

.bench_cell_label <- function(rows) {
  sprintf(
    "%s %s%s n=%s p=%s order=%s families=%s splits=%s %s threads=%s",
    rows$op,
    rows$task,
    ifelse(is.na(rows$loss), "", paste0("/", rows$loss)),
    rows$n,
    rows$p,
    rows$max_interaction,
    rows$n_families,
    rows$n_splits,
    rows$split_structure,
    rows$nthreads
  )
}

.bench_fmt_pct <- function(x) ifelse(is.na(x), "n/a", sprintf("%+.1f%%", x))

bench_report <- function(rows, baseline = "dev") {
  baseline <- .bench_baseline(rows, baseline)
  refs <- unique(rows$ref)
  metric <- unique(rows$metric)
  out <- c(
    "# randomPlantedForest benchmark",
    "",
    paste0("- Refs: ", paste(refs, collapse = ", "), " (baseline: ", baseline, ")"),
    paste0("- Tier: ", paste(unique(rows$tier), collapse = ", ")),
    paste0("- Memory metric: ", paste(metric, collapse = ", ")),
    paste0("- Host: ", paste(unique(rows$host), collapse = ", ")),
    ""
  )
  minor <- unique(stats::na.omit(sub("^(\\d+\\.\\d+).*$", "\\1", rows$r_version)))
  if (length(minor) > 1L) {
    out <- c(
      out,
      "## R versions differ",
      "",
      paste0("These rows come from R ", paste(minor, collapse = " and "), ", so their numbers are not comparable."),
      ""
    )
  }
  if (all(metric == "none")) {
    out <- c(out, "Memory was not measured (no cgroup v2 or RSS).", "")
  }
  if (length(refs) < 2L) {
    out <- c(out, "Only one distinct ref was measured, so there is nothing to compare.", "")
  } else {
    d <- bench_deltas(rows, baseline)
    out <- c(
      out,
      "## Deltas",
      "",
      "| cell | ref | time | verdict | memory | verdict | output |",
      "|---|---|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s |",
        .bench_cell_label(d),
        d$ref,
        .bench_fmt_pct(d$time_pct),
        ifelse(is.na(d$time_verdict), "", d$time_verdict),
        ifelse(d$mem_below, "below resolution", .bench_fmt_pct(d$mem_pct)),
        ifelse(is.na(d$mem_verdict), "", d$mem_verdict),
        ifelse(
          !is.na(d$skipped),
          "skipped",
          ifelse(is.na(d$output_differs), "n/a", ifelse(d$output_differs, "OUTPUT DIFFERS", "same"))
        )
      ),
      ""
    )
    diff <- d[!is.na(d$output_differs) & d$output_differs, ]
    if (nrow(diff)) {
      out <- c(out, "## OUTPUT DIFFERS", "", sprintf("- %s: %s vs %s", .bench_cell_label(diff), diff$ref, baseline), "")
    }
  }
  broken <- rows[!is.na(rows$bound_mb) & !is.na(rows$peak_delta_mb) & rows$peak_delta_mb > rows$bound_mb, ]
  if (nrow(broken)) {
    out <- c(
      out,
      "## Memory bound exceeded",
      "",
      sprintf(
        "- %s, %s: measured %.0f MB exceeds the bound of %.0f MB",
        .bench_cell_label(broken),
        broken$ref,
        broken$peak_delta_mb,
        broken$bound_mb
      ),
      ""
    )
  }
  skipped <- rows[!is.na(rows$skipped), ]
  if (nrow(skipped)) {
    out <- c(
      out,
      "## Skipped",
      "",
      sprintf("- %s, %s: %s", .bench_cell_label(skipped), skipped$ref, skipped$skipped),
      ""
    )
  }
  missing <- attr(rows, "missing")
  if (!is.null(missing) && nrow(missing)) {
    out <- c(
      out,
      "## Missing jobs",
      "",
      "These cells have no results, for any ref:",
      "",
      sprintf("- job %d, %s %s n=%s: %s", missing$job.id, missing$op, missing$task, missing$n, missing$status),
      ""
    )
  }
  failed <- rows[!is.na(rows$error), ]
  if (nrow(failed)) {
    out <- c(out, "## Errors", "", sprintf("- %s, %s: %s", .bench_cell_label(failed), failed$ref, failed$error), "")
  }
  paste(out, collapse = "\n")
}

bench_write_results <- function(rows, dir = "bench/results", baseline = "dev") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(dir, paste0("rows-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"))
  utils::write.csv(rows, csv, row.names = FALSE)
  md <- file.path(dir, "report.md")
  writeLines(bench_report(rows, baseline), md)
  c(csv, md)
}

# Anchor rows only, newest per cell and ref: the trend across releases.
bench_history_append <- function(rows, anchors, path = "bench/history.csv") {
  missing <- attr(rows, "missing")
  if (!is.null(missing) && nrow(missing)) {
    message("Results are incomplete (", nrow(missing), " job(s) missing), so the history is not updated.")
    return(invisible(path))
  }
  new <- rows[rows$ref %in% anchors, ]
  if (!nrow(new)) {
    return(invisible(path))
  }
  all <- if (file.exists(path)) rbind(utils::read.csv(path, stringsAsFactors = FALSE)[names(new)], new) else new
  all <- all[order(all$timestamp, decreasing = TRUE), ]
  all <- all[!duplicated(paste(.bench_cell_key(all), all$ref)), ]
  utils::write.csv(all, path, row.names = FALSE)
  invisible(path)
}
