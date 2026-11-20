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
  refs <- unique(rows$ref)
  if (!baseline %in% refs) {
    stop(
      sprintf("Baseline '%s' was not measured. Use one of %s.", baseline, paste(refs, collapse = ", ")),
      call. = FALSE
    )
  }
  baseline
}

bench_deltas <- function(rows, baseline = "dev") {
  baseline <- .bench_baseline(rows, baseline)
  rows$key <- .bench_cell_key(rows)
  base <- rows[rows$ref == baseline, ]
  other <- rows[rows$ref != baseline, ]
  b <- base[match(other$key, base$key), ]
  out <- other[c("ref", BENCH_CELL_COLS)]
  out$baseline <- rep(baseline, nrow(out))
  out$time_base <- b$time_median
  out$time_ref <- other$time_median
  out$mem_base <- b$peak_delta_mb
  out$mem_ref <- other$peak_delta_mb
  out$time_pct <- .bench_pct(other$time_median, b$time_median)
  out$mem_pct <- .bench_pct(other$peak_delta_mb, b$peak_delta_mb)
  # peaks within the allocator's reuse of setup memory or the sampling interval
  out$mem_below <- !is.na(other$peak_delta_mb) & !is.na(b$peak_delta_mb) & other$peak_delta_mb < 1 & b$peak_delta_mb < 1
  out$time_verdict <- .bench_verdict(out$time_pct, BENCH_TIME_THRESHOLD, c("slower", "faster", "inconclusive"))
  out$mem_verdict <- .bench_verdict(out$mem_pct, BENCH_MEM_THRESHOLD, c("more", "less", "same"))
  out$mem_verdict[out$mem_below] <- NA_character_
  # NA when either side has no output: unknown, not unchanged
  out$output_differs <- ifelse(is.na(other$digest) | is.na(b$digest), NA, other$digest != b$digest)
  # either side could not run the op
  out$skipped <- ifelse(is.na(other$skipped), b$skipped, other$skipped)
  rownames(out) <- NULL
  out
}

.bench_cell_label <- function(rows) {
  predicted <- ifelse(rows$op %in% c("predict", "predict_components"), paste0(" predicted=", rows$n_new), "")
  paste0(
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
    ),
    predicted
  )
}

# A rows file, with columns added since it was written filled in: n_new was
# 1000 before it was a setting, raw values were not kept.
bench_read_rows <- function(path) {
  rows <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"n_new" %in% names(rows)) {
    rows$n_new <- rep(1000L, nrow(rows))
  }
  for (col in setdiff(bench_schema(), names(rows))) {
    rows[[col]] <- rep(NA, nrow(rows))
  }
  rows[bench_schema()]
}

.bench_fmt_pct <- function(x) ifelse(is.na(x), "n/a", sprintf("%+.1f%%", x))

# "baseline → ref (change)", with the unit, so a percentage never stands alone
.bench_fmt_change <- function(base, ref, pct, unit) {
  ifelse(
    is.na(pct),
    "n/a",
    sprintf("%s %s → %s %s (%s)", format(signif(base, 3)), unit, format(signif(ref, 3)), unit, .bench_fmt_pct(pct))
  )
}

# One line per ref, op and reason, with the number of cells.
.bench_skipped_lines <- function(skipped) {
  groups <- unique(skipped[c("ref", "op", "skipped")])
  vapply(
    seq_len(nrow(groups)),
    function(i) {
      g <- groups[i, ]
      k <- sum(skipped$ref == g$ref & skipped$op == g$op & skipped$skipped == g$skipped)
      sprintf("- %s cannot run %s (%d cell%s): %s", g$ref, g$op, k, if (k == 1L) "" else "s", g$skipped)
    },
    character(1)
  )
}

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
    "",
    "## How to read this",
    "",
    "- A cell is one operation on simulated data. Every ref runs it in the same job, in alternating order, so refs are compared on the same machine.",
    "- `n` is the number of training rows. Prediction ops (`predict`, `predict_components`) predict `predicted` rows, a separate set of rows.",
    "- Time is the median of seconds per call of the operation alone. Its setup, such as fitting the forest before `predict`, is not timed.",
    "- Memory in this table is what the operation adds to the peak of the R process above the level after setup. An operation can reuse memory its setup freed, which this misses, so for ops with a setup (`predict`, `purify`, `predict_components`) the number is a lower bound and can even fall as `n` grows. `below resolution` means less than 1 MB was seen. The HTML report plots the absolute peak instead.",
    sprintf(
      "- Verdicts: time changes below %d%% and memory changes below %d%% count as noise.",
      BENCH_TIME_THRESHOLD,
      BENCH_MEM_THRESHOLD
    ),
    "- `OUTPUT DIFFERS`: predictions or components are not identical to the baseline's.",
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
    # cells a ref could not run are listed under Skipped
    d <- d[is.na(d$skipped), ]
    out <- c(
      out,
      "## Deltas",
      "",
      paste0("Each change is ", baseline, " → ref. Cells a ref could not run are listed under Skipped."),
      "",
      "| cell | ref | time | verdict | memory above setup | verdict | output |",
      "|---|---|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s |",
        .bench_cell_label(d),
        d$ref,
        mapply(.bench_fmt_change, d$time_base, d$time_ref, d$time_pct, MoreArgs = list(unit = "s")),
        ifelse(is.na(d$time_verdict), "", d$time_verdict),
        ifelse(
          d$mem_below,
          "below resolution",
          mapply(.bench_fmt_change, d$mem_base, d$mem_ref, d$mem_pct, MoreArgs = list(unit = "MB"))
        ),
        ifelse(is.na(d$mem_verdict), "", d$mem_verdict),
        ifelse(is.na(d$output_differs), "n/a", ifelse(d$output_differs, "OUTPUT DIFFERS", "same"))
      ),
      ""
    )
    compared <- d[!is.na(d$output_differs), ]
    diff <- compared[compared$output_differs, ]
    if (nrow(diff) && nrow(diff) == nrow(compared)) {
      out <- c(
        out,
        "## OUTPUT DIFFERS",
        "",
        sprintf("Every compared cell (%d) gives output that differs from %s's.", nrow(diff), baseline),
        "That is expected when fits changed on purpose since the baseline, and a red flag for a change that claims identical fits.",
        ""
      )
    } else if (nrow(diff)) {
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
      .bench_skipped_lines(skipped),
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
  all <- if (file.exists(path)) rbind(bench_read_rows(path)[names(new)], new) else new
  all <- all[order(all$timestamp, decreasing = TRUE), ]
  all <- all[!duplicated(paste(.bench_cell_key(all), all$ref)), ]
  utils::write.csv(all, path, row.names = FALSE)
  invisible(path)
}

# The report for a rows file: Markdown on stdout, or a self-contained HTML
# page rendered by Quarto with this R and library. Returns the HTML path.
bench_render_report <- function(rows_csv, baseline = "dev", format = "md", dir = "bench/results") {
  bench_check_choice(format, c("md", "html"), "format")
  if (!file.exists(rows_csv)) {
    stop("Rows file ", rows_csv, " does not exist.", call. = FALSE)
  }
  if (format == "md") {
    cat(bench_report(bench_read_rows(rows_csv), baseline), sep = "\n")
    return(invisible(NULL))
  }
  if (!nzchar(Sys.which("quarto"))) {
    stop("Quarto is not on the PATH. Install it, or on the BIPS cluster run `module load quarto` first.", call. = FALSE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- paste0("report-", sub("^rows-", "", sub("\\.csv$", ".html", basename(rows_csv))))
  root <- normalizePath(".")
  # render a copy, so concurrent renders and the intermediate files stay out of bench/
  work <- tempfile("render")
  dir.create(work)
  file.copy("bench/report.qmd", work)
  # Quarto's knitr runs in this R, whose R_LIBS point to rig's library
  old <- Sys.getenv("QUARTO_R", unset = NA)
  Sys.setenv(QUARTO_R = R.home("bin"))
  on.exit(if (is.na(old)) Sys.unsetenv("QUARTO_R") else Sys.setenv(QUARTO_R = old), add = TRUE)
  status <- system2(
    "quarto",
    c(
      "render",
      shQuote(file.path(work, "report.qmd")),
      "-P",
      shQuote(paste0("rows:", normalizePath(rows_csv))),
      "-P",
      shQuote(paste0("baseline:", baseline)),
      "-P",
      shQuote(paste0("root:", root))
    )
  )
  if (status != 0L) {
    stop("Quarto failed to render the report (exit status ", status, "), see its output above.", call. = FALSE)
  }
  path <- file.path(dir, out)
  file.copy(file.path(work, "report.html"), path, overwrite = TRUE)
  normalizePath(path)
}
