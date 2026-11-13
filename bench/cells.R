# The benchmark grid as data: one row per cell, all refs measured per cell.

BENCH_CELL_COLS <- c(
  "op",
  "task",
  "loss",
  "n",
  "p",
  "max_interaction",
  "n_families",
  "n_splits",
  "split_structure",
  "nthreads",
  "iters",
  "mem_reps"
)

.bench_grid <- function(...) {
  g <- expand.grid(..., stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  defaults <- list(
    loss = NA_character_,
    p = 10L,
    max_interaction = 2L,
    n_families = 50L,
    n_splits = 100L,
    split_structure = "leaf",
    nthreads = 1L,
    iters = 5L,
    mem_reps = 5L
  )
  for (nm in setdiff(names(defaults), names(g))) {
    g[[nm]] <- defaults[[nm]]
  }
  g[BENCH_CELL_COLS]
}

# A few minutes on a laptop for three refs.
.bench_cells_quick <- function() {
  .bench_grid(
    op = c("fit_fixed", "fit_oob", "predict", "purify"),
    task = c("regr", "binary"),
    n = c(1000L, 5000L),
    p = 5L,
    n_families = 20L,
    n_splits = 30L,
    iters = 3L,
    mem_reps = 3L
  )
}

# Focused sub-grids rather than one product, which would be thousands of cells.
.bench_cells_full <- function() {
  classif <- expand.grid(task = c("binary", "multiclass"), loss = c("logit", "exponential"), stringsAsFactors = FALSE)
  rbind(
    # scaling in n and order
    .bench_grid(
      op = c("fit_fixed", "fit_oob", "predict", "purify", "predict_components"),
      task = "regr",
      n = c(10000L, 50000L),
      max_interaction = c(2L, 3L)
    ),
    # split structures
    .bench_grid(
      op = "fit_fixed",
      task = "regr",
      n = 10000L,
      max_interaction = 3L,
      split_structure = c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")
    ),
    # classification losses
    do.call(
      rbind,
      lapply(seq_len(nrow(classif)), function(i) {
        .bench_grid(
          op = c("fit_fixed", "fit_oob", "predict"),
          task = classif$task[i],
          loss = classif$loss[i],
          n = 10000L
        )
      })
    ),
    # threads
    .bench_grid(
      op = c("fit_fixed", "predict"),
      task = "regr",
      n = 50000L,
      max_interaction = 3L,
      n_families = 200L,
      nthreads = c(1L, 4L, 16L)
    ),
    # the default call
    .bench_grid(op = "fit_default", task = c("regr", "binary"), n = 10000L, max_interaction = 3L, n_families = 200L)
  )
}

bench_cells <- function(tier = c("quick", "full")) {
  tier <- match.arg(tier)
  cells <- switch(tier, quick = .bench_cells_quick(), full = .bench_cells_full())
  # purify() works in place, see bench_ops()
  cells$iters[cells$op == "purify"] <- 1L
  cells$tier <- tier
  cells <- unique(cells)
  rownames(cells) <- NULL
  cells
}

bench_select_cells <- function(cells, ops = NA, tasks = NA, sizes = NA, max_cells = NA) {
  keep <- rep(TRUE, nrow(cells))
  if (!all(is.na(ops))) {
    keep <- keep & cells$op %in% ops
  }
  if (!all(is.na(tasks))) {
    keep <- keep & cells$task %in% tasks
  }
  if (!all(is.na(sizes))) {
    keep <- keep & cells$n %in% as.integer(sizes)
  }
  cells <- cells[keep, , drop = FALSE]
  if (!is.na(max_cells)) {
    cells <- utils::head(cells, max_cells)
  }
  rownames(cells) <- NULL
  cells
}
