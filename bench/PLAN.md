# Performance benchmark (`bench/`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `bench/` performance regression suite for randomPlantedForest that measures runtime, peak memory and output hashes of rpf's expensive steps across git refs, locally (`quick`) and on the BIPS Slurm cluster (`full`).

**Architecture:** A port of `bench/` from bips-hb/arf PR #76 (head `fffd285`), adapted to rpf. A launcher script run by `rig run` declares the shared dependencies inline and hands off to a Rapp command-line app. Each ref is installed into its own library holding only rpf, every measurement runs in a fresh `callr` process, and all refs of a cell are measured in one job with interleaved rounds. batchtools runs cells locally or on Slurm, and results are plain tables.

**Tech Stack:** R, rig 0.11 inline script metadata, Rapp 0.4, callr, ps, batchtools, testthat (standalone, not the package suite).

**Spec:** `bench/DESIGN.md`. Reference code: arf PR #76, readable with `gh api "repos/bips-hb/arf/contents/bench/<file>?ref=fffd28543ef5b1916ea139c72ec9f976bee4923d" --jq .content | base64 -d`. Port its logic where a task says so, and drop what is arf-specific (CRAN refs, foreach/psock/mirai backends, data.table thread handling, evidence-row variants).

## Global Constraints

- Everything lives under `bench/`. The only change outside it is the line `^bench$` in `.Rbuildignore`. No package code, tests, CI or Makefile changes.
- Do not commit. The maintainer makes all commits. Each task ends with its changes left in the working tree and a short note of what changed.
- Branch: `bench`, based on `origin/dev`. The default base ref is `dev` until 0.6.0 is tagged.
- Measurements and registries return plain data (vectors, data frames), never R6 objects or fitted models.
- All refs of a cell are measured in the same process tree on the same machine, rounds interleaved across refs.
- BLAS and OpenMP threads are pinned to 1 in every measurement child (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `R_DATATABLE_NUM_THREADS` = 1). rpf's `nthreads` comes only from the cell.
- Prose (README, messages) follows the repo owner's style: one sentence per line, no semicolons joining clauses.
- R code is formatted with `air format bench/`.

## Review Focus

1. A ref fails to install (compile error, missing dependency): the run stops with an error naming the ref, and no measurement ever loads the session's rpf instead. Test in Task 2.
2. A cell needs a newer API than a ref has (e.g. `n_splits = NULL` on `v0.5.0`): that ref is skipped for that cell with a reason in the results, not an error. Test in Task 4.
3. No memory accounting (macOS without `ps`, a container without cgroup v2): the run warns once, records time with `metric = "none"` and `NA` memory, and the report says memory was not measured. Test in Task 3.
4. A dirty working tree: the run refuses unless the ref under test is `WORKTREE`, which measures the uncommitted state on purpose. Test in Task 2.
5. Only one distinct ref after deduplication (a branch that touches no package file): the run measures it, says no comparison is possible, and the report has no delta rows instead of failing. Test in Task 6.

---

### Task 1: Scaffold, launcher and test runner

**Files:**
- Create: `bench/rpf-bench` (launcher, executable), `bench/cli.R` (Rapp app), `bench/tests.R`, `bench/.gitignore`
- Modify: `.Rbuildignore` (add `^bench$`)

**Interfaces:**
- Produces: `rig run bench/rpf-bench <subcommand> [options]` from the repo root. `bench/cli.R` sources the other `bench/*.R` files, and later tasks add subcommands to it. `rig run bench/tests.R` runs the suite's own tests.

- [ ] **Step 1: Write `bench/rpf-bench`.** Shebang `#!/usr/bin/env -S rig run`, a `# /// script` block declaring `R = ">= 4.4"` and the shared dependencies: `Rapp`, `callr`, `ps`, `batchtools`, `rlang`, `Rcpp`, `checkmate`, `cli`, `data.table`, `hardhat`. Body: stop unless `bench/cli.R` exists relative to the working directory ("Run from the repository root"), then `Rapp::run("bench/cli.R", commandArgs(TRUE))`. rpf's `Imports` and `LinkingTo` must be covered so ref libraries can hold rpf alone.
- [ ] **Step 2: Write `bench/cli.R`** as a Rapp app (see arf's `arf-bench` for the `#|` annotation style) with `tier` (`quick`/`full`, short `t`) and `cluster` (`local`/`slurm`, short `c`) options and a `switch(command <- "", ...)` holding a placeholder `info` subcommand that prints `.libPaths()` and R's version.
- [ ] **Step 3: Write `bench/tests.R`** with its own `# /// script` block (the launcher's dependencies plus `testthat`), sourcing `bench/*.R` except `cli.R`, and a first test `test_that("bench runs in its own library", expect_match(.libPaths()[1], "rig"))`.
- [ ] **Step 4: `bench/.gitignore`** ignoring `lib/`, `results/`, `registry/`, `worktrees/`. Add `^bench$` to `.Rbuildignore`.
- [ ] **Step 5: Verify.** `chmod +x bench/rpf-bench && rig run bench/rpf-bench info` prints a rig script library first. `rig run bench/rpf-bench --help` lists `info`. `rig run bench/tests.R` passes. `R CMD build .` produces a tarball without `bench/` (`tar tzf randomPlantedForest_*.tar.gz | grep -c bench/` is 0).
- [ ] **Step 6: Checkpoint.** Leave changes uncommitted and report.

### Task 2: Refs, installation and isolation

**Files:**
- Create: `bench/refs.R`
- Test: `bench/tests.R`

**Interfaces:**
- Produces:
  - `bench_parse_ref(spec) -> list(type = "git" | "worktree", rev = character(1))`. `"WORKTREE"` means the working tree with uncommitted changes, everything else is a git rev (`HEAD`, `dev`, `v0.5.0`, a SHA, a branch).
  - `bench_ref_id(spec) -> character(1)`: content id over `DESCRIPTION`, `NAMESPACE`, `R/`, `src/` (excluding `*.o`, `*.so`), `inst/`. For git revs from `git rev-parse <rev>:<path>`, for `WORKTREE` from file hashes.
  - `bench_resolve_refs(specs) -> character`: specs with duplicate content ids dropped (first kept), with a message naming what was dropped.
  - `bench_install_ref(spec, root = "bench/lib") -> list(spec, id, lib, version, commit)`: installs into `<root>/<id>`, reusing an existing complete install. Git revs go through `git worktree add --detach bench/worktrees/<id> <rev>` and `R CMD INSTALL --no-test-load -l <lib>`, then `git worktree remove`. Stops with the ref name and the install log tail on failure.
  - `bench_assert_lib(lib)`: in any process, stops unless `find.package("randomPlantedForest")` resolves inside `lib`.
  - `bench_assert_clean_tree(refs)`: stops on modified tracked files unless `"WORKTREE"` is among `refs`.
- Port from arf `refs.R`, dropping CRAN refs.

- [ ] **Step 1: Write failing tests** against a throwaway git repo built in `tempdir()` that holds a minimal R-only package named `randomPlantedForest` (DESCRIPTION, NAMESPACE exporting `rpf`, `R/rpf.R` returning its arguments). Tests:
  - `bench_parse_ref("WORKTREE")$type == "worktree"`, `bench_parse_ref("dev")$type == "git"`.
  - Two commits that differ only in `README.md` give the same `bench_ref_id()`, a change in `R/` gives a different one.
  - `bench_resolve_refs(c("HEAD", "HEAD~1"))` returns one spec when only `README.md` changed.
  - `bench_install_ref()` puts the package in its own library, and a second call reuses it (install log not rewritten).
  - A commit with a syntax error in `R/` makes `bench_install_ref()` fail with an error that contains the ref spec (Review Focus 1).
  - `bench_assert_lib()` errors in a `callr::r()` child whose `.libPaths()` lacks the ref library (Review Focus 1).
  - With a modified tracked file, `bench_assert_clean_tree(c("HEAD", "dev"))` errors and `bench_assert_clean_tree(c("WORKTREE", "dev"))` does not (Review Focus 4).
- [ ] **Step 2: Run `rig run bench/tests.R`**, expect these tests to fail with "could not find function".
- [ ] **Step 3: Implement `bench/refs.R`** with the interfaces above.
- [ ] **Step 4: Run `rig run bench/tests.R`**, expect all to pass. Then install real refs once: `bench_install_ref("dev")` and `bench_install_ref("v0.5.0")` from a script, both succeed (each compiles rpf, a few minutes).
- [ ] **Step 5: Checkpoint.**

### Task 3: Measurement in a fresh process

**Files:**
- Create: `bench/measure.R`
- Test: `bench/tests.R`

**Interfaces:**
- Produces:
  - `bench_metric() -> "cgroup" | "rss" | "none"`. `"cgroup"` when the process has a cgroup v2 directory with `memory.stat` that is dedicated to this run (inside a Slurm job, or inside a scope this suite started). `"rss"` on macOS when `ps` works. `"none"` otherwise, with one warning per run that accurate memory measurement needs cgroup v2 on Linux.
  - `bench_ensure_scope(args) -> invisible()`: on Linux outside Slurm and outside a scope this suite started, re-executes the CLI under `systemd-run --user --scope --quiet --setenv=RPF_BENCH_SCOPED=1` when `systemd-run` exists, so local cgroup numbers are not polluted by the desktop session. Does nothing elsewhere.
  - `bench_measure(lib, op, cell, data_path, metric) -> list(time = numeric(cell$iters), peak_mb, floor_mb, digest, error)`: runs in a `callr::r()` child with `libpath = c(lib, .libPaths())` and the pinned thread variables, asserts the library (Task 2), loads the data, runs `op$setup` untimed, then times `cell$iters` calls of `op$run` and returns the hash of the last output. Memory: the peak over a sampler at a 5 ms interval during one extra untimed call, minus `floor_mb`, the same sampler's level after setup. For `"cgroup"` the sample is `anon + shmem` of the cgroup (port arf's `.bench_cgroup_anon_kb` and `.bench_cgroup_dir`), for `"rss"` the child's RSS via `ps::ps_memory_info()`. Errors inside the child come back in `error`, not as a crash.
  - `bench_digest(x) -> character(1)`: `rlang::hash()` of `x` after dropping attributes that carry environments or timings.

- [ ] **Step 1: Write failing tests.** `bench_measure()` on the fake package from Task 2 with an op that allocates a 200 MB numeric vector returns `peak_mb - floor_mb` between 150 and 400 when `bench_metric() != "none"`, skipped otherwise. With `metric = "none"`, `peak_mb` and `floor_mb` are `NA` and `time` has `iters` entries (Review Focus 3). An op that calls `stop("boom")` returns `error` containing "boom". `bench_digest()` is equal for two identical numeric vectors and differs after changing one element.
- [ ] **Step 2: Run, expect failures.**
- [ ] **Step 3: Implement `bench/measure.R`.** The sampler runs in the parent and polls while the child works (`callr::r_bg()`), as arf's `.bench_run_sampled` does.
- [ ] **Step 4: Run, expect passes** on Linux. On macOS, the cgroup test is skipped and the RSS path is exercised.
- [ ] **Step 5: Checkpoint.**

### Task 4: Operations, data and cells

**Files:**
- Create: `bench/ops.R`, `bench/cells.R`
- Test: `bench/tests.R`

**Interfaces:**
- Produces:
  - `bench_make_data(task = c("regr", "binary", "multiclass"), n, p, seed = 1L) -> list(train = data.frame, test = data.frame)`: `p` numeric predictors from `runif`, plus one factor predictor with 5 levels, an outcome `y` from a fixed function with main effects and one two-way and one three-way interaction, Gaussian noise for `regr`, and classes by thresholding (binary) or binning into 3 classes (multiclass). `test` has 1,000 rows, enough to pass rpf's 512-row leaf-index threshold in `predict()`.
  - `bench_ops() -> named list` of ops, each `list(min_version, setup = function(data, cell), run = function(state, cell))`. Ops: `fit_fixed` (`rpf(n_splits = cell$n_splits, prune = FALSE)`, output: predictions on `test`), `fit_oob` (`n_splits = NULL, prune = FALSE`), `fit_default` (the default call), `purify`, `predict` (setup fits, run predicts `test`), `predict_components`. `min_version` is `"0.5.0"` for `fit_fixed`, `purify`, `predict` and `predict_components`, and `"0.5.0.9000"` for `fit_oob` and `fit_default`.
  - `bench_rpf_args(args, version) -> list`: for versions below `0.5.0.9000`, renames `n_families` → `ntrees`, `n_splits` → `splits`, `n_try_split_points` → `split_try`, `share_candidates` → `t_try`, and maps `split_structure` values back (`"leaf"` → `"leaves"`, `"leaf_binned"` → `"hist"`, `"tree_then_leaf"` → `"cur_trees_1"`, `"tree_then_all_leaves"` → `"cur_trees_2"`, `"target_tree_then_leaf"` → `"res_trees"`). Unchanged otherwise.
  - `bench_cells(tier = c("quick", "full")) -> data.frame` with columns `op, task, loss, n, p, max_interaction, n_families, n_splits, split_structure, nthreads, iters, mem_reps, tier`. `loss` is `NA` for `regr` and passed to `rpf()` otherwise. `quick`: ops `fit_fixed`, `fit_oob`, `predict`, `purify`, tasks `regr` and `binary`, `n` 1,000 and 5,000, `p` 5, `max_interaction` 2, `loss` the default (`NA`, so `rpf()` picks it), `n_families` 20, `n_splits` 30, `split_structure` `"leaf"`, `nthreads` 1, `iters` 3, `mem_reps` 3, sized so the tier finishes in about 10 minutes for three refs on a laptop. `full`: all ops and tasks, `loss` `"logit"` and `"exponential"` for classification, `n` 10,000 and 50,000, `p` 10, `max_interaction` 2 and 3, `n_families` 50 and 200, `n_splits` 100, every `split_structure` for `fit_fixed`, `nthreads` 1, 4 and 16, `iters` 5, `mem_reps` 5. Grids beyond these are the builder's call.
  - `bench_select_cells(cells, ops = NA, tasks = NA, sizes = NA, max_cells = NA) -> data.frame`.
  - `bench_ref_supports(op, version) -> logical(1)`.
- Port the structure of arf `cells.R`, with rpf's dimensions.

- [ ] **Step 1: Write failing tests.** `bench_make_data("multiclass", 200, 4)$train$y` is a factor with 3 levels and `nrow(...$test) == 1000`. `bench_rpf_args(list(n_families = 5, split_structure = "leaf"), "0.5.0")` gives `list(ntrees = 5, split_structure = "leaves")` and is unchanged for `"0.5.0.9000"`. `bench_ref_supports("fit_oob", "0.5.0")` is `FALSE` and `bench_ref_supports("fit_fixed", "0.5.0")` is `TRUE` (Review Focus 2). Every row of `bench_cells("quick")` and `bench_cells("full")` names an op in `bench_ops()`. `bench_select_cells(bench_cells("quick"), ops = "predict")` has only `predict` rows.
- [ ] **Step 2: Run, expect failures.**
- [ ] **Step 3: Implement `bench/ops.R` and `bench/cells.R`.**
- [ ] **Step 4: Run, expect passes.** Then run each op once against the installed `dev` ref through `bench_measure()` on a quick cell and check that none returns an `error`.
- [ ] **Step 5: Checkpoint.**

### Task 5: One cell across refs

**Files:**
- Create: `bench/run-cell.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_install_ref()` results (Task 2), `bench_measure()` and `bench_metric()` (Task 3), `bench_ops()`, `bench_make_data()`, `bench_ref_supports()`, `bench_rpf_args()` (Task 4).
- Produces:
  - `bench_schema() -> character`: `ref, rpf_version, commit, op, task, loss, n, p, max_interaction, n_families, n_splits, split_structure, nthreads, iters, metric, peak_mb, floor_mb, peak_delta_mb, mem_reps, mem_limit_mb, time_median, time_min, time_max, digest, skipped, error, bound_mb, tier, host, kernel, r_version, job_id, timestamp`.
  - `bench_run_cell(cell, refs, data_dir) -> data.frame` with exactly `bench_schema()` columns and one row per ref: writes the cell's data once to `data_dir`, then runs `cell$mem_reps` rounds, each measuring every supporting ref once in a shuffled order. `time_*` come from all timed calls of all rounds, `peak_delta_mb` is the median over rounds of `peak_mb - floor_mb`. Unsupported refs get `skipped = "needs rpf >= <min_version>"` and `NA` measurements. `bound_mb` is the estimate from the ref's `rpf_resources()` if its namespace exports one (maximum over the returned steps, in MB), otherwise `NA` (Claude's call on the name, to be matched in the estimator's PR).

- [ ] **Step 1: Write failing tests** with the fake package installed as two refs whose `rpf()` differ in an allocation size. One quick cell gives two rows with all `bench_schema()` columns, different `digest` when outputs differ, `skipped` set and `NA` time for a ref whose version is below the op's `min_version` (Review Focus 2), and `bound_mb` filled only for the ref whose fake package exports `rpf_resources()`.
- [ ] **Step 2: Run, expect failures.**
- [ ] **Step 3: Implement `bench/run-cell.R`.**
- [ ] **Step 4: Run, expect passes.**
- [ ] **Step 5: Checkpoint.**

### Task 6: Collation and report

**Files:**
- Create: `bench/collate.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: rows with `bench_schema()` columns (Task 5).
- Produces:
  - `bench_deltas(rows, baseline = "dev") -> data.frame`: per cell key (all cell columns) and non-baseline ref, `time_pct` and `mem_pct` against the baseline, verdicts `faster`/`slower`/`inconclusive` (time threshold 10%) and `less`/`more`/`same` (memory threshold 3%), `output_differs` (digest differs from the baseline), `bound_ok` (`peak_delta_mb <= bound_mb`, `NA` without a bound). Port arf's `bench_deltas`.
  - `bench_report(rows, baseline = "dev") -> character`: Markdown with a header (refs, metric, tier, host), the delta table, a list of `OUTPUT DIFFERS` cells, a list of cells where `bound_ok` is `FALSE`, skipped cells with reasons, and the line "Memory was not measured (no cgroup v2 or RSS)" when the metric is `"none"`.
  - `bench_write_results(rows, dir = "bench/results") -> character` (paths written: `rows-<timestamp>.csv`, `report.md`).
  - `bench_history_append(rows, anchors, path = "bench/history.csv")`: appends rows of anchor refs, keeping the newest per cell and ref.

- [ ] **Step 1: Write failing tests** on hand-made rows: a 20% slower time gives `slower`, 5% gives `inconclusive`, a digest mismatch appears under `OUTPUT DIFFERS`, `peak_delta_mb = 120, bound_mb = 100` appears under the bound violations. Rows from a single ref give a report with no delta rows and the sentence "Only one distinct ref was measured, so there is nothing to compare." (Review Focus 5). Rows with `metric = "none"` give the "Memory was not measured" line (Review Focus 3).
- [ ] **Step 2: Run, expect failures.**
- [ ] **Step 3: Implement `bench/collate.R`.**
- [ ] **Step 4: Run, expect passes.**
- [ ] **Step 5: Checkpoint.**

### Task 7: Registry, submission and the CLI

**Files:**
- Create: `bench/registry.R`
- Modify: `bench/cli.R` (replace `info` with the real subcommands)
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `bench_make_registry(dir, cluster = c("local", "slurm"))`: batchtools registry with `source` set to `bench/*.R` (absolute paths, not `cli.R`). `local` uses `makeClusterFunctionsInteractive()` so cells run one after another. `slurm` loads the site configuration (`/etc/xdg/batchtools/config.R`).
  - `bench_submit(reg, cells, refs, resources = list())`: one job per cell, each calling `bench_run_cell()`, results saved as data frames. Per-cell resources on Slurm: `ncpus = max(cell$nthreads, 1) + 1`, memory from a conservative per-op default (builder's call), walltime default 8 h.
  - `bench_collect(reg) -> data.frame` of all finished jobs' rows.
  - `bench_resubmit(reg, mem_factor = 2)`: resubmits expired or failed jobs with memory times `mem_factor`.
  - CLI subcommands in `bench/cli.R`: `run` (options `--refs` default `WORKTREE,dev` with anchors from `bench/anchors.csv` added on `full`, `--ops`, `--tasks`, `--sizes`, `--max-cells`, `--iters`, `--mem-reps`), `collect`, `resubmit`, `report`. `run` calls `bench_ensure_scope()`, checks the tree and dependencies, installs refs, resolves duplicates, and on `local` waits and writes the report. On `slurm` it submits and exits.
  - `bench/anchors.csv` with one column `ref` and the row `v0.5.0`.

- [ ] **Step 1: Write a failing test**: a local registry in `tempdir()` with two fake-package refs and two quick cells, `bench_submit()` then `bench_collect()` gives 4 rows with `bench_schema()` columns.
- [ ] **Step 2: Run, expect failure.**
- [ ] **Step 3: Implement `bench/registry.R` and the CLI subcommands.**
- [ ] **Step 4: Run, expect pass.** Then end to end: `rig run bench/rpf-bench run --refs WORKTREE,dev,v0.5.0` finishes on the quick tier in under 30 minutes, writes `bench/results/report.md`, shows `fit_oob` skipped for `v0.5.0`, and shows no `OUTPUT DIFFERS` between `WORKTREE` and `dev` on an unmodified tree (they are deduplicated into one ref then, which the report states).
- [ ] **Step 5: Checkpoint.**

### Task 8: README and handover

**Files:**
- Create: `bench/README.md`

- [ ] **Step 1: Write `bench/README.md`**: requirements (rig 0.11, Linux with cgroup v2 for memory, `systemd-run` for clean local numbers, macOS RSS fallback), the two tiers, how to compare a PR branch against `dev` (`rig run bench/rpf-bench run --refs <branch>,dev`), how to run `full` on the BIPS cluster (`-t full -c slurm run`, then `collect`, `resubmit --mem-factor 2`, `report`), how to read a report (verdict thresholds, `OUTPUT DIFFERS`, bound violations, skipped cells), and where results go. Point to `DESIGN.md` for the why.
- [ ] **Step 2: Verify** every command in the README by running it (the `slurm` ones only to the point of a dry `--help`).
- [ ] **Step 3: Checkpoint.** Report what is ready for the maintainer to validate on the cluster: the `full` tier end to end, and Slurm memory and walltime defaults.
