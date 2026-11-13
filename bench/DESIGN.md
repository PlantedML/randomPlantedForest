# Performance benchmark (`bench/`): design

Status: design, not yet built.
Written 2026-11-13.

## Why

randomPlantedForest changes fast, and many pull requests claim speed or memory gains with "bitwise the same fits".
Checking such a claim should take minutes on a development machine, not a cluster setup.

`bench/` is a performance regression test.
It measures runtime and peak memory of rpf's expensive steps for several versions of the package in the same session, and it checks that their outputs are identical where they should be.

It has a second purpose.
rpf should be able to say, before fitting, how much memory an `rpf()` call will need at most.
That estimate is derived analytically and lives in the package.
`bench/` checks it against measured memory, so wrong assumptions in the math show up as failing cells.

Predictive performance (OpenML suites, comparison with other learners) is out of scope here.
It is a separate project, a rewrite of the PlantedML `rpf_benchmark` repo.

## Prior art

The design follows `bench/` in bips-hb/arf PR #76.
Similar benchmark setups exist for arf and xplainfi, and shared tooling across the three packages is a long-term soft goal.
Code here should keep the generic machinery (installing refs, running measurements, sampling memory, collecting, reporting) apart from the rpf-specific parts (cells, data, the operations measured), so it could be extracted later.
No shared layer is built now, and `bench/` may drift from arf's version where rpf needs it.

## Goals

- Compare versions ("refs") of rpf: the working tree or a branch (e.g. an open PR), the `dev` branch, and an anchor, a tagged release.
- Measure per cell and ref: runtime, peak memory, and a hash of the output.
- Report deltas against a base ref with a verdict per cell, and flag cells whose outputs differ between refs.
- Record the package's predicted memory bound next to the measured peak, and flag cells where the measurement exceeds the bound.
- Two tiers:
  - `quick`: runs on the developer's machine (Linux or macOS) in at most about 30 minutes, ideally a few minutes.
  - `full`: runs on the BIPS cluster (Slurm, via batchtools) and covers the whole grid.
- Usable by the maintainers (Lukas Burk, Munir Hiabu) and by coding agents working in this repo.

## Non-goals

- Predictive performance, tuning, comparison with other learners.
- CI integration, for now.
- A third tier or a third machine context: the two tiers cover local development and the cluster.

## Anti-goals

- Interfering with ongoing package development.
  `bench/` lives in its own pull request, is listed in `.Rbuildignore`, and touches no package code, tests, CI or Makefile targets.
  The memory estimator is package code and comes in a separate pull request.
- Storing R6 objects or fitted models in registries or results: every measurement returns plain data.
- Comparisons across machines or cluster nodes: all refs of a cell are measured in the same job.

## How it works

These parts are decided.

**Refs and isolation (as in arf).**
A ref is a git ref (branch, tag, commit, or the working tree).
rpf is not on CRAN, so anchors are git tags.
Each ref is installed into its own library containing only rpf, from a `git worktree`.
All dependencies come from one shared library, so refs differ only in rpf.
Refs with identical package content are measured once.
A run refuses to start on a dirty tracked tree, unless the working tree itself is the ref under test (builder's choice how to mark that).

**Runner environment.**
The runner is a command-line tool built with Rapp, as in arf.
Its dependencies, including Rapp, are declared inline with rig 0.11's script metadata (`# /// script` block, `rig run`), which gives a locked, cached environment.
That environment is the shared library for all refs.
rig is available on the BIPS cluster, so local and cluster runs use the same mechanism.

**Measurement.**
Every measurement runs in a fresh R process (`callr`).
Rounds are interleaved across refs within a cell, and setup (data generation, loading) is not timed.
Runtime is the median of several timed calls.
Memory is sampled from the cgroup v2 memory accounting of the job, as in arf, and reported as peak minus a baseline taken in the same round.
On macOS, peak RSS of the measurement process is the fallback if it is simple to obtain.
Otherwise the runner warns that accurate memory measurement needs cgroup v2 on Linux and measures time only.
Each cell hashes its output (predictions, components) per ref.

**Threads.**
BLAS and OpenMP threads are pinned to 1.
rpf's own `nthreads` is a cell dimension, not ambient state.

**Tiers.**
`quick` runs cells one after another on the local machine.
`full` submits one batchtools job per cell to Slurm, using the site configuration.
Cells on one machine are not run concurrently, since they would disturb each other's timing and memory.
On the cluster, jobs may share nodes with other jobs, as in arf: every ref of a cell is measured in the same job with interleaved rounds, so shared-node noise affects all refs alike, and the cgroup memory accounting is per job (Claude's call).
The `full` tier is validated on the cluster by the maintainers, since the builder has no cluster access (Claude's call).

**Results.**
Raw results are plain tables (CSV or similar), gitignored.
Results for anchor refs can be appended to a committed history file for a trend across releases, as in arf.

## What is measured

The operations, on simulated data:

- growing the forest with a fixed `n_splits` and `prune = FALSE`, per `split_structure`
- growing with `n_splits = NULL` (out-of-bag choice), `prune = FALSE`
- the default call (out-of-bag `n_splits`, pruning, display)
- `purify()`
- `predict()` on an unpurified forest, including row counts above the leaf-index threshold
- `predict_components()`
- classification with `loss = "logit"` and `"exponential"`, binary and multiclass

Cell dimensions include n, p, `max_interaction`, `n_families`, `n_splits` and `nthreads`.
`quick` uses a small grid, `full` the whole one.

## The memory estimator (separate pull request)

An experimental function in rpf (`lifecycle::badge("experimental")`) that takes the arguments of `rpf()`, data and settings, and returns an upper bound on peak memory without fitting.
A time estimate may follow.
Its shape follows the futureverse resources proposal (<https://www.futureverse.org/roadmap/resources.html>), `resources(fn) <- function(...) memory(...)`, so that once futureverse implements it, registering the estimator as rpf's resource method is a one-line change and rpf can serve as a showcase.

Known difficulties:

- Memory is driven by the number of leaves.
  With a given `n_splits` the bound can be tight.
  With `n_splits = NULL` the only bound before fitting comes from the cap of 25,600 splits and is loose.
- The steps peak separately (growing, out-of-bag tracking, pruning, display), so the bound is the largest of them, not their sum.
- Allocator overhead: measured peak RSS can far exceed live data (see issue #129, `adaptive_pause_gb`).
  The bound needs a factor for this, calibrated with `bench/`.

The contract `bench/` relies on (Claude's calls, to be confirmed in the estimator's pull request):

- The bound is the additional peak memory of the call in bytes, above the state just before the call, with the input data already in memory.
  This is what `bench/` measures (peak minus a baseline taken after setup), so both sides compare the same quantity.
- It returns the peak per step (growing, out-of-bag tracking, pruning, display) and their maximum, so a failing check points at the step whose model is wrong.
- It covers `rpf()` first. `purify()` and `predict()` may follow.
- `bench/` detects it by name in each ref's namespace, records it for refs that provide it, and skips the check for refs that don't.
- Results keep the measured/bound ratio per cell and step, which is the data for calibrating the allocator factor.

## Open points

- **Anchor compatibility.**
  0.5.0 has the old argument names and lacks `prune`, `n_splits = NULL` and other features, so cells written for the current API cannot all run on it.
  Proposal (builder's call to refine): each cell declares the minimum version it needs, refs below it skip the cell, and the first release with the current API (0.6.0) becomes the main anchor.
  Until it is tagged, the `dev` branch serves as the base: development now goes through `dev`, and `master` gets fewer direct merges.
- **Expected output differences.**
  Some changes alter fits on purpose (e.g. the fast candidate draw).
  The report flags output differences, and whether a difference is expected is noted by the person reading it, not encoded in the suite.

## Left to the builder

- File layout under `bench/`, CLI subcommands and their names (arf has run, plan, collect, resubmit, viz, diag).
- The exact grids for both tiers, within the `quick` time budget.
- Data generators for regression and classification.
- Number of timing iterations and memory rounds, the sampling interval, verdict thresholds (arf starts from 3% for memory and 10% for time).
- Report format (Markdown table, Quarto page, or both).
- Walltime and memory requests per cell on the cluster.
