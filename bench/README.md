# Performance benchmark

`bench/` measures runtime, peak memory and output of randomPlantedForest's expensive steps for several versions ("refs") of the package and compares them.
It answers questions like "is this pull request really faster, and does it still give the same fits?".
Why it is built this way is in [DESIGN.md](DESIGN.md).

It is not part of the package: `bench/` is in `.Rbuildignore`.

## Requirements

- [rig](https://github.com/r-lib/rig) 0.11 or later. The runner declares its dependencies inline, and `rig run` installs them into a cached library of its own on first use.
- Run everything from the repository root.
- For memory: Linux with cgroup v2. Local runs re-start themselves inside a `systemd-run --user --scope`, so other programs don't count. On macOS, peak RSS is measured instead. Without either, only time is measured and the report says so.

## Comparing refs locally

```sh
rig run bench/rpf-bench run --refs my-branch,dev
```

A ref is a branch, tag or commit, or `WORKTREE` for the working tree with uncommitted changes.
A ref that is not a local branch is looked up on `origin`, so `dev` works without a local checkout.
Refs with identical package files are measured once.
Without `--refs` the run compares `WORKTREE` with `dev`.
`--baseline` (default `dev`) sets the ref the others are compared against and must be one of the refs.

The quick tier runs 16 cells on simulated data and takes a few minutes per ref.
Narrow it down with `--ops`, `--tasks`, `--sizes` and `--max-cells`, see `rig run bench/rpf-bench run --help`.
An unknown value in any of them stops the run before anything is installed.
`--dry-run` prints the selected cells and refs and exits.

Each ref is installed once into `bench/lib/<id>` (compiling rpf takes one to two minutes) and reused afterwards.
A run refuses to start when package files have uncommitted changes, unless `WORKTREE` is one of the refs.

## The full tier on the BIPS cluster

```sh
rig run bench/rpf-bench run -t full -b slurm
# once the jobs are done:
rig run bench/rpf-bench collect
# for jobs that ran out of memory (or time: --walltime-factor), list them, then submit:
rig run bench/rpf-bench resubmit --mem-factor 2
rig run bench/rpf-bench resubmit --mem-factor 2 --yes
```

Each cell is one Slurm job, which measures all refs one after another so that comparisons never cross nodes.
The full tier adds the anchors in `anchors.csv` (released versions) to the refs, and `collect` appends their rows to `history.csv` for a trend across releases.
Jobs run the R that submitted them and find rig's library through `R_LIBS`, which Slurm passes on by default. `/etc/xdg/batchtools/config.R` provides the site configuration.

## Reading a report

`run` (locally) and `collect` write `bench/results/report.md` and the raw rows as CSV.
`rig run bench/rpf-bench report` prints the newest report again.

```sh
rig run bench/rpf-bench report --format html --baseline v0.5.0
```

renders `bench/report.qmd` with Quarto into a single page under `bench/results/`: the text report, one tab per op with boxplots of every timed call and every round's memory per cell and ref, and scaling over the number of rows.
It needs Quarto on the `PATH` (on the BIPS cluster: `module load quarto`).
`--rows` picks another rows file, for example one copied from the cluster.
Results collected before the raw values were kept show only median, minimum and maximum, drawn as crosses without a box.

- Time verdicts: `faster` or `slower` from a change of 10% in the median time, `inconclusive` below that.
- Memory verdicts: `more` or `less` from a change of 3% in peak memory above the level after setup, `same` below that.
- Memory `below resolution`: the operation used less than 1 MB beyond the level after setup. This is common for `predict`, `purify` and `predict_components`: their setup fits a forest, and they reuse memory the allocator kept from that fit. Very short operations can also fall between two samples (every 5 ms). Their memory is not measurable this way.
- For the same reason, the memory above setup can even fall as `n` grows. The HTML report therefore plots the whole peak of the R process, with the level after setup marked, and the report itself explains how to read each number.
- `n` is the number of training rows. `predict` and `predict_components` predict `predicted` rows (1000 unless the cell says otherwise), and the full tier varies that number for them separately.
- `OUTPUT DIFFERS`: the predictions or components differ from the baseline's. That is expected for changes that alter fits on purpose, and a red flag for changes that claim identical fits.
- Memory bound exceeded: the ref's own estimate (`rpf_resources()`, once rpf has one) was lower than the measured memory.
- Skipped: the ref's `rpf()` lacks an argument the operation needs, for example `fit_oob` needs `oob_path`, which came with `n_splits = NULL`.

## Tests

```sh
rig run bench/rpf-bench test
```

The suite's own tests use a small fake package, so they run without compiling rpf.
