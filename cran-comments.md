## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The note is CRAN's incoming feasibility note ("New submission"). The spell check may also list the names Hiabu and
Mammen (authors of the method's paper, cited in the Description); they are spelled correctly.

## Test environments

* local macOS (arm64), R 4.5.2, Apple clang 21: `R CMD check --as-cran` with the settings of CRAN's incoming checks
  (`_R_CHECK_CRAN_INCOMING_=true` and the `_R_CHECK_*` variables of the incoming checks), and again with
  `_R_CHECK_DEPENDS_ONLY_=true` (tests, examples and vignette without the suggested packages).
* Linux (Ubuntu 24.04, x86_64), R 4.4.2, GCC 14.1: the same two checks, with `-Wall -pedantic`; the notes and the
  warning there were the machine's (no `inconsolata.sty` for the PDF manual, the toolchain's `-march` flag, an R older
  than the one that built the tarball).
* win-builder (R-devel, R-release) and the macOS builder: <!-- to be filled in by the maintainer -->

## Compiled code

* Every source file compiles without a warning under `-Wall -pedantic` with GCC 14 and 16, Apple clang 21 and LLVM
  clang 22 and 23, and without a warning with GCC's link-time optimization (`-flto -Wodr -Wlto-type-mismatch`).
* The test suite and every example run under AddressSanitizer and UndefinedBehaviorSanitizer (Apple clang), and under
  GCC's UndefinedBehaviorSanitizer on Linux, without a report.

## Threads, time and platforms

* `rpf()` uses one thread by default (`nthreads = 1`). Examples, the vignette and the tests that run on CRAN use at
  most 2 threads.
* On the macOS machine above the tests take about 45 s, the examples 3 s in all (each under 1 s) and the vignette 10 s.
* Tests that pin what one example fit does on the development platform (macOS arm64) are skipped on CRAN; the other
  tests check what holds for every fit.

## Reverse dependencies

None on CRAN (this is the package's first submission).
