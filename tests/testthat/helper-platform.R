# The platform the tests' example fits were chosen on. A fit depends on the platform's random integers (libc++'s
# std::uniform_int_distribution on macOS, libstdc++'s on Linux and Windows) and floating point, so an expectation about
# what one example fit does (which components it keeps, whether it coarsens) holds there; elsewhere the tests check
# what holds for every fit. Not on CRAN (NOT_CRAN unset): CRAN's macOS machines have another compiler and R, so their
# fits need not be ours.
on_macos_arm64 <- function() {
  isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))) &&
    Sys.info()[["sysname"]] == "Darwin" &&
    R.version$arch == "aarch64"
}
