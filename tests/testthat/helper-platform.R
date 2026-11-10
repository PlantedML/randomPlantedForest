# The platform the tests' example fits were chosen on. A fit depends on the platform's random integers (libc++'s
# std::uniform_int_distribution on macOS, libstdc++'s on Linux and Windows) and floating point, so an expectation about
# what one example fit does (which components it keeps, whether it coarsens, that two fits which could round an exact
# tie in loss apart have none) holds there; elsewhere the tests check what holds for every fit.
on_macos_arm64 <- function() {
  Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64"
}
