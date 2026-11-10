# Digests of small fits, to show that a change leaves every fit bitwise as it was (test-fit-digests.R).
#
# A digest is the md5 of the numbers as little-endian doubles, so -0 and +0 differ:
# - a fit: every family's tree variables, leaf values and leaf intervals (get_model()), then its predictions at the
#   first 40 rows (type "prob" for classification);
# - a purified fit: the predict_components() matrix (columns in the returned order) and intercept; separately the
#   component names; separately predict() after purify().
#
# Every argument whose default differs from master's is given explicitly (max_interaction, n_families, splits, n_try_split_points,
# share_candidates, loss, prune = FALSE), and options(randomPlantedForest.exact_draw = TRUE), so the fits are master's settings;
# the deterministic fits use n_families = 1. Digests hold per platform only (libm and
# std::uniform_int_distribution differ between platforms), so the test pins them for macOS arm64.

port_digest_data <- function() {
  set.seed(120)
  n <- 200
  d <- data.frame(
    x1 = stats::runif(n),
    x2 = stats::rnorm(n),
    x3 = sample(0:1, n, TRUE),
    x4 = round(stats::rnorm(n), 1),
    x5 = sample(c(-1, 0, 0.5, 1), n, TRUE)
  )
  # every other zero of x5 becomes -0, computed: a literal -0 in a function can turn into +0 when it is byte-compiled
  # (the compiler's constant pool merges constants that are identical(), and identical(-0, 0))
  z <- which(d$x5 == 0)[c(TRUE, FALSE)]
  d$x5[z] <- -d$x5[z]
  stopifnot(any(d$x5 == 0 & 1 / d$x5 < 0), any(d$x5 == 0 & 1 / d$x5 > 0))
  d$y <- sin(2 * pi * d$x1) + d$x2 * d$x3 + 0.5 * d$x4 + (d$x5 > 0.25) + stats::rnorm(n, sd = 0.3)
  d$cl <- factor(d$y > stats::median(d$y))
  # three classes (added for the port; uses no random numbers, so the columns above are as in T120's helper)
  d$cl3 <- cut(d$y, stats::quantile(d$y, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("a", "b", "c"))
  d
}

port_md5 <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(as.double(x), f, size = 8, endian = "little")
  unname(tools::md5sum(f))
}

port_md5_text <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(as.character(x), f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

# the numbers of a fit: its model, then its predictions on newdata
port_fit_numbers <- function(fit, newdata, type = NULL) {
  m <- fit$fit$get_model()
  nums <- unlist(lapply(m, function(fam) c(unlist(fam$variables), unlist(fam$values), unlist(fam$intervals))))
  pred <- if (is.null(type)) predict(fit, newdata) else predict(fit, newdata, type = type)
  c(as.double(nums), as.double(as.matrix(pred)))
}

# predict_components() of a purified fit: numbers (matrix, then intercept), names, and predict() after purify()
port_purified_numbers <- function(fit, newdata, type = NULL) {
  # a fresh copy each call: master's predict_components() converts its new_data to a data.table in place (F4)
  nd <- newdata[seq_len(nrow(newdata)), , drop = FALSE]
  pc <- predict_components(fit, nd)
  nd <- newdata[seq_len(nrow(newdata)), , drop = FALSE]
  pred <- if (is.null(type)) predict(fit, nd) else predict(fit, nd, type = type)
  list(
    components = c(as.double(as.matrix(pc$m)), as.double(pc$intercept)),
    names = c(names(pc$m), names(pc$intercept)),
    predict = as.double(as.matrix(pred))
  )
}

# all numbers, by digest name
port_fit_numbers_all <- function() {
  # the exact draw of "leaf" candidates, as master (its default may change)
  old <- options(randomPlantedForest.exact_draw = TRUE)
  on.exit(options(old), add = TRUE)
  d <- port_digest_data()
  nd <- d[1:40, ]
  fo <- y ~ x1 + x2 + x3 + x4 + x5
  fc <- cl ~ x1 + x2 + x3 + x4 + x5
  fc3 <- cl3 ~ x1 + x2 + x3 + x4 + x5
  out <- list()
  # every split structure: regression (random at 1 and 2 threads, deterministic) and binary classification. The
  # digests are named by the structures' names of 0.5.0, so that they compare with earlier versions' digests
  structures <- c(
    leaves = "leaf",
    cur_trees_1 = "tree_then_leaf",
    cur_trees_2 = "tree_then_all_leaves",
    res_trees = "target_tree_then_leaf",
    hist = "leaf_binned"
  )
  for (lab in names(structures)) {
    s <- structures[[lab]]
    for (cfg in list(c(det = 0, nt = 1), c(det = 0, nt = 2), c(det = 1, nt = 1))) {
      set.seed(7)
      f <- rpf(
        fo,
        data = d,
        max_interaction = 3,
        n_families = if (cfg[["det"]] == 1) 1 else 4,
        n_splits = 40,
        split_structure = s,
        n_try_split_points = 5,
        share_candidates = 0.5,
        prune = FALSE,
        nthreads = cfg[["nt"]],
        deterministic = cfg[["det"]] == 1
      )
      out[[sprintf("reg_%s_det%d_nt%d", lab, cfg[["det"]], cfg[["nt"]])]] <- port_fit_numbers(f, nd)
    }
    set.seed(9)
    f <- rpf(
      fc,
      data = d,
      max_interaction = 2,
      n_families = 4,
      n_splits = 30,
      split_structure = s,
      n_try_split_points = 5,
      share_candidates = 0.5,
      prune = FALSE,
      loss = "exponential"
    )
    out[[sprintf("cls_%s", lab)]] <- port_fit_numbers(f, nd, type = "prob")
  }
  # every classification loss, binary and three classes ("leaf")
  for (k in c(2, 3)) {
    for (l in c("L1", "L2", "logit", "exponential")) {
      set.seed(11)
      f <- rpf(
        if (k == 2) fc else fc3,
        data = d,
        max_interaction = 2,
        n_families = 4,
        n_splits = 30,
        split_structure = "leaf",
        n_try_split_points = 5,
        share_candidates = 0.5,
        prune = FALSE,
        loss = l
      )
      out[[sprintf("loss%d_%s", k, l)]] <- port_fit_numbers(f, nd, type = "prob")
    }
  }
  # purified fits: predict_components() and predict() after purify()
  purified <- function(name, fo, type, mi, purify_args = list(), nthreads = 1) {
    set.seed(13)
    f <- rpf(
      fo,
      data = d,
      max_interaction = mi,
      n_families = 4,
      n_splits = 40,
      split_structure = "leaf",
      n_try_split_points = 5,
      share_candidates = 0.5,
      prune = FALSE,
      loss = if (is.null(type)) NULL else "exponential",
      nthreads = nthreads
    )
    do.call(purify, c(list(f), purify_args))
    p <- port_purified_numbers(f, nd, type)
    out[[sprintf("%s_components", name)]] <<- p$components
    out[[sprintf("%s_names", name)]] <<- p$names
    out[[sprintf("%s_predict", name)]] <<- p$predict
  }
  purified("pur_reg", fo, NULL, 3)
  purified("pur_reg_nt2", fo, NULL, 3, list(nthreads = 2), nthreads = 2)
  purified("pur_reg_mode1", fo, NULL, 3, list(mode = 1))
  purified("pur_reg_cap2", fo, NULL, 3, list(maxp_interaction = 2))
  purified("pur_cls", fc, "prob", 2)
  purified("pur_cls3", fc3, "prob", 2)
  out
}

port_fit_digests <- function(numbers = port_fit_numbers_all()) {
  vapply(
    names(numbers),
    function(nm) if (grepl("_names$", nm)) port_md5_text(numbers[[nm]]) else port_md5(numbers[[nm]]),
    character(1)
  )
}
