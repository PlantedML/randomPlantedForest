# The display of prune_components(max_panels = P): every kept component of order >= 3 whose own cut points need more
# than P heat maps is coarsened into its panel form, then the same rule prunes again on the kept set with the
# coarsened components.
#
# The panel form of a component f_S (|S| >= 3; C++ coarsen_select(), src/include/coarsen.hpp): |S| - 2 coarse features
# C are cut into equal-mass intervals of the training rows (q_k intervals, ties merged), with a maximal q vector of
# product <= P (q_k = 1 averages feature k out entirely; with |S| - 2 >= 4 and P = 10 some q_k are always 1); the two
# other features stay at full resolution, so each of the at most P coarse cells is one heat map. Within a cell the
# component is averaged with product-marginal weights: each rectangle's indicator on a coarse axis k is replaced by
# mu_k(B & I) / mu_k(I) for the cell's interval I (mu_k: the training marginal of feature k). Averaging over a coarse
# feature k under mu_k commutes with this (the mean along k of the coarsened component is the coarsened mean along k),
# but the means along the two full-resolution features and the mean over the training rows change. purify() identifies
# the components under the joint distribution of the training rows, not under the product of the marginals, so the
# coarsened component is not purified; it is used as it is (nothing is purified again or moved into the intercept).
# Among all coarse sets C (in lexicographic order) and maximal q vectors (lexicographic, listed directly: never
# P^(|S| - 2) candidates), the one keeping the most of the forest component's variance on the training rows is chosen
# (ties: the first; a later one must be larger by more than 1e-15). A component whose own cut points (the leaves of all
# trees containing S) already need at most P heat maps (the product of its |S| - 2 smallest interval counts) stays as
# it is; its panels (those |S| - 2 features and their own cut points) are kept for a plot ($display$asis), as
# $display$coarse keeps a coarsened one's.
#
# The choice needs the component's values at the training rows (the forest's decomposition, prune_columns()) and its
# rectangles merged over the families (C++ merged_components(): a purification of the families restricted to the
# components being chosen), not the purified store. The components are chosen in batches of prune_coarsen_batch: a
# batch's rectangles and candidates are held only while it is chosen, and the choice of one component does not depend
# on the others.
#
# The second pruning:
# - chosen out of bag: the same candidates, folds and correction, on each family's coarsened component at its
#   out-of-bag rows (coarsening is linear in the rectangles, so every family's rectangles of S are coarsened with the
#   same cells: C++ coarsen_oob(), a purification of the families restricted to the coarsened components); R's
#   coefficients out of bag;
# - a given rule R fitted out of bag: the given rule and tau again, R's coefficients on the coarsened OOB components;
# - otherwise (chosen by CV, a given D, a given R fitted in sample): the chosen or given rule and tau again, R's
#   coefficients in sample (a CV choice would need the display in every fold's refit).
# If no kept component needs coarsening, the model is the first pruning's.
#
# first: the first pruning's model; dec: the forest's decomposition on the training rows; sums: the forest's OOB sums
# (prune_columns(); NULL unless select_by = "oob"); cand: the candidates.
prune_display <- function(first, object, tr, dec, sums, cand, max_panels, nthreads) {
  tick <- function() proc.time()[["elapsed"]]
  secs <- list()
  t0 <- tick()
  n <- length(tr$y)
  # the predictor names as the components are named (prune_columns())
  fnames <- colnames(preprocess_predictors_predict(object, tr$predictors[1, , drop = FALSE]))

  # ---- the kept components of order >= 3 in their panel form
  kept <- names(first$beta)
  kept_sets <- first$sets
  hi <- which(lengths(kept_sets) >= 3)
  sel <- list()
  if (length(hi) > 0) {
    batches <- split(hi, (seq_along(hi) - 1L) %/% prune_coarsen_batch)
    sel <- unlist(
      lapply(batches, function(b) {
        object$fit$coarsen_select(
          lapply(kept_sets[b], as.integer),
          max_panels,
          nthreads,
          dec$own[, kept[b], drop = FALSE]
        )
      }),
      recursive = FALSE,
      use.names = FALSE
    )
  }
  names(sel) <- kept[hi]
  co <- names(sel)[vapply(sel, function(s) s$form == "coarsened", logical(1))]
  secs$coarsen <- tick() - t0
  forms <- data.frame(
    component = kept[hi],
    order = lengths(kept_sets[hi]),
    form = vapply(sel, function(s) s$form, ""),
    forest_panels = vapply(sel, function(s) s$forest_panels, numeric(1)),
    heat_maps = vapply(
      sel,
      function(s) if (s$form == "coarsened") as.numeric(s$cells) else s$forest_panels,
      numeric(1)
    ),
    var_kept = vapply(sel, function(s) if (s$form == "coarsened") s$var_kept else 1, numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(forms) <- NULL
  coarse <- lapply(stats::setNames(co, co), function(nm) {
    s <- sel[[nm]]
    list(
      S = fnames[s$S],
      C = fnames[s$C],
      q = s$q,
      cuts = stats::setNames(s$cuts, fnames[s$C]),
      cells = s$cells,
      var_kept = s$var_kept,
      forest_panels = s$forest_panels,
      candidates = s$candidates,
      stored = s$stored,
      seconds = s$seconds
    )
  })
  # the components kept as they are: their panels, the |S| - 2 predictors C with the fewest own cut points (the leaf
  # bounds of all trees containing S inside the training range; the first of ties) cut at these; as in coarse, S is
  # the whole set and the heat maps are over S minus C
  as_is <- names(sel)[vapply(sel, function(s) s$form == "asis", logical(1))]
  asis <- lapply(stats::setNames(as_is, as_is), function(nm) {
    s <- sel[[nm]]
    list(
      S = fnames[s$S],
      C = fnames[s$C],
      cuts = stats::setNames(s$cuts, fnames[s$C]),
      forest_panels = s$forest_panels
    )
  })

  # ---- the second pruning, on the kept set with the coarsened components
  second <- length(co) > 0
  if (!second) {
    res <- first
  } else {
    t1 <- tick()
    j <- match(co, kept)
    own <- dec$own[, kept, drop = FALSE]
    own[, j] <- vapply(co, function(nm) sel[[nm]]$g, numeric(n))
    dimnames(own) <- NULL
    dec2 <- prune_dec_matrix(own, kept_sets, dec$intercept, n, fnames)
    rule <- first$rule
    tau <- first$tau
    oob <- NULL
    oob_fit <- NULL
    if (identical(first$select_by, "oob")) {
      t2 <- tick()
      specs <- lapply(co, function(nm) list(S = sel[[nm]]$S, C = sel[[nm]]$C, cuts = sel[[nm]]$cuts))
      cs <- object$fit$coarsen_oob(specs, nthreads)
      secs$coarsen_oob <- tick() - t2
      # the forest's OOB sums of the kept components, the coarsened ones replaced
      idx <- match(kept, prune_set_names(lapply(sums$sets, as.integer), fnames))
      half <- function(s) {
        out <- matrix(0, n, length(kept))
        out[, !is.na(idx)] <- s[, idx[!is.na(idx)], drop = FALSE]
        out
      }
      S0 <- half(sums$sum0)
      S1 <- half(sums$sum1)
      S0[, j] <- cs$sum0
      S1[, j] <- cs$sum1
      o2 <- list(
        sets = kept_sets,
        sum0 = S0,
        sum1 = S1,
        intercept0 = sums$intercept0,
        intercept1 = sums$intercept1,
        count = sums$count
      )
      rm(S0, S1, cs)
      choose <- !is.null(first$oob)
      sc <- prune_oob_score(
        o2,
        dec2,
        fnames,
        tr$y,
        if (choose) cand,
        first$fold_id,
        object$params$ntrees,
        deterministic = isTRUE(object$params$deterministic)
      )
      oob_fit <- sc$oob_fit
      if (choose) {
        oob <- sc$table
        choice <- oob[oob[[c(min = "min", "1se" = "one_se")[[first$select]]]], ]
        rule <- choice$rule
        tau <- choice$tau
      }
      rm(sc, o2)
    }
    rf <- prune_fit_rule(rule, tau, dec2, tr$y, oob = oob_fit)
    res <- prune_result(
      object,
      rule,
      tau,
      dec2,
      rf,
      first$cv,
      oob,
      first$select,
      first$select_by,
      first$fold_id,
      first$seed
    )
    secs$second <- tick() - t1
  }

  # the display's columns: a component's form and heat maps (NA for a dropped one of order >= 3 without a form)
  cmp <- res$components
  m <- match(cmp$component, forms$component)
  cmp$form <- ifelse(cmp$order <= 2, "full", NA_character_)
  cmp$form[!is.na(m)] <- forms$form[m[!is.na(m)]]
  cmp$form[cmp$component %in% co] <- "coarsened"
  cmp$heat_maps <- ifelse(cmp$order <= 2, 1, NA_real_)
  cmp$heat_maps[!is.na(m)] <- forms$heat_maps[m[!is.na(m)]]
  res$components <- cmp
  for (nm in co) {
    coarse[[nm]]$kept <- nm %in% names(res$beta)
  }
  for (nm in as_is) {
    asis[[nm]]$kept <- nm %in% names(res$beta)
  }
  first$fit <- NULL
  secs$total <- tick() - t0
  res$display <- list(
    max_panels = max_panels,
    first = first,
    coarse = coarse,
    asis = asis,
    forms = forms,
    second = second,
    seconds = secs,
    # the coarsened components' built indices (prune_display_built())
    cache = new.env(parent = emptyenv())
  )
  res
}

# The components the display chooses the panel forms of in one restricted pass (prune_display()): fewer passes are
# faster, but each component chosen at once holds its rectangles and candidates (about 0.5-1 GB per three-way
# interaction of California housing at order 8, 200 families: all eleven at once peaked 5 GB higher)
prune_coarsen_batch <- 2L

# The built form of coarsened component `nm` (k-d index over its stored rectangles), kept in the display's cache
# environment; rebuilt if missing or if the pointer is gone (a deserialized object)
prune_display_built <- function(object, nm) {
  cache <- object$display$cache
  xp <- if (is.environment(cache)) cache[[nm]] else NULL
  if (is.null(xp) || !rpf_coarse_valid(xp)) {
    xp <- rpf_coarse_build(object$display$coarse[[nm]]$stored)
    if (is.environment(cache)) {
      assign(nm, xp, envir = cache)
    }
  }
  xp
}
