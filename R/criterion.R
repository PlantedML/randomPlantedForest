# The out-of-bag criterion of a forest, per loss. Every out-of-bag choice (the adaptive stop along the splits, the
# pruning's rule and tau, the tuner) scores a candidate by the mean, over the training rows with an out-of-bag (OOB)
# family, of its per-row criterion, corrected for the number of OOB families.
#
# - Regression, and classification with L1 / L2: the squared error of the score, summed over the score columns (the
#   Brier score of the raw score for classification: one 0/1 column for a binary outcome, the one-hot columns for K
#   classes). With one column this is the regression criterion, and the regression code runs (oob_score_y()).
# - Classification with logit / exponential: the log-loss of the probabilities that predict(type = "prob") makes:
#   binary p = plogis(g); multiclass logit the softmax of the K - 1 scores and the first level's score 0; multiclass
#   exponential the softmax of the K scores.
#
# The correction. An OOB mean averages m families, the forest all B, so it scatters more; to second order the OOB
# criterion exceeds the forest's by (1/m - 1/B) tr(H V) / 2, with H the Hessian of the per-row criterion in the score
# and V the across-family covariance, estimated from the halves of even and odd families: V = d d' / (1/m0 + 1/m1),
# d the difference of their means. For the squared error H = 2I (the regression correction); for the log-loss
# H = diag(p) - pp' over the free score columns. A row whose OOB families are all even or all odd gets, where the
# caller imputes (the pruning), its own tr(H Vbar) / 2 with Vbar the mean of the other rows' V; otherwise (the path, the
# stop) none, as in regression. The adaptive fit's C++ bound uses the same criterion (oob_row_error()).

# The criterion of a fit of `mode` and `loss` with the outcome y (numeric for regression, the factor for
# classification): a list of class "rpf_crit" with type ("mse", "brier", "logloss"), link ("identity", "sigmoid",
# "softmax_ref", "softmax"), V (score columns), K (classes), y (one 0/1 or numeric column, else NULL), cls (the class
# indices), Y (the 0/1 score columns of L1 / L2 multiclass), centre ("none", "ref", "mean"; crit_centred_var()) and n
oob_criterion <- function(mode, loss, y) {
  if (mode == "regression") {
    return(structure(
      list(
        type = "mse",
        link = "identity",
        V = 1L,
        K = 1L,
        y = as.numeric(y),
        cls = NULL,
        Y = NULL,
        centre = "none",
        n = length(y)
      ),
      class = "rpf_crit"
    ))
  }
  if (!is.factor(y)) {
    cli::cli_abort("A classification criterion needs the factor outcome.", .internal = TRUE)
  }
  K <- nlevels(y)
  cls <- as.integer(y)
  if (loss %in% c("L1", "L2")) {
    Y <- if (K == 2) matrix(as.numeric(cls == 2L), ncol = 1) else outer(cls, seq_len(K), "==") * 1
    return(structure(
      list(
        type = "brier",
        link = "identity",
        V = ncol(Y),
        K = K,
        y = if (K == 2) Y[, 1] else NULL,
        cls = cls,
        Y = Y,
        centre = "none",
        n = length(cls)
      ),
      class = "rpf_crit"
    ))
  }
  structure(
    list(
      type = "logloss",
      link = if (K == 2) {
        "sigmoid"
      } else if (loss == "logit") {
        "softmax_ref"
      } else {
        "softmax"
      },
      V = if (K == 2) {
        1L
      } else if (loss == "logit") {
        K - 1L
      } else {
        K
      },
      K = K,
      y = NULL,
      cls = cls,
      Y = NULL,
      centre = if (K == 2) {
        "none"
      } else if (loss == "logit") {
        "ref"
      } else {
        "mean"
      },
      n = length(cls)
    ),
    class = "rpf_crit"
  )
}

# What the out-of-bag choices of a fit score: the numeric outcome (regression, and a binary L1 / L2 classification
# forest, whose 0/1 column is scored as a regression response, by the same code), or the classification criterion
oob_score_y <- function(mode, loss, y) {
  if (mode == "regression") {
    return(as.numeric(y))
  }
  crit <- oob_criterion(mode, loss, y)
  if (crit_is_scalar_mse(crit)) crit$y else crit
}

# Whether the criterion is the regression one (one column, squared error)
crit_is_scalar_mse <- function(crit) {
  crit$V == 1L && crit$type != "logloss"
}

# The criterion restricted to the rows `idx`
crit_subset <- function(crit, idx) {
  crit$y <- if (!is.null(crit$y)) crit$y[idx]
  crit$cls <- if (!is.null(crit$cls)) crit$cls[idx]
  crit$Y <- if (!is.null(crit$Y)) crit$Y[idx, , drop = FALSE]
  crit$n <- length(idx)
  crit
}

# The classification level. The tree families are averaged on the score (logit) scale, and the mean of their
# probabilities lies above the probability of their mean score where the probability is convex in the score (a rare
# class: Jensen), more so the more the families disagree, i.e. with more splits. So a forest predicts a rare class too
# rarely. Its level is the constant c (one value per score column) that the forest's score gets after the fit (the C++
# shift_level(), in every family's null leaf, so part of the intercept): the maximum-likelihood intercept of the
# log-loss given the forest's out-of-bag scores, whose equations are the balance property
#   sum_i p_k(c + g_i) - Delta_ik = n_k   (every free class k: the mean probability is the class's frequency)
# over the rows `idx` with an out-of-bag family, g the OOB mean scores (G; the half means G0, G1, half counts m0, m1),
# B the families. Delta is the second-order term of the OOB correction (criterion.R, above): an OOB mean of m families
# scatters more than the forest's mean of B, which moves its expected probability by
#   Delta_k = (1/m - 1/B) / 2 * sum_uv d2 p_k / dg_u dg_v V_uv = (1/m - 1/B) / 2 * w p_k ((d_k - p'd)^2 - (p'd^2 - (p'd)^2))
# with V = w d d' the across-family covariance from the halves (d over all K classes, the reference class's 0); rows
# without both halves get none. Solved by Newton's method from c = 0 with the Jacobian of the first term,
# sum_i diag(p_i) - p_i p_i' (positive definite), halving the step while the residual's norm does not fall; it stops at
# a largest residual below 1e-10 of a row (usually after 3 to 5 steps). Free classes: binary the second class;
# multiclass logit the K - 1 classes after the first (the reference, score 0); multiclass exponential (K scores, the
# softmax) the classes after the first, the first's score then fixed at 0, and c is centred at the end (the softmax
# ignores what all classes share; SAMME's scores are centred). 0 where a free class or the reference class has no
# out-of-bag row (no finite level). Returns V values, with the attribute `iterations`.
#
# Its arithmetic is plain double arithmetic in a fixed order (level_colsums(), level_rowsums(), level_solve()), so the
# level's bits depend on its inputs and on exp() and log() only: no BLAS or LAPACK, whose crossprod() and solve() round
# differently by library (OpenBLAS's triangular solve multiplies by the pivot's reciprocal where the reference BLAS
# divides by it), and no sum(), colSums() or rowSums(), which add in long double (80 bits on x86_64, 64 on arm64 macOS).
# Then the stopping test and the step halving take the same branches on every machine too.
classif_level <- function(crit, idx, G, G0, G1, m0, m1, B) {
  G <- as.matrix(G)
  V <- crit$V
  n <- nrow(G)
  cls <- crit$cls[idx]
  K <- crit$K
  zero <- structure(numeric(V), iterations = 0L)
  if (crit$type != "logloss" || n == 0) {
    return(zero)
  }
  n_k <- tabulate(cls, K)
  if (any(n_k == 0)) {
    return(zero)
  }
  # the free classes and their score columns (exponential's softmax: the first score fixed)
  free <- seq_len(K)[-1]
  free_col <- if (crit$link == "softmax") free else free - 1L
  D <- as.matrix(G0) - as.matrix(G1)
  m <- m0 + m1
  both <- m0 > 0 & m1 > 0
  coef <- (1 / m - 1 / B) * 0.5 / (1 / pmax(m0, 1) + 1 / pmax(m1, 1))
  coef[!both] <- 0
  # the full probabilities (n x K) and the score differences of the halves over all K classes (reference 0)
  full <- function(M) {
    if (crit$link == "softmax") M else cbind(0, M)
  }
  Df <- full(D)
  probs <- function(cv) {
    Gc <- G + rep(cv, each = n)
    if (crit$link == "sigmoid") {
      p <- stats::plogis(Gc[, 1])
      return(cbind(1 - p, p))
    }
    Gf <- full(Gc)
    # log-sum-exp of the rows (crit_lse()), the sum in a fixed order
    mx <- apply(Gf, 1, max)
    exp(Gf - (mx + log(level_rowsums(exp(Gf - mx)))))
  }
  resid <- function(Pf) {
    pd <- level_rowsums(Pf * Df)
    vp <- level_rowsums(Pf * Df^2) - pd^2
    Delta <- coef * Pf * ((Df - pd)^2 - vp)
    n_k[free] - level_colsums(Pf[, free, drop = FALSE] - Delta[, free, drop = FALSE])
  }
  # the Jacobian sum_i diag(p_i) - p_i p_i' of the free classes' probabilities P: entry (u, v) sums p_u (1[u = v] - p_v)
  jacobian <- function(P) {
    nf <- ncol(P)
    uv <- which(upper.tri(diag(nf), diag = TRUE), arr.ind = TRUE)
    terms <- vapply(
      seq_len(nrow(uv)),
      function(j) {
        u <- uv[j, 1]
        v <- uv[j, 2]
        P[, u] * (as.numeric(u == v) - P[, v])
      },
      numeric(nrow(P))
    )
    s <- level_colsums(matrix(terms, nrow(P)))
    J <- matrix(0, nf, nf)
    J[uv] <- s
    J[uv[, 2:1, drop = FALSE]] <- s
    J
  }
  cv <- numeric(V)
  Pf <- probs(cv)
  r <- resid(Pf)
  it <- 0L
  while (max(abs(r)) > 1e-10 * n && it < 50L) {
    it <- it + 1L
    step <- level_solve(jacobian(Pf[, free, drop = FALSE]), r)
    if (is.null(step) || !all(is.finite(step))) {
      break
    }
    t <- 1
    improved <- FALSE
    for (h in seq_len(30)) {
      cv_new <- cv
      cv_new[free_col] <- cv[free_col] + t * step
      Pf_new <- probs(cv_new)
      r_new <- resid(Pf_new)
      if (all(is.finite(r_new)) && level_sum(r_new^2) < level_sum(r^2)) {
        improved <- TRUE
        break
      }
      t <- t / 2
    }
    if (!improved) {
      break
    }
    cv <- cv_new
    Pf <- Pf_new
    r <- r_new
  }
  if (!all(is.finite(cv))) {
    return(zero)
  }
  if (crit$link == "softmax") {
    cv <- cv - level_sum(cv) / length(cv)
  }
  structure(cv, iterations = it)
}

# The level's sums (classif_level()), in double arithmetic in a fixed order, the same on every machine. Column sums
# of an n x k matrix: pairwise, the first half of the rows plus the second half (the odd row out kept last) until one
# row is left.
level_colsums <- function(M) {
  M <- as.matrix(M)
  if (nrow(M) == 0L) {
    return(numeric(ncol(M)))
  }
  while (nrow(M) > 1L) {
    h <- nrow(M) %/% 2L
    top <- M[seq_len(h), , drop = FALSE] + M[h + seq_len(h), , drop = FALSE]
    M <- if (nrow(M) %% 2L == 1L) rbind(top, M[nrow(M), , drop = FALSE]) else top
  }
  M[1, ]
}

# The sum of a vector (level_colsums())
level_sum <- function(x) {
  level_colsums(matrix(as.numeric(x), ncol = 1L))
}

# Row sums of an n x k matrix: the columns added from the first
level_rowsums <- function(M) {
  s <- M[, 1]
  for (j in seq_len(ncol(M))[-1]) {
    s <- s + M[, j]
  }
  s
}

# The solution of J x = r for the level's Newton step (J the k x k Jacobian, positive definite): Gaussian elimination
# with the diagonal as pivots in their order, in double arithmetic. NULL if a pivot is not positive and finite (J
# numerically singular), as solve() fails.
level_solve <- function(J, r) {
  k <- length(r)
  A <- J
  b <- as.numeric(r)
  for (p in seq_len(k)) {
    if (!is.finite(A[p, p]) || A[p, p] <= 0) {
      return(NULL)
    }
    for (i in seq_len(k - p) + p) {
      l <- A[i, p] / A[p, p]
      for (j in seq_len(k - p) + p) {
        A[i, j] <- A[i, j] - l * A[p, j]
      }
      b[i] <- b[i] - l * b[p]
    }
  }
  x <- numeric(k)
  for (i in rev(seq_len(k))) {
    s <- b[i]
    for (j in seq_len(k - i) + i) {
      s <- s - A[i, j] * x[j]
    }
    x[i] <- s / A[i, i]
  }
  x
}

# Whether a fit with these params gets the classification level (classif_level()): rpf()'s argument `level` (TRUE
# by default) is stored in params$level; a forest fitted before the level existed has no params$level and none.
# Classification with loss = "logit" or "exponential".
level_on <- function(params, mode) {
  mode == "classification" && isTRUE(params$level) && isTRUE(params$loss %in% c("logit", "exponential"))
}

# The level of the fitted forest `fit` (the C++ object, before shift_level()) with the classification criterion
# `crit` of its training rows: from its OOB path's last checkpoint, at its last split (the fit tracks it,
# rpf_impl()); for a forest without out-of-bag rows (deterministic) the balance on its own predictions of the
# training rows X. B: the families.
classif_fit_level <- function(fit, crit, X, B, deterministic) {
  in_sample <- function() {
    G <- matrix(fit$predict_matrix(X, 0, 0L), ncol = crit$V)
    n <- nrow(G)
    classif_level(crit, seq_len(n), G, G, G, rep(1L, n), rep(0L, n), 1)
  }
  if (deterministic) {
    return(in_sample())
  }
  path <- fit$oob_path()
  m0 <- path$count[, 1]
  m1 <- path$count[, 2]
  ok <- m0 + m1 > 0
  if (!any(ok)) {
    return(in_sample())
  }
  V <- crit$V
  cols <- (length(path$checkpoints) - 1L) * V + seq_len(V)
  s0 <- path$sum0[ok, cols, drop = FALSE]
  s1 <- path$sum1[ok, cols, drop = FALSE]
  m0 <- m0[ok]
  m1 <- m1[ok]
  classif_level(crit, which(ok), (s0 + s1) / (m0 + m1), s0 / pmax(m0, 1), s1 / pmax(m1, 1), m0, m1, B)
}

# The level of a forest with fixed splits, after its fit (rpf_impl() tracked its OOB path at its last split):
# classif_fit_level(), then classif_level_finish(). Returns the level.
classif_level_fixed <- function(fit, params, crit, X) {
  level <- classif_fit_level(fit, crit, X, params$ntrees, isTRUE(params$deterministic))
  classif_level_finish(fit, params, level)
  level
}

# A forest after its fit with the level: shift_level(), then the purification that rpf(purify = TRUE) asked for (the
# C++ fit does not purify such a forest itself, as the level comes after it; rpf_impl())
classif_level_finish <- function(fit, params, level) {
  fit$shift_level(as.numeric(level))
  if (isTRUE(params$purify)) {
    fit$purify_threads(0L, 0L, 2L)
  }
  invisible(fit)
}

# Log-sum-exp of the rows of the score columns G (n x V), with the reference score 0 joined if `ref`
crit_lse <- function(G, ref) {
  mx <- apply(G, 1, max)
  if (ref) {
    mx <- pmax(mx, 0)
  }
  s <- rowSums(exp(G - mx))
  if (ref) {
    s <- s + exp(-mx)
  }
  mx + log(s)
}

# The probabilities of the free score columns (n x V) and the per-row log-loss of the classes `cls`
crit_logloss <- function(crit, G, cls) {
  if (crit$link == "sigmoid") {
    g <- G[, 1]
    return(list(
      P = matrix(stats::plogis(g), ncol = 1),
      e = ifelse(cls == 2L, -stats::plogis(g, log.p = TRUE), -stats::plogis(-g, log.p = TRUE))
    ))
  }
  ref <- crit$link == "softmax_ref"
  lse <- crit_lse(G, ref)
  P <- exp(G - lse)
  # the class's score: logit's first level has the score 0, its other levels column cls - 1
  gy <- if (ref) {
    ifelse(cls == 1L, 0, G[cbind(seq_along(cls), pmax(cls - 1L, 1L))])
  } else {
    G[cbind(seq_along(cls), cls)]
  }
  list(P = P, e = lse - gy)
}

# Per-row raw (e) and corrected (corr) criterion of the scores G (n x V; the OOB means of the rows `idx` of the
# criterion's rows), the half means G0, G1, the half counts m0, m1 and B families. impute: rows without both halves get
# the mean second-order term of the others (the pruning), otherwise none (the path and the stop, as for regression).
crit_rows <- function(crit, idx, G, G0, G1, m0, m1, B, impute) {
  G <- as.matrix(G)
  D <- as.matrix(G0) - as.matrix(G1)
  m <- m0 + m1
  both <- m0 > 0 & m1 > 0
  w <- 1 / (1 / pmax(m0, 1) + 1 / pmax(m1, 1))
  if (crit$type != "logloss") {
    Y <- if (crit$V == 1L) matrix(crit$y[idx], ncol = 1) else crit$Y[idx, , drop = FALSE]
    e <- rowSums((Y - G)^2)
    q <- rowSums(D^2) * w
    q[!both] <- if (impute && any(both)) mean(q[both]) else 0
  } else {
    ll <- crit_logloss(crit, G, crit$cls[idx])
    e <- ll$e
    P <- ll$P
    q <- 0.5 * (rowSums(P * D^2) - rowSums(P * D)^2) * w
    if (impute && any(both) && any(!both)) {
      Dw <- D[both, , drop = FALSE] * sqrt(w[both])
      Vbar <- crossprod(Dw) / sum(both)
      Pn <- P[!both, , drop = FALSE]
      q[!both] <- 0.5 * (as.numeric(Pn %*% diag(Vbar)) - rowSums((Pn %*% Vbar) * Pn))
    } else {
      q[!both] <- 0
    }
  }
  list(e = e, corr = e - (1 / m - 1 / B) * q)
}

# The variance of the centred score columns of M (n x V; one component's or the prediction's), summed: the size a
# component's share is measured in. logit's K - 1 scores get the reference score 0 joined and exponential's K scores
# are centred across the classes (the softmax ignores what all classes share); L1 / L2 one-hot components already sum
# to 0 across the classes. The C++ coarsening centres the same way (RandomPlantedForest::score_centring()).
crit_centred_var <- function(crit, M) {
  M <- as.matrix(M)
  if (crit$centre == "ref") {
    M <- cbind(0, M)
  }
  if (crit$centre != "none") {
    M <- M - rowMeans(M)
  }
  var_pop <- function(v) mean((v - mean(v))^2)
  sum(apply(M, 2, var_pop))
}

# The columns of the components ki (indices) in a component-major n x (K * V) matrix
crit_blk <- function(ki, V) {
  as.integer(outer(seq_len(V), (ki - 1L) * V, "+"))
}

# Component-major columns F -> the sum over the components ki, per score column (n x V)
crit_sum_comp <- function(F, ki, V) {
  out <- matrix(0, nrow(F), V)
  for (v in seq_len(V)) {
    out[, v] <- rowSums(F[, (ki - 1L) * V + v, drop = FALSE])
  }
  out
}

# Rule R's coefficients under a classification criterion (other than one squared-error column, which keeps the
# regression least squares): one intercept per score column (exponential's softmax: the first fixed at 0) and one beta
# per kept component, shared by the classes. F: rows x (K * V) component-major; ki: the kept components; idx: the rows
# fitted on. Squared error: least squares stacked over (row, score column). Log-loss: maximum likelihood by Newton's
# method with step halving from beta = 1 and alpha = start_alpha (rule D's model); it stops at a relative change of
# the deviance below 1e-10, after 50 iterations, or where the deviance no longer falls (finite on separable data).
# Returns list(alpha = V values, beta = length(ki) values).
crit_fit_R <- function(crit, F, ki, idx, start_alpha = NULL) {
  V <- crit$V
  n <- length(idx)
  k <- length(ki)
  Fv <- lapply(seq_len(V), function(v) F[idx, (ki - 1L) * V + v, drop = FALSE])
  if (crit$type != "logloss") {
    Y <- if (V == 1L) matrix(crit$y[idx], ncol = 1) else crit$Y[idx, , drop = FALSE]
    X <- do.call(
      rbind,
      lapply(seq_len(V), function(v) {
        E <- matrix(0, n, V)
        E[, v] <- 1
        cbind(E, Fv[[v]])
      })
    )
    coef <- prune_lstsq_qr(X, as.numeric(Y))
    return(list(alpha = coef[seq_len(V)], beta = coef[-seq_len(V)]))
  }
  cls <- crit$cls[idx]
  free_a <- if (crit$link == "softmax") seq_len(V)[-1] else seq_len(V)
  na <- length(free_a)
  p <- na + k
  Xv <- lapply(seq_len(V), function(v) {
    E <- matrix(0, n, na)
    j <- match(v, free_a)
    if (!is.na(j)) {
      E[, j] <- 1
    }
    cbind(E, Fv[[v]])
  })
  eta <- function(th) {
    G <- matrix(0, n, V)
    for (v in seq_len(V)) {
      G[, v] <- as.numeric(Xv[[v]] %*% th)
    }
    G
  }
  dev <- function(th) sum(crit_logloss(crit, eta(th), cls)$e)
  # the indicator of the class among the free score columns (logit's first level: none)
  Yc <- matrix(0, n, V)
  col <- if (crit$link == "softmax_ref") {
    cls - 1L
  } else if (crit$link == "sigmoid") {
    ifelse(cls == 2L, 1L, 0L)
  } else {
    cls
  }
  Yc[cbind(which(col > 0), col[col > 0])] <- 1
  a0 <- if (is.null(start_alpha)) rep(0, V) else start_alpha
  if (crit$link == "softmax") {
    a0 <- a0 - a0[1]
  }
  th <- c(a0[free_a], rep(1, k))
  d_old <- dev(th)
  if (!is.finite(d_old)) {
    th <- c(rep(0, na), rep(1, k))
    d_old <- dev(th)
  }
  for (it in seq_len(50)) {
    P <- crit_logloss(crit, eta(th), cls)$P
    R <- P - Yc
    grad <- numeric(p)
    H <- matrix(0, p, p)
    U <- matrix(0, n, p)
    for (v in seq_len(V)) {
      grad <- grad + as.numeric(crossprod(Xv[[v]], R[, v]))
      H <- H + crossprod(Xv[[v]] * P[, v], Xv[[v]])
      U <- U + Xv[[v]] * P[, v]
    }
    H <- H - crossprod(U)
    step <- tryCatch(solve(H, grad), error = function(e) prune_lstsq(H, grad))
    if (!all(is.finite(step))) {
      break
    }
    t <- 1
    improved <- FALSE
    for (h in seq_len(30)) {
      th_new <- th - t * step
      d_new <- dev(th_new)
      if (is.finite(d_new) && d_new <= d_old) {
        improved <- TRUE
        break
      }
      t <- t / 2
    }
    if (!improved) {
      break
    }
    rel <- (d_old - d_new) / max(abs(d_old), 1e-300)
    th <- th_new
    d_old <- d_new
    if (rel < 1e-10) {
      break
    }
  }
  alpha <- numeric(V)
  alpha[free_a] <- th[seq_len(na)]
  list(alpha = alpha, beta = th[na + seq_len(k)])
}

# The scores alpha + sum_S beta_S f_S (rows x V) of rule R's coefficients `coef` on the component-major components F
crit_predict_R <- function(coef, F, ki, V) {
  G <- matrix(0, nrow(F), V)
  for (v in seq_len(V)) {
    G[, v] <- coef$alpha[v] + as.numeric(F[, (ki - 1L) * V + v, drop = FALSE] %*% coef$beta)
  }
  G
}

# The raw per-row criterion of the scores G (rows idx of the criterion's rows), as of a whole forest (no correction)
crit_loss_raw <- function(crit, idx, G) {
  G <- as.matrix(G)
  n <- nrow(G)
  crit_rows(crit, idx, G, G, G, rep(1, n), rep(1, n), 1, impute = FALSE)$e
}
