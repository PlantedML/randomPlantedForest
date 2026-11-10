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
