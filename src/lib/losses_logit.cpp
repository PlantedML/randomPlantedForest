// Classification losses: Logit family variants. Extracted from cpf.cpp.
//
// A logit loss reads the weights of the split leaf's individuals only. Each loss keeps exp(W) and exp(W_new) of those
// individuals in two row-major n x value_size buffers, n = |I_s| + |I_b|, with I_s's rows first. A candidate split
// costs O(leaf size), as for the other losses: copying the whole weight matrix (every individual of the sample) twice
// per candidate made the logit losses cost O(sample size) per candidate.
//
// `logit` (the default classification loss) also reads the terms of the current weights, which do not depend on the
// candidate, from LogitTerms (cpf.hpp), and skips the log() of terms whose y is 0: per individual and candidate it
// takes K exponentials and one logarithm (one-hot Y) instead of 2K exponentials and 2K + 2 logarithms.
#include "cpf.hpp"

namespace {

// The split's individuals, I_s then I_b: the row order of the buffers.
std::vector<int> split_individuals(const Split &split)
{
  std::vector<int> individuals(split.I_s);
  individuals.insert(individuals.end(), split.I_b.begin(), split.I_b.end());
  return individuals;
}

// The sum of row j of a row-major matrix with `cols` columns, summed from the first column on.
double row_sum(const std::vector<double> &mat, size_t j, size_t cols)
{
  auto first = mat.begin() + j * cols;
  return std::accumulate(first, first + cols, 0.0);
}

} // namespace

thread_local LogitTerms *tls_logit_terms_ptr = nullptr;

void LogitTerms::reset(
    const std::vector<std::vector<double>> &W_, const std::vector<std::vector<double>> &Y_, size_t K_
)
{
  W = &W_;
  Y = &Y_;
  K = K_;
  const size_t n = W_.size();
  terms.assign(n * (K + 1), 0.0);
  e.assign(K, 0.0);
  y_class.assign(n, -1);
  irregular.assign(n, 0);
  n_irregular = 0;
  for (size_t i = 0; i < n; ++i) {
    // the one column (K: the reference class) where row i's y is not 0, or -1
    const std::vector<double> &y = Y_[i];
    const double y_ref = 1 - std::accumulate(y.begin(), y.end(), 0.0);
    int cls = y_ref != 0 ? (int)K : -1, nonzero = y_ref != 0;
    for (size_t p = 0; p < K; ++p)
      if (y[p] != 0) {
        cls = (int)p;
        ++nonzero;
      }
    y_class[i] = nonzero == 1 ? cls : -1;
    update_row(i);
  }
}

// The same expressions as logit_loss_direct()'s terms of the current weights, so the same values.
void LogitTerms::update_row(size_t i)
{
  const std::vector<double> &w = (*W)[i], &y = (*Y)[i];
  for (size_t p = 0; p < K; ++p)
    e[p] = exp(w[p]);
  const double den = 1 + std::accumulate(e.begin(), e.end(), 0.0);
  double *t = &terms[i * (K + 1)];
  for (size_t p = 0; p < K; ++p)
    t[p] = y[p] * log(e[p] / den);
  const double y_ref = 1 - std::accumulate(y.begin(), y.end(), 0.0);
  t[K] = y_ref * log(1 / den);
  // regular: one-hot, and every term outside the class is +-0 (not NaN)
  bool regular = y_class[i] >= 0;
  for (size_t p = 0; p <= K && regular; ++p)
    if ((int)p != y_class[i] && t[p] != 0)
      regular = false;
  const char irr = regular ? 0 : 1;
  n_irregular += (long)irr - (long)irregular[i];
  irregular[i] = irr;
}

// The logit loss of a candidate: the log-likelihood of the leaf's individuals under the current weights minus that
// under the weights the split would give. The terms of the current weights come from LogitTerms; the others are
// computed here.
//
// min_sum gets the same additions in the same order as in logit_loss_direct(), except additions of +-0: min_sum starts
// at +0, and a sum of doubles is -0 only if both operands are -0, so min_sum is never -0 and adding or subtracting +-0
// leaves it unchanged. A term y * log(r) with y = 0 and a finite r > 0 is +-0. So when every individual's y is one-hot,
// every cached term outside its class is +-0, and every r is finite and > 0 (checked here, else the general loop),
// only each individual's class column contributes: one log() per individual, in the original order (class 0's
// individuals in leaf order, then class 1's, ..., then the reference class's).
void ClassificationRPF::logit_loss(Split &split)
{
  LogitTerms *lt = tls_logit_terms_ptr;
  if (!lt || lt->W != split.W || lt->Y != split.Y || lt->K != value_size) {
    logit_loss_direct(split);
    return;
  }
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::min(std::max(delta, M), 1 - delta); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::min(std::max(delta, M), 1 - delta); });
  double M_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double M_bp = std::min(std::max(delta, split.M_bp), 1 - delta);
  const std::vector<std::vector<double>> &Wm = *split.W, &Ym = *split.Y;
  const std::vector<int> &I_s = split.I_s, &I_b = split.I_b;
  const size_t n_s = I_s.size(), n = n_s + I_b.size(), K = value_size, K1 = K + 1;
  const double *terms = lt->terms.data();
  std::vector<int> &inds = lt->inds;
  inds.resize(n);
  std::copy(I_s.begin(), I_s.end(), inds.begin());
  std::copy(I_b.begin(), I_b.end(), inds.begin() + n_s);
  std::vector<double> &W_new = lt->W_new, &W_new_den = lt->W_new_den;
  W_new.resize(n * K);
  W_new_den.resize(n);
  for (size_t p = 0; p < K; ++p) {
    // the mean weight of each side, as calcMean() sums and divides it (without its index checks)
    double W_s_mean = 0, W_b_mean = 0;
    for (size_t j = 0; j < n_s; ++j)
      W_s_mean += Wm[inds[j]][p];
    for (size_t j = n_s; j < n; ++j)
      W_b_mean += Wm[inds[j]][p];
    W_s_mean /= (int)n_s;
    W_b_mean /= (int)(n - n_s);
    const double c_s = log(M_s[p] / M_sp), c_b = log(M_b[p] / M_bp);
    for (size_t j = 0; j < n_s; ++j)
      W_new[j * K + p] = exp(Wm[inds[j]][p] + c_s - W_s_mean);
    for (size_t j = n_s; j < n; ++j)
      W_new[j * K + p] = exp(Wm[inds[j]][p] + c_b - W_b_mean);
  }
  // 1 + the sum of each individual's exponentials
  for (size_t j = 0; j < n; ++j)
    W_new_den[j] = 1 + row_sum(W_new, j, K);
  double min_sum = 0;

  // Every r = W_new / den and 1 / den is finite and > 0 if every W_new >= 1e-250 and den <= 1e50 (r >= 1e-300)
  // (a NaN fails both comparisons)
  bool fast = lt->n_irregular == 0;
  for (size_t q = 0; q < n * K; ++q)
    fast &= W_new[q] >= 1e-250;
  for (size_t j = 0; j < n; ++j)
    fast &= W_new_den[j] <= 1e50;

  if (fast) {
    // the leaf's positions grouped by class (counting sort: leaf order within a class)
    std::vector<size_t> &start = lt->class_start;
    std::vector<int> &by_class = lt->by_class;
    start.assign(K1 + 1, 0);
    by_class.resize(n);
    const int *y_class = lt->y_class.data();
    for (size_t j = 0; j < n; ++j)
      ++start[(size_t)y_class[inds[j]] + 1];
    for (size_t c = 0; c < K1; ++c)
      start[c + 1] += start[c];
    std::vector<size_t> &next = lt->class_next;
    next.assign(start.begin(), start.end() - 1);
    for (size_t j = 0; j < n; ++j)
      by_class[next[(size_t)y_class[inds[j]]]++] = (int)j;
    for (size_t p = 0; p < K; ++p) {
      for (size_t a = start[p]; a < start[p + 1]; ++a) {
        const size_t j = (size_t)by_class[a], ind = (size_t)inds[j];
        min_sum += terms[ind * K1 + p];
        min_sum -= Ym[ind][p] * log(W_new[j * K + p] / W_new_den[j]);
      }
    }
    for (size_t a = start[K]; a < start[K1]; ++a) {
      const size_t j = (size_t)by_class[a], ind = (size_t)inds[j];
      const std::vector<double> &Y_j = Ym[ind];
      const double y_ref = 1 - std::accumulate(Y_j.begin(), Y_j.end(), 0.0);
      min_sum += terms[ind * K1 + K];
      min_sum -= y_ref * log(1 / W_new_den[j]);
    }
  } else {
    // every term (the y = 0 ones too), as logit_loss_direct()
    for (size_t p = 0; p < K; ++p) {
      for (size_t j = 0; j < n; ++j) {
        const size_t ind = (size_t)inds[j];
        const double y = Ym[ind][p];
        min_sum += terms[ind * K1 + p];
        min_sum -= y * log(W_new[j * K + p] / W_new_den[j]);
      }
    }
    for (size_t j = 0; j < n; ++j) {
      const size_t ind = (size_t)inds[j];
      const std::vector<double> &Y_j = Ym[ind];
      const double y_ref = 1 - std::accumulate(Y_j.begin(), Y_j.end(), 0.0);
      min_sum += terms[ind * K1 + K];
      min_sum -= y_ref * log(1 / W_new_den[j]);
    }
  }
  split.min_sum = std::isnan(min_sum) ? INF : min_sum;
}

void ClassificationRPF::logit_loss_direct(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::min(std::max(delta, M), 1 - delta); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::min(std::max(delta, M), 1 - delta); });
  double M_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double M_bp = std::min(std::max(delta, split.M_bp), 1 - delta);
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  const std::vector<int> individuals = split_individuals(split);
  const size_t n = individuals.size(), n_s = split.I_s.size(), K = value_size;
  std::vector<double> W(n * K), W_new(n * K);
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n_s; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_s[p] / M_sp) - W_s_mean[p]);
    }
    for (size_t j = n_s; j < n; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_b[p] / M_bp) - W_b_mean[p]);
    }
  }
  // 1 + the sum of each individual's exponentials
  std::vector<double> W_den(n), W_new_den(n);
  for (size_t j = 0; j < n; ++j) {
    W_den[j] = 1 + row_sum(W, j, K);
    W_new_den[j] = 1 + row_sum(W_new, j, K);
  }
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n; ++j) {
      const double y = (*split.Y)[individuals[j]][p];
      split.min_sum += y * log(W[j * K + p] / W_den[j]);
      split.min_sum -= y * log(W_new[j * K + p] / W_new_den[j]);
    }
  }
  for (size_t j = 0; j < n; ++j) {
    const std::vector<double> &Y_j = (*split.Y)[individuals[j]];
    const double y_ref = 1 - std::accumulate(Y_j.begin(), Y_j.end(), 0.0);
    split.min_sum += y_ref * log(1 / W_den[j]);
    split.min_sum -= y_ref * log(1 / W_new_den[j]);
  }
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

void ClassificationRPF::logit_loss_2(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::vector<double> M_s2 = split.M_s, M_b2 = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(delta, 1 - M); });
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  const std::vector<int> individuals = split_individuals(split);
  const size_t n = individuals.size(), n_s = split.I_s.size(), K = value_size;
  std::vector<double> W(n * K), W_new(n * K);
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n_s; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_s[p] / M_s2[p]) - W_s_mean[p]);
    }
    for (size_t j = n_s; j < n; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_b[p] / M_b2[p]) - W_b_mean[p]);
    }
  }
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n; ++j) {
      const double y = (*split.Y)[individuals[j]][p];
      split.min_sum += y * log(W[j * K + p] / (1 + W[j * K + p]));
      split.min_sum -= y * log(W_new[j * K + p] / (1 + W_new[j * K + p]));
    }
  }
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

void ClassificationRPF::logit_loss_3(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_s.begin(), M_s.end(), [&](double &M) { M = log(M); });
  std::for_each(M_b.begin(), M_b.end(), [&](double &M) { M = log(M); });
  double M_sp = std::max(delta, split.M_sp);
  double M_bp = std::max(delta, split.M_bp);
  M_sp = log(M_sp);
  M_bp = log(M_bp);
  double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
  double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  const std::vector<int> individuals = split_individuals(split);
  const size_t n = individuals.size(), n_s = split.I_s.size(), K = value_size;
  std::vector<double> W(n * K), W_new(n * K);
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n_s; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + M_s[p] - sum_s - W_s_mean[p]);
    }
    for (size_t j = n_s; j < n; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + M_b[p] - sum_b - W_b_mean[p]);
    }
  }
  // the sum of each individual's exponentials
  std::vector<double> W_den(n), W_new_den(n);
  for (size_t j = 0; j < n; ++j) {
    W_den[j] = row_sum(W, j, K);
    W_new_den[j] = row_sum(W_new, j, K);
  }
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n; ++j) {
      const double y = (*split.Y)[individuals[j]][p];
      split.min_sum += y * log(W[j * K + p] / W_den[j]);
      split.min_sum -= y * log(W_new[j * K + p] / W_new_den[j]);
    }
  }
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

void ClassificationRPF::logit_loss_4(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::vector<double> M_s2 = split.M_s, M_b2 = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(delta, 1 - M); });
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  const std::vector<int> individuals = split_individuals(split);
  const size_t n = individuals.size(), n_s = split.I_s.size(), K = value_size;
  std::vector<double> W(n * K), W_new(n * K);
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n_s; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_s[p] / M_s2[p]) - W_s_mean[p]);
    }
    for (size_t j = n_s; j < n; ++j) {
      const double w = (*split.W)[individuals[j]][p];
      W[j * K + p] = exp(w);
      W_new[j * K + p] = exp(w + log(M_b[p] / M_b2[p]) - W_b_mean[p]);
    }
  }
  for (size_t p = 0; p < K; ++p) {
    for (size_t j = 0; j < n; ++j) {
      const double y = (*split.Y)[individuals[j]][p];
      split.min_sum += y * log(W[j * K + p] / (1 + W[j * K + p]));
      split.min_sum -= y * log(W_new[j * K + p] / (1 + W_new[j * K + p]));
    }
  }
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}
