// Classification losses: Logit family variants. Extracted from cpf.cpp.
//
// A logit loss reads the weights of the split leaf's individuals only. Each loss keeps exp(W) and exp(W_new) of those
// individuals in two row-major n x value_size buffers, n = |I_s| + |I_b|, with I_s's rows first. A candidate split
// costs O(leaf size), as for the other losses: copying the whole weight matrix (every individual of the sample) twice
// per candidate made the logit losses cost O(sample size) per candidate.
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

void ClassificationRPF::logit_loss(Split &split)
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
