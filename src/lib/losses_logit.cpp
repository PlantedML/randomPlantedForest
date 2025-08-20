// Classification losses: Logit family variants. Extracted from cpf.cpp.
#include "cpf.hpp"

void ClassificationRPF::logit_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M){ M = std::min(std::max(delta, M), 1 - delta); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M){ M = std::min(std::max(delta, M), 1 - delta); });
  double M_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double M_bp = std::min(std::max(delta, split.M_bp), 1 - delta);
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_sp) - W_s_mean[p]); }
    for (auto individual : split.I_b) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_bp) - W_b_mean[p]); }
  }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
    for (auto individual : split.I_b) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
  }
  for (auto individual : split.I_s) { split.min_sum += (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
  for (auto individual : split.I_b) { split.min_sum += (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::logit_loss_2(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::vector<double> M_s2 = split.M_s, M_b2 = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M){ M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M){ M = std::max(delta, 1 - M); });
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_s2[p]) - W_s_mean[p]); }
    for (auto individual : split.I_b) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_b2[p]) - W_b_mean[p]); }
  }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p])); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); }
    for (auto individual : split.I_b) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p])); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); }
  }
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::logit_loss_3(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_s.begin(), M_s.end(), [&](double &M){ M = log(M); });
  std::for_each(M_b.begin(), M_b.end(), [&](double &M){ M = log(M); });
  double M_sp = std::max(delta, split.M_sp);
  double M_bp = std::max(delta, split.M_bp);
  M_sp = log(M_sp);
  M_bp = log(M_bp);
  double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
  double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { W_new[individual][p] = W_new[individual][p] + M_s[p] - sum_s - W_s_mean[p]; }
    for (auto individual : split.I_b) { W_new[individual][p] = W_new[individual][p] + M_b[p] - sum_b - W_b_mean[p]; }
  }
  std::vector<double> W_sp, W_bp, W_sp_new, W_bp_new, Y_sp, Y_bp;
  for (auto individual : split.I_s) { W_sp.push_back(-accumulate(W[individual].begin(), W[individual].end(), 0.0)); W_sp_new.push_back(-accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0)); Y_sp.push_back(1 - accumulate(Y[individual].begin(), Y[individual].end(), 0.0)); }
  for (auto individual : split.I_b) { W_bp.push_back(-accumulate(W[individual].begin(), W[individual].end(), 0.0)); W_bp_new.push_back(-accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0)); Y_bp.push_back(1 - accumulate(Y[individual].begin(), Y[individual].end(), 0.0)); }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p]); }
    for (auto individual : split.I_b) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p]); }
  }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
    for (auto individual : split.I_b) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (std::accumulate(W[individual].begin(), W[individual].end(), 0.0))); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); }
  }
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::logit_loss_4(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  std::vector<double> M_s = split.M_s, M_b = split.M_b;
  std::vector<double> M_s2 = split.M_s, M_b2 = split.M_b;
  std::for_each(M_s.begin(), M_s.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M){ M = std::max(delta, M); });
  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M){ M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M){ M = std::max(delta, 1 - M); });
  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);
  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_s2[p]) - W_s_mean[p]); }
    for (auto individual : split.I_b) { W[individual][p] = exp(W[individual][p]); W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_b2[p]) - W_b_mean[p]); }
  }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p])); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); }
    for (auto individual : split.I_b) { split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p])); split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); }
  }
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}


