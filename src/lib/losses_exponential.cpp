// Classification losses: Exponential family variants. Extracted from cpf.cpp.
#include "cpf.hpp"

void ClassificationRPF::exponential_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = std::vector<double>(value_size, 0);
  split.M_b = std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0), W_b_sum(value_size, 0), sum_s(value_size, 0), sum_b(value_size, 0);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) W_s_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_b) W_b_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_s) sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    for (auto individual : split.I_b) sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    split.M_s[p] = sum_s[p]; split.M_b[p] = sum_b[p];
    sum_s[p] = std::min(std::max(delta, sum_s[p]), 1 - delta);
    sum_b[p] = std::min(std::max(delta, sum_b[p]), 1 - delta);
  }
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  double sum_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double sum_bp = std::min(std::max(delta, split.M_bp), 1 - delta);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_s[p] / sum_sp));
    for (auto individual : split.I_b) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_b[p] / sum_bp));
    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }
  for (const auto &s : W_s_sum) if (s == 0) split.min_sum = INF;
  for (const auto &s : W_b_sum) if (s == 0) split.min_sum = INF;
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::exponential_loss_2(Split &split)
{
  split.min_sum = 0;
  std::vector<double> W_s_sum(value_size, 0), W_b_sum(value_size, 0), sum_s(value_size, 0), sum_b(value_size, 0), sum_s2(value_size, 0), sum_b2(value_size, 0);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) W_s_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_b) W_b_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_s) sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    for (auto individual : split.I_b) sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    split.M_s[p] = sum_s[p]; split.M_b[p] = sum_b[p];
    sum_s2[p] = std::max(delta, 1 - sum_s[p]); sum_b2[p] = std::max(delta, 1 - sum_b[p]);
    sum_s[p] = std::max(delta, sum_s[p]); sum_b[p] = std::max(delta, sum_b[p]);
  }
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_s[p] / sum_s2[p]));
    for (auto individual : split.I_b) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_b[p] / sum_b2[p]));
    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }
  for (const auto &s : W_s_sum) if (s == 0) split.min_sum = INF;
  for (const auto &s : W_b_sum) if (s == 0) split.min_sum = INF;
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::exponential_loss_3(Split &split)
{
  split.min_sum = 0;
  split.M_s = std::vector<double>(value_size, 0);
  split.M_b = std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0), W_b_sum(value_size, 0), sum_s(value_size, 0), sum_b(value_size, 0);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) W_s_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_b) W_b_sum[p] += (*split.W)[individual][p];
    for (auto individual : split.I_s) sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    for (auto individual : split.I_b) sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    split.M_s[p] = sum_s[p]; split.M_b[p] = sum_b[p];
    sum_s[p] = std::max(delta, sum_s[p]); sum_b[p] = std::max(delta, sum_b[p]);
    sum_s[p] = log(sum_s[p]); sum_b[p] = log(sum_b[p]);
  }
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  double sum_sp = std::max(delta, split.M_sp), sum_bp = std::max(delta, split.M_bp);
  sum_sp = log(sum_sp); sum_bp = log(sum_bp);
  sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
  sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);
  sum_sp = sum_sp / (sum_s.size() + 1); sum_bp = sum_bp / (sum_b.size() + 1);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * (sum_s[p] - sum_sp));
    for (auto individual : split.I_b) split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * (sum_b[p] - sum_bp));
    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }
  for (const auto &s : W_s_sum) if (s == 0) split.min_sum = INF;
  for (const auto &s : W_b_sum) if (s == 0) split.min_sum = INF;
  if (std::isnan(split.min_sum)) split.min_sum = INF;
}


