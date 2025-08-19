
#include "cpf.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <limits>
#include <unordered_set>

// ----------------- rpf subclass for classification -----------------

/**
 * \brief Create a prediction model based on Random Forests for classification data sets.
 */


void ClassificationRPF::L1_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_s[p]) - std::fabs((*split.Y)[individual][p]);
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_b[p]) - std::fabs((*split.Y)[individual][p]);
    }
  }
}

void ClassificationRPF::median_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = calcMedian(*split.Y, split.I_s);
  split.M_b = calcMedian(*split.Y, split.I_b);

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_s[p]) - std::fabs((*split.Y)[individual][p]);
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_b[p]) - std::fabs((*split.Y)[individual][p]);
    }
  }
}

void ClassificationRPF::logit_loss(Split &split)
{

  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);

  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;

  std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                { M = std::min(std::max(delta, M), 1 - delta); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                { M = std::min(std::max(delta, M), 1 - delta); });

  double M_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double M_bp = std::min(std::max(delta, split.M_bp), 1 - delta);

  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);

  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_sp) - W_s_mean[p]);
    }
    for (auto individual : split.I_b)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_bp) - W_b_mean[p]);
    }
  }

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));             // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));             // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
    }
  }

  for (auto individual : split.I_s)
  {
    split.min_sum += (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));         // ~ R_old
    split.min_sum -= (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
  }
  for (auto individual : split.I_b)
  {
    split.min_sum += (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));         // ~ R_old
    split.min_sum -= (1 - std::accumulate((*split.Y)[individual].begin(), (*split.Y)[individual].end(), 0.0)) * log(1 / (1 + std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
  }

  if (std::isnan(split.min_sum))
  {
    split.min_sum = INF;
  }
}

void ClassificationRPF::logit_loss_2(Split &split)
{

  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();

  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;

  std::vector<double> M_s2 = split.M_s;
  std::vector<double> M_b2 = split.M_b;

  std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                { M = std::max(delta, M); });

  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                { M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                { M = std::max(delta, 1 - M); });

  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);

  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_s2[p]) - W_s_mean[p]);
    }
    for (auto individual : split.I_b)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_b2[p]) - W_b_mean[p]);
    }
  }

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p]));         // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); // ~ R_new
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p]));         // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); // ~ R_new
    }
  }

  if (std::isnan(split.min_sum))
  {
    split.min_sum = INF;
  }
}

void ClassificationRPF::logit_loss_3(Split &split)
{

  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);

  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;

  std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                { M = std::max(delta, M); });

  std::for_each(M_s.begin(), M_s.end(), [&](double &M)
                { M = log(M); });
  std::for_each(M_b.begin(), M_b.end(), [&](double &M)
                { M = log(M); });

  double M_sp = std::max(delta, split.M_sp);
  double M_bp = std::max(delta, split.M_bp);

  M_sp = log(M_sp);
  M_bp = log(M_bp);

  double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
  double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);

  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);

  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;

  // std::vector<std::vector<double>> Y_s = split.Y_s;
  // std::vector<std::vector<double>> Y_b = split.Y_b;

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W_new[individual][p] = W_new[individual][p] + M_s[p] - sum_s - W_s_mean[p];
    }
    for (auto individual : split.I_b)
    {
      W_new[individual][p] = W_new[individual][p] + M_b[p] - sum_b - W_b_mean[p];
    }
  }

  std::vector<double> W_sp;
  std::vector<double> W_bp;
  std::vector<double> W_sp_new;
  std::vector<double> W_bp_new;

  std::vector<double> Y_sp;
  std::vector<double> Y_bp;

  for (auto individual : split.I_s)
  {
    W_sp.push_back(-accumulate(W[individual].begin(), W[individual].end(), 0.0));
    W_sp_new.push_back(-accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0));
    Y_sp.push_back(1 - accumulate(Y[individual].begin(), Y[individual].end(), 0.0));
  }

  for (auto individual : split.I_b)
  {
    W_bp.push_back(-accumulate(W[individual].begin(), W[individual].end(), 0.0));
    W_bp_new.push_back(-accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0));
    Y_bp.push_back(1 - accumulate(Y[individual].begin(), Y[individual].end(), 0.0));
  }

  /*
   W_s = transpose(W_s);
   W_s.push_back(W_sp);
   W_s = transpose(W_s);
   W_b = transpose(W_b);
   W_b.push_back(W_bp);
   W_b = transpose(W_b);
   W_s_new = transpose(W_s_new);
   W_s_new.push_back(W_sp_new);
   W_s_new = transpose(W_s_new);
   W_b_new = transpose(W_b_new);
   W_b_new.push_back(W_bp_new);
   W_b_new = transpose(W_b_new);
   Y_s=transpose(Y_s);
   Y_s.push_back(Y_sp);
   Y_s = transpose(Y_s);
   Y_b = transpose(Y_b);
   Y_b.push_back(Y_bp);
   Y_b = transpose(Y_b);
   */

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p]);
    }
    for (auto individual : split.I_b)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p]);
    }
  }

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));             // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (std::accumulate(W[individual].begin(), W[individual].end(), 0.0)));             // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (std::accumulate(W_new[individual].begin(), W_new[individual].end(), 0.0))); // ~ R_new
    }
  }

  if (std::isnan(split.min_sum))
  {
    split.min_sum = INF;
  }
}

void ClassificationRPF::logit_loss_4(Split &split)
{

  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();

  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;

  std::vector<double> M_s2 = split.M_s;
  std::vector<double> M_b2 = split.M_b;

  std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                { M = std::max(delta, M); });

  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                { M = std::max(delta, 1 - M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                { M = std::max(delta, 1 - M); });

  std::vector<double> W_s_mean = calcMean(*split.W, split.I_s);
  std::vector<double> W_b_mean = calcMean(*split.W, split.I_b);

  std::vector<std::vector<double>> W = *split.W, W_new = *split.W;

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_s[p] / M_s2[p]) - W_s_mean[p]);
    }
    for (auto individual : split.I_b)
    {
      W[individual][p] = exp(W[individual][p]);
      W_new[individual][p] = exp(W_new[individual][p] + log(M_b[p] / M_b2[p]) - W_b_mean[p]);
    }
  }

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p]));         // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); // ~ R_new
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.Y)[individual][p] * log(W[individual][p] / (1 + W[individual][p]));         // ~ R_old
      split.min_sum -= (*split.Y)[individual][p] * log(W_new[individual][p] / (1 + W_new[individual][p])); // ~ R_new
    }
  }

  if (std::isnan(split.min_sum))
  {
    split.min_sum = INF;
  }
}

void ClassificationRPF::exponential_loss(Split &split)
{

  split.min_sum = 0;
  split.M_s = std::vector<double>(value_size, 0);
  split.M_b = std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      W_s_sum[p] += (*split.W)[individual][p];
    }
    for (auto individual : split.I_b)
    {
      W_b_sum[p] += (*split.W)[individual][p];
    }
    for (auto individual : split.I_s)
    {
      sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    }
    for (auto individual : split.I_b)
    {
      sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    }

    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];

    sum_s[p] = std::min(std::max(delta, sum_s[p]), 1 - delta);
    sum_b[p] = std::min(std::max(delta, sum_b[p]), 1 - delta);
  }

  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);

  double sum_sp = std::min(std::max(delta, split.M_sp), 1 - delta);
  double sum_bp = std::min(std::max(delta, split.M_bp), 1 - delta);

  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_s[p] / sum_sp));
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_b[p] / sum_bp));
    }

    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }

  // check if valid result
  for (const auto &s : W_s_sum)
    if (s == 0)
      split.min_sum = INF;
  for (const auto &s : W_b_sum)
    if (s == 0)
      split.min_sum = INF;
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

void ClassificationRPF::exponential_loss_2(Split &split)
{

  split.min_sum = 0;
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);
  std::vector<double> sum_s2(value_size, 0);
  std::vector<double> sum_b2(value_size, 0);

  for (size_t p = 0; p < value_size; ++p)
  {

    for (auto individual : split.I_s)
    {
      W_s_sum[p] += (*split.W)[individual][p];
    }
    for (auto individual : split.I_b)
    {
      W_b_sum[p] += (*split.W)[individual][p];
    }

    for (auto individual : split.I_s)
    {
      sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    }
    for (auto individual : split.I_b)
    {
      sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    }

    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];

    sum_s2[p] = std::max(delta, 1 - sum_s[p]);
    sum_b2[p] = std::max(delta, 1 - sum_s[p]);

    sum_s[p] = std::max(delta, sum_s[p]);
    sum_b[p] = std::max(delta, sum_b[p]);
  }

  for (size_t p = 0; p < value_size; ++p)
  {

    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_s[p] / sum_s2[p]));
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * log(sum_b[p] / sum_b2[p]));
    }

    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }

  // check if valid result
  for (const auto &s : W_s_sum)
    if (s == 0)
      split.min_sum = INF;
  for (const auto &s : W_b_sum)
    if (s == 0)
      split.min_sum = INF;
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

void ClassificationRPF::exponential_loss_3(Split &split)
{

  split.min_sum = 0;
  split.M_s = std::vector<double>(value_size, 0);
  split.M_b = std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);

  for (size_t p = 0; p < value_size; ++p)
  {

    for (auto individual : split.I_s)
    {
      W_s_sum[p] += (*split.W)[individual][p];
    }
    for (auto individual : split.I_b)
    {
      W_b_sum[p] += (*split.W)[individual][p];
    }

    for (auto individual : split.I_s)
    {
      sum_s[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_s_sum[p]);
    }
    for (auto individual : split.I_b)
    {
      sum_b[p] += (((*split.Y)[individual][p] + 1) / 2) * ((*split.W)[individual][p] / W_b_sum[p]);
    }

    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];
    sum_s[p] = std::max(delta, sum_s[p]);
    sum_b[p] = std::max(delta, sum_b[p]);
    sum_s[p] = log(sum_s[p]);
    sum_b[p] = log(sum_b[p]);
  }

  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);

  double sum_sp = std::max(delta, split.M_sp);
  double sum_bp = std::max(delta, split.M_bp);

  sum_sp = log(sum_sp);
  sum_bp = log(sum_bp);

  sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
  sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);

  sum_sp = sum_sp / (sum_s.size() + 1);
  sum_bp = sum_bp / (sum_b.size() + 1);

  for (size_t p = 0; p < value_size; ++p)
  {

    for (auto individual : split.I_s)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * (sum_s[p] - sum_sp));
    }
    for (auto individual : split.I_b)
    {
      split.min_sum += (*split.W)[individual][p] * exp(-0.5 * (*split.Y)[individual][p] * (sum_b[p] - sum_bp));
    }

    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }

  // check if valid result
  for (const auto &s : W_s_sum)
    if (s == 0)
      split.min_sum = INF;
  for (const auto &s : W_b_sum)
    if (s == 0)
      split.min_sum = INF;
  if (std::isnan(split.min_sum))
    split.min_sum = INF;
}

// constructor with parameters split_try, t_try, purify_forest, deterministic, nthreads
ClassificationRPF::ClassificationRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                                     const String loss, const NumericVector parameters)
    : RandomPlantedForest( 
        samples_Y, 
        samples_X, 
        // pass first 13 parameters to base (includes split_structure)
        parameters.size() >= 13 ? parameters[Rcpp::Range(0, 12)] : parameters[Rcpp::Range(0, 11)] 
      )
{

  // Ensure correct Rcpp RNG state
  Rcpp::RNGScope scope;

  // initialize class members
  std::vector<double> pars = to_std_vec(parameters);
  if (loss == "L1")
  {
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  }
  else if (loss == "L2")
  {
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  }
  else if (loss == "median")
  {
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  }
  else if (loss == "logit")
  {
    this->loss = LossType::logit;
    this->calcLoss = &ClassificationRPF::logit_loss;
  }
  else if (loss == "logit_2")
  {
    this->loss = LossType::logit_2;
    this->calcLoss = &ClassificationRPF::logit_loss_2;
  }
  else if (loss == "logit_3")
  {
    this->loss = LossType::logit_3;
    this->calcLoss = &ClassificationRPF::logit_loss_3;
  }
  else if (loss == "logit_4")
  {
    this->loss = LossType::logit_4;
    this->calcLoss = &ClassificationRPF::logit_loss_4;
  }
  else if (loss == "exponential")
  {
    this->loss = LossType::exponential;
    this->calcLoss = &ClassificationRPF::exponential_loss;
  }
  else if (loss == "exponential_2")
  {
    this->loss = LossType::exponential_2;
    this->calcLoss = &ClassificationRPF::exponential_loss_2;
  }
  else if (loss == "exponential_3")
  {
    this->loss = LossType::exponential_3;
    this->calcLoss = &ClassificationRPF::exponential_loss_3;
  }
  else
  {
    Rcout << "Unkown loss function, set to default (L2)." << std::endl;
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  }
  if (pars.size() != 15)
  {
    Rcout << "Wrong number of parameters - set to default." << std::endl;
    this->max_interaction = 1;
    this->n_trees = 50;
    this->n_splits = 30;
    this->split_try = 10;
    this->t_try = 0.4;
    this->purify_forest = 0;
    this->deterministic = 0;
    this->nthreads = 1;
    this->cross_validate = 0;
    this->split_decay_rate_ = 0.1;
    this->max_candidates_   = 50;
    this->delete_leaves   = 1;
    this->delta = 0.1;
    this->epsilon = 0;
  }
  else
  {
    this->max_interaction = pars[0];
    this->n_trees = pars[1];
    this->n_splits = pars[2];
    this->split_try = pars[3];
    this->t_try = pars[4];
    this->purify_forest = pars[5];
    this->deterministic = pars[6];
    this->nthreads = pars[7];
    this->cross_validate = pars[8];
    this->split_decay_rate_ = pars[9];
    this->max_candidates_   = static_cast<size_t>(pars[10]);
    this->delete_leaves   =  pars[11];
    // pars[12] is split_structure for base; already consumed by base
    this->delta = pars[13];
    this->epsilon = pars[14];   
  }

  // set data and data related members
  this->set_data(samples_Y, samples_X);
}

// Mode 1: cur_trees_2 (classification variant)
Split ClassificationRPF::calcOptimalSplit_curTrees2(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{

  Split curr_split, min_split;
  min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y;
  curr_split.W = &weights;
  std::set<int> tree_dims;
  std::vector<double> unique_samples;
  int k;
  unsigned int n = 0;
  double leaf_size;

  // sample possible splits
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));  
  std::vector<size_t> split_candidates;

      // 1) Build weights = exp(-decay_rate * age)
    std::vector<double> weights_vec(possible_splits.size());
    for (size_t i = 0; i < possible_splits.size(); ++i) {
      weights_vec[i] = std::exp(-split_decay_rate_ * possible_splits[i].age);
    }

    // 2) Sample n_candidates indices *without* replacement
    std::vector<size_t> sample_idxs;
    sample_idxs.reserve(n_candidates);

    if (!deterministic) {
      // fully-qualified random_device:
      std::mt19937 gen{ std::random_device{}() };
      std::discrete_distribution<size_t> dist(weights_vec.begin(), weights_vec.end());
      std::vector<bool> used(possible_splits.size(), false);

      // draw until we have n_candidates distinct picks
      while (sample_idxs.size() < n_candidates) {
        size_t idx = dist(gen);
        if (!used[idx]) {
          used[idx] = true;
          sample_idxs.push_back(idx);
        }
      }
    } else {
      // deterministic fallback: first n_candidates
      for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i)
        sample_idxs.push_back(i);
    }

    split_candidates = sample_idxs;

  // track which one gave us the best split
  size_t chosen_idx = std::numeric_limits<size_t>::max();

  // consider a fraction of possible splits
  while (n < n_candidates)
  {

    // since size of possible splits changes, check if candidate in range
    if (possible_splits.empty())
      break;
    if (split_candidates[n] >= 0 && (size_t)split_candidates[n] >= possible_splits.size())
      continue;

    auto candidate = possible_splits.begin();
    std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
    k = candidate->dim - 1;                     // split dim of  candidate, converted to index starting at 0
    leaf_size = n_leaves[k];

    // Test if splitting in the  tree w.r.t. the coordinate "k" is an element of candidate tree
    tree_dims = candidate->tree->split_dims;
    tree_dims.erase(k + 1);
    tree_dims.erase(0);

    std::vector<std::shared_ptr<DecisionTree>> curr_trees;
    if (tree_dims.size() == 0)
      curr_trees.push_back(curr_family[std::set<int>{0}]);
    if (curr_family.find(tree_dims) != curr_family.end())
      curr_trees.push_back(curr_family[tree_dims]);
    if (curr_family.find(candidate->tree->split_dims) != curr_family.end())

      // go through all trees in current family
      for (auto &curr_tree : curr_trees)
      {

        // skip if tree has no leaves
        if (curr_tree->leaves.size() == 0)
          continue;

        // go through all leaves of current tree
     /*    for (auto &leaf : curr_tree->leaves)
        {

          std::vector<double> tot_sum(value_size, 0);

          // extract sample points according to individuals from X and Y
          unique_samples = std::vector<double>(leaf.individuals.size());
          for (unsigned int i = 0; i < leaf.individuals.size(); ++i)
          {
            unique_samples[i] = X[leaf.individuals[i]][k];
          }
          std::sort(unique_samples.begin(), unique_samples.end());
          unique_samples.erase(std::unique(unique_samples.begin(), unique_samples.end()), unique_samples.end());

          // check if number of sample points is within limit
          if (unique_samples.size() < 2 * leaf_size)
            continue;

          // consider split_try-number of samples
          std::vector<int> samples;
          if (deterministic)
          { // sequential samples if deterministic
            samples = std::vector<int>(std::min((int)unique_samples.size() - 1, 9));
            std::iota(samples.begin(), samples.end(), 1);
          }
          else
          { // randomly picked samples otherwise
            samples = std::vector<int>(split_try);
            for (size_t i = 0; i < samples.size(); ++i)
              samples[i] = R::runif(leaf_size, unique_samples.size() - leaf_size);
            std::sort(samples.begin(), samples.end());
          }

          // go through samples
          for (size_t sample_pos = 0; sample_pos < samples.size(); ++sample_pos)
          {

            // get samplepoint
            sample_point = unique_samples[samples[sample_pos]];

            // clear current split
            {
              curr_split.I_s.clear();
              curr_split.I_b.clear();
              curr_split.I_s.reserve(leaf.individuals.size());
              curr_split.I_b.reserve(leaf.individuals.size());
              curr_split.M_s = std::vector<double>(value_size, 0);
              curr_split.M_b = std::vector<double>(value_size, 0);
            }

            // get samples greater/smaller than samplepoint
            if (sample_pos == 0)
            {
              curr_split.sum_s = std::vector<double>(value_size, 0);
              curr_split.sum_b = std::vector<double>(value_size, 0);

              for (int individual : leaf.individuals)
              {
                if (X[individual][k] < sample_point)
                {
                  curr_split.I_s.push_back(individual);
                  curr_split.sum_s += Y[individual];
                }
                else
                {
                  curr_split.I_b.push_back(individual);
                  curr_split.sum_b += Y[individual];
                }
              }

              tot_sum = curr_split.sum_s + curr_split.sum_b;
            }
            else
            {

              for (int individual : leaf.individuals)
              {
                if (X[individual][k] < sample_point)
                {
                  if (X[individual][k] >= unique_samples[samples[sample_pos - 1]])
                  {
                    curr_split.sum_s += Y[individual];
                  }
                  curr_split.I_s.push_back(individual);
                }
                else
                {
                  curr_split.I_b.push_back(individual);
                }
              }

              curr_split.sum_b = tot_sum - curr_split.sum_s;
            }

            // accumulate squared mean and get mean
            (this->*ClassificationRPF::calcLoss)(curr_split);

            // update split if squared sum is smaller
            if (curr_split.min_sum < min_split.min_sum)
            {
              min_split = curr_split;
              min_split.tree_index = curr_tree;
              min_split.leaf_index = &leaf;
              min_split.split_coordinate = k + 1;
              min_split.split_point = sample_point;
              chosen_idx = split_candidates[n];
            }
          }
        } */

        // new: exactly split_try random (leaf, cut‐point) trials per tree
        for (size_t trial = 0; trial < split_try; ++trial) {
        // 1) pick a random leaf
        int leaf_idx = static_cast<int>(R::runif(0, curr_tree->leaves.size()));
        auto &leaf = curr_tree->leaves[leaf_idx];

        // 2) collect unique feature values in that leaf
        std::vector<double> unique_samples(leaf.individuals.size());
        for (size_t i = 0; i < leaf.individuals.size(); ++i)
            unique_samples[i] = X[leaf.individuals[i]][k];
        std::sort(unique_samples.begin(), unique_samples.end());
        unique_samples.erase(
            std::unique(unique_samples.begin(), unique_samples.end()),
            unique_samples.end()
        );

        // if too few distinct values, retry this trial
        if (unique_samples.size() < 2 * leaf_size) {
            --trial;
            continue;
        }

        // 3) pick one random cut‐point in [leaf_size, unique_samples.size()-leaf_size)
        int s_idx = static_cast<int>(
            R::runif(leaf_size, unique_samples.size() - leaf_size)
        );
        double sample_point = unique_samples[s_idx];

        // 4) partition individuals and accumulate sums
        curr_split.I_s.clear();
        curr_split.I_b.clear();
        curr_split.I_s.reserve(leaf.individuals.size());
        curr_split.I_b.reserve(leaf.individuals.size());
        curr_split.sum_s.assign(value_size, 0);
        curr_split.sum_b.assign(value_size, 0);

        for (int ind : leaf.individuals) {
            if (X[ind][k] < sample_point) {
                curr_split.I_s.push_back(ind);
                curr_split.sum_s += Y[ind];
            } else {
                curr_split.I_b.push_back(ind);
                curr_split.sum_b += Y[ind];
            }
        }

        // 5) compute your L1/L2/logit/etc. loss
        (this->*ClassificationRPF::calcLoss)(curr_split);

        // 6) update best‐so‐far split
        if (curr_split.min_sum < min_split.min_sum) {
            min_split = curr_split;
            min_split.tree_index       = curr_tree;
            min_split.leaf_index       = &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point      = sample_point;
            chosen_idx                 = split_candidates[n];
        }
    }

      }

    ++n;
  }

  for (size_t idx : split_candidates) {
    if (idx == chosen_idx) {
     possible_splits[idx].age = 0.0;        // reset for the winner
    } else {
     possible_splits[idx].age += 1.0;      // age everyone else
    }
  }

  return min_split;
}

// Mode 3: leaves (classification variant)
Split ClassificationRPF::calcOptimalSplit_leaves(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  Split curr_split, min_split; min_split.min_sum = std::numeric_limits<double>::infinity();
  curr_split.Y = &Y; curr_split.W = &weights;
  unsigned int raw_candidates = static_cast<unsigned int>(std::ceil(t_try * possible_splits.size()));
  unsigned int upper = std::min<size_t>(max_candidates_, possible_splits.size());
  unsigned int n_candidates = std::max<unsigned int>(1u, std::min<unsigned int>(raw_candidates, upper));
  std::vector<double> weights_vec(possible_splits.size());
  for (size_t i = 0; i < possible_splits.size(); ++i) weights_vec[i] = std::exp(-split_decay_rate_ * possible_splits[i].age);
  std::vector<size_t> sample_idxs; sample_idxs.reserve(n_candidates);
  if (!deterministic) {
    std::mt19937 gen{ std::random_device{}() };
    std::discrete_distribution<size_t> dist(weights_vec.begin(), weights_vec.end());
    std::vector<bool> used(possible_splits.size(), false);
    while (sample_idxs.size() < n_candidates) { size_t idx = dist(gen); if (!used[idx]) { used[idx] = true; sample_idxs.push_back(idx);} }
  } else {
    for (size_t i = 0; i < n_candidates && i < possible_splits.size(); ++i) sample_idxs.push_back(i);
  }
  int best_idx = -1;
  for (size_t idx : sample_idxs) {
    auto it = possible_splits.begin(); std::advance(it, idx);
    int k = it->dim - 1; int leaf_size = n_leaves[k];
    auto treePtr = it->tree; if (treePtr->leaves.empty() || it->leaf_idx >= treePtr->leaves.size()) continue;
    Leaf* leafPtr = &treePtr->leaves[it->leaf_idx];
    std::vector<double> unique; unique.reserve(leafPtr->individuals.size());
    for (int ind : leafPtr->individuals) unique.push_back(X[ind][k]);
    std::sort(unique.begin(), unique.end()); unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    int left = (int)leaf_size; int right = (int)unique.size() - (int)leaf_size; if (right <= left) continue;
    size_t window = (size_t)(right - left); size_t draws = std::min((size_t)split_try, window);
    std::unordered_set<int> used_pos;
    for (size_t t=0; t<draws; ++t) {
      int s_idx = -1;
      if (deterministic) {
        int range = right - left; int guess = left + (int)(((double)used_pos.size()+0.5)/((double)range+0.5)*range);
        if (guess < left) guess = left; if (guess >= right) guess = right - 1; int lo=guess, hi=guess;
        while (lo>=left || hi<right) { if (lo>=left && !used_pos.count(lo)) { s_idx=lo; break; } if (hi<right && !used_pos.count(hi)) { s_idx=hi; break; } --lo; ++hi; }
        if (s_idx == -1) for (int p=left;p<right;++p) if (!used_pos.count(p)) { s_idx=p; break; }
      } else { do { s_idx = (int)R::runif(left, right); } while (used_pos.count(s_idx)); }
      used_pos.insert(s_idx);
      double sp = unique[(size_t)s_idx];
      curr_split.I_s.clear(); curr_split.I_b.clear();
      curr_split.M_s.assign(value_size, 0); curr_split.M_b.assign(value_size, 0);
      curr_split.sum_s.assign(value_size, 0); curr_split.sum_b.assign(value_size, 0);
      for (int ind : leafPtr->individuals) { if (X[ind][k] < sp) { curr_split.I_s.push_back(ind); curr_split.sum_s += Y[ind]; } else { curr_split.I_b.push_back(ind); curr_split.sum_b += Y[ind]; } }
      (this->*ClassificationRPF::calcLoss)(curr_split);
      if (curr_split.min_sum < min_split.min_sum) { min_split = curr_split; min_split.tree_index = treePtr; min_split.leaf_index = leafPtr; min_split.split_coordinate = k + 1; min_split.split_point = sp; best_idx = (int)idx; }
    }
  }
  for (size_t idx : sample_idxs) { if ((int)idx != best_idx) possible_splits[idx].age += 1.0; else possible_splits[idx].age = 0.0; }
  return min_split;
}

// Mode 2: cur_trees_1 (classification variant)
Split ClassificationRPF::calcOptimalSplit_curTrees1(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  // reuse current implementation by sampling per-leaf candidates across predecessor/current trees
  // We delegate to the old flow by temporarily constructing the same sampling but using loss with W
  // For brevity, call the curTrees2 variant which already samples leaves within available trees
  return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family, weights);
}

// Mode 0: res_trees (classification variant)
Split ClassificationRPF::calcOptimalSplit_resTrees(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<RandomPlantedForest::ResultingTreeCandidate> &possible_trees, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  // Classification loss evaluation on res_trees follows the base structure; to keep changes minimal here,
  // we adopt the cur_trees_1 sampling over the trees in possible_trees' dims using our calcLoss and W.
  // Construct a transient SplitCandidate view equivalent and reuse curTrees1.
  std::vector<SplitCandidate> proxy;
  for (auto &c : possible_trees) {
    for (int k_dim : c.tree->split_dims) {
      proxy.emplace_back(k_dim, c.tree, (size_t)0);
    }
  }
  return this->calcOptimalSplit_curTrees1(Y, X, proxy, curr_family, weights);
}

// Dispatcher selecting by split_structure_mode_
Split ClassificationRPF::calcOptimalSplit(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                           std::vector<SplitCandidate> &possible_splits, TreeFamily &curr_family, std::vector<std::vector<double>> &weights)
{
  if (split_structure_mode_ == 3) return this->calcOptimalSplit_leaves(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 2) return this->calcOptimalSplit_curTrees1(Y, X, possible_splits, curr_family, weights);
  if (split_structure_mode_ == 1) return this->calcOptimalSplit_curTrees2(Y, X, possible_splits, curr_family, weights);
  return Split{};
}

void ClassificationRPF::create_tree_family(std::vector<Leaf> initial_leaves, size_t n)
{

  TreeFamily curr_family;
  curr_family.insert(std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))); // save tree with one leaf in the beginning

  // Seed per mode
  std::vector<SplitCandidate> possible_splits;
  std::vector<ResultingTreeCandidate> possible_trees;
  if (split_structure_mode_ == 0) {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_trees.emplace_back(treePtr);
    }
  } else if (split_structure_mode_ == 3) {
    auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
        std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
        if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
        if (!this->leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
      }
    };
    auto null_tree = curr_family[{0}];
    if (!null_tree->leaves.empty()) add_leaf_candidates(null_tree, 0);
  } else {
    for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
      auto treePtr = std::make_shared<DecisionTree>(DecisionTree({feature_dim}));
      curr_family.insert({{feature_dim}, treePtr});
      possible_splits.emplace_back(feature_dim, treePtr, (size_t)0);
    }
  }

  // sample data points with replacement
  int sample_index;
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;

  // deterministic
  if (deterministic)
  {
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  }
  else
  {
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);

    // bagging/subsampling
    for (size_t i = 0; i < sample_size; ++i)
    {
      sample_index = R::runif(0, sample_size - 1);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }

  // initialize weights
  std::vector<std::vector<double>> weights;
  switch (this->loss)
  {
  case LossType::logit:
  case LossType::logit_2:
  case LossType::logit_3:
  case LossType::logit_4:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 0);
    break;
  case LossType::exponential:
  case LossType::exponential_2:
  case LossType::exponential_3:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 1);
    break;
  default:
    weights = std::vector<std::vector<double>>(sample_size);
    for (auto &W : weights)
      W = std::vector<double>(value_size, 0);
  }

  // modify existing or add new trees through splitting
  Split curr_split;
  for (int split_count = 0; split_count < n_splits; ++split_count)
  {

    // find optimal split
    if (split_structure_mode_ == 0) curr_split = this->calcOptimalSplit_resTrees(samples_Y, samples_X, possible_trees, curr_family, weights);
    else curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family, weights);

    // continue only if we get a significant result
    if (!std::isinf(curr_split.min_sum))
    {

      // update pools by mode
      if (split_structure_mode_ == 0) {
        std::set<int> Dprime = curr_split.tree_index->split_dims; Dprime.insert(curr_split.split_coordinate); Dprime.erase(0);
        if (!this->resultingTreeExists(possible_trees, Dprime)) { if (auto found = treeExists(Dprime, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({Dprime, std::make_shared<DecisionTree>(DecisionTree(Dprime))}); possible_trees.emplace_back(curr_family[Dprime]); } }
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
          std::set<int> U = Dprime; U.insert(feature_dim); if (U.size() == Dprime.size()) continue; if (max_interaction >= 0 && U.size() > (size_t)max_interaction) continue; if (this->resultingTreeExists(possible_trees, U)) continue; if (auto found = treeExists(U, curr_family)) possible_trees.emplace_back(found); else { curr_family.insert({U, std::make_shared<DecisionTree>(DecisionTree(U))}); possible_trees.emplace_back(curr_family[U]); }
        }
      } else if (split_structure_mode_ == 3) {
        auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
            if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
            if (!this->leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
          }
        };
        // added after leaf construction below (we need indices)
      } else {
        if (curr_split.tree_index->split_dims.count(curr_split.split_coordinate) == 0) {
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
            std::set<int> curr_dims = curr_split.tree_index->split_dims; curr_dims.insert(curr_split.split_coordinate); curr_dims.insert(feature_dim); curr_dims.erase(0);
            if (possibleExists(feature_dim, possible_splits, curr_dims)) continue;
            if (max_interaction >= 0 && curr_dims.size() > (size_t)max_interaction) continue;
            if (auto found = treeExists(curr_dims, curr_family)) possible_splits.emplace_back(feature_dim, found, (size_t)0);
            else { curr_family.insert({curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))}); possible_splits.emplace_back(feature_dim, curr_family[curr_dims], (size_t)0); }
          }
        }
      }

      // update values of individuals of split interval
      std::vector<double> update_s = curr_split.M_s, update_b = curr_split.M_b;
      switch (this->loss)
      {
      case LossType::L1:
      case LossType::L2:
      case LossType::median:
      {
        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            samples_Y[individual] -= update_s;
          }
          else
          {
            samples_Y[individual] -= update_b;
          }
        }
        break;
      }
      case LossType::logit:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::min(std::max(epsilon, M), 1 - epsilon); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::min(std::max(epsilon, M), 1 - epsilon); });

        double M_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double M_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_sp) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_bp) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_2:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_3:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s.begin(), M_s.end(), [&](double &M)
                      { M = log(M); });
        std::for_each(M_b.begin(), M_b.end(), [&](double &M)
                      { M = log(M); });

        double M_sp = std::max(epsilon, curr_split.M_sp);
        double M_bp = std::max(epsilon, curr_split.M_bp);

        M_sp = log(M_sp);
        M_bp = log(M_bp);

        double sum_s = (std::accumulate(M_s.begin(), M_s.end(), 0.0) + M_sp) / (M_s.size() + 1);
        double sum_b = (std::accumulate(M_b.begin(), M_b.end(), 0.0) + M_bp) / (M_b.size() + 1);

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (unsigned int p = 0; p < M_s.size(); ++p)
        {
          update_s[p] = M_s[p] - sum_s - W_s_mean[p];
          update_b[p] = M_b[p] - sum_b - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::logit_4:
      {

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;

        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;

        std::for_each(M_s.begin(), M_s.end(), [this](double &M)
                      { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M)
                      { M = std::max(epsilon, M); });

        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M)
                      { M = std::max(epsilon, 1 - M); });

        std::vector<double> W_s_mean = calcMean(*curr_split.W, curr_split.I_s);
        std::vector<double> W_b_mean = calcMean(*curr_split.W, curr_split.I_b);

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
          {
            weights[individual] += update_s;
          }
          else
          {
            weights[individual] += update_b;
          }
        }

        break;
      }
      case LossType::exponential:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::min(std::max(epsilon, S), 1 - epsilon); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::min(std::max(epsilon, S), 1 - epsilon); });

        double sum_sp = std::min(std::max(epsilon, curr_split.M_sp), 1 - epsilon);
        double sum_bp = std::min(std::max(epsilon, curr_split.M_bp), 1 - epsilon);

        for (unsigned int p = 0; p < sum_s.size(); ++p)
        {
          update_s[p] = log(sum_s[p] / sum_sp);
          update_b[p] = log(sum_b[p] / sum_bp);
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (unsigned int p = 0; p < update_s.size(); ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_2:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;
        std::vector<double> sum_s2 = curr_split.M_s;
        std::vector<double> sum_b2 = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::max(epsilon, S); });

        std::for_each(sum_s2.begin(), sum_s2.end(), [this](double &S)
                      { S = std::max(epsilon, 1 - S); });
        std::for_each(sum_b2.begin(), sum_b2.end(), [this](double &S)
                      { S = std::max(epsilon, 1 - S); });

        for (size_t p = 0; p < value_size; ++p)
        {
          update_s[p] = log(sum_s[p] / sum_s2[p]);
          update_b[p] = log(sum_b[p] / sum_b2[p]);
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (size_t p = 0; p < value_size; ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      case LossType::exponential_3:
      {

        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;

        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S)
                      { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S)
                      { S = std::max(epsilon, S); });

        std::for_each(sum_s.begin(), sum_s.end(), [&](double &S)
                      { S = log(S); });
        std::for_each(sum_b.begin(), sum_b.end(), [&](double &S)
                      { S = log(S); });

        double sum_sp = std::max(epsilon, curr_split.M_sp);
        double sum_bp = std::max(epsilon, curr_split.M_bp);

        sum_sp = log(sum_sp);
        sum_bp = log(sum_bp);

        sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
        sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);

        sum_sp = sum_sp / (sum_s.size() + 1);
        sum_bp = sum_bp / (sum_b.size() + 1);

        for (size_t p = 0; p < sum_s.size(); ++p)
        {
          update_s[p] = sum_s[p] - sum_sp;
          update_b[p] = sum_b[p] - sum_bp;
        }

        for (int individual : curr_split.leaf_index->individuals)
        {
          for (size_t p = 0; p < update_s.size(); ++p)
          {
            if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }
            else
            {
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }

        break;
      }
      }

      // construct new leaves
      Leaf leaf_s, leaf_b;
      {
        leaf_s.individuals = curr_split.I_s;
        leaf_b.individuals = curr_split.I_b;

        leaf_s.value = update_s;
        leaf_b.value = update_b;

        // initialize interval with split interval
        leaf_s.intervals = curr_split.leaf_index->intervals;
        leaf_b.intervals = curr_split.leaf_index->intervals;

        // interval of leaf with smaller individuals has new upper bound in splitting dimension
        leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point;
        // interval of leaf with bigger individuals has new lower bound in splitting dimension
        leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
      }

      // construct split_dims of resulting tree when splitting in split_coordinate
      std::set<int> resulting_dims = curr_split.tree_index->split_dims;
      resulting_dims.insert(curr_split.split_coordinate);
      resulting_dims.erase(0);

      // check if resulting tree already exists in family
      std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);

      // determine which tree is modified
      if ((curr_split.tree_index->split_dims.count(curr_split.split_coordinate))&& delete_leaves)
      { // if split variable is already in tree to be split
        // change values
        {
          leaf_s.value += curr_split.leaf_index->value;
          leaf_b.value += curr_split.leaf_index->value;
        }
        *curr_split.leaf_index = leaf_b;                 // replace old interval
        curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
        if (split_structure_mode_ == 3) {
          size_t idx_b = (size_t)(curr_split.leaf_index - &curr_split.tree_index->leaves[0]);
          size_t idx_s = curr_split.tree_index->leaves.size() - 1;
          auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
            for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
              std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
              if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
              if (!this->leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
            }
          };
          add_leaf_candidates(curr_split.tree_index, idx_b);
          add_leaf_candidates(curr_split.tree_index, idx_s);
        }
      }
      else
      {                                       // otherwise
        found_tree->leaves.push_back(leaf_s); // append new leaves
        found_tree->leaves.push_back(leaf_b);
        if (split_structure_mode_ == 3) {
          size_t idx_s = found_tree->leaves.size() - 2; size_t idx_b = found_tree->leaves.size() - 1;
          auto add_leaf_candidates = [&](const std::shared_ptr<DecisionTree>& T, size_t li) {
            for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim) {
              std::set<int> res_dims = T->split_dims; res_dims.insert(feature_dim); res_dims.erase(0);
              if (max_interaction >= 0 && res_dims.size() > (size_t)max_interaction) continue;
              if (!this->leafCandidateExists(possible_splits, T, li, feature_dim)) possible_splits.emplace_back(feature_dim, T, li);
            }
          };
          add_leaf_candidates(found_tree, idx_s);
          add_leaf_candidates(found_tree, idx_b);
        }
      }
    }
  }

  // remove empty trees & clear individuals of each tree
  auto keys = getKeys(curr_family);
  for (auto &key : keys)
  {
    if (curr_family[key]->leaves.size() == 0)
    {
      curr_family.erase(key);
      continue;
    }
    for (auto &leaf : curr_family[key]->leaves)
    {
      leaf.individuals.clear();
    }
  }

  tree_families[n] = curr_family;
}

// fit forest to new data
void ClassificationRPF::fit()
{

  // setup initial set of individuals
  std::vector<int> initial_individuals(sample_size);
  std::iota(initial_individuals.begin(), initial_individuals.end(), 0);

  // initialize intervals with lower and upper bounds
  std::vector<Interval> initial_intervals(feature_size);
  for (int i = 0; i < feature_size; ++i)
    initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};

  // set properties of first leaf
  Leaf initial_leaf;
  {
    initial_leaf.value = std::vector<double>(value_size, 0);
    initial_leaf.individuals = initial_individuals;
    initial_leaf.intervals = initial_intervals;
  }
  std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf

  // initialize tree families
  this->tree_families = std::vector<TreeFamily>(n_trees);

  // Loop over number of tree families and dispatch threads in batches
  // of nhreads at once. Note: R RNG used in sampling is not thread-safe.
  unsigned int threads_to_use = static_cast<unsigned int>(nthreads);
  if (threads_to_use > 1 && !deterministic) {
    Rcout << "Non-deterministic mode with >1 threads is not supported (R RNG not thread-safe). Falling back to 1 thread.\n";
    threads_to_use = 1;
  }
  if (threads_to_use > 1)
  {
    if (threads_to_use > std::thread::hardware_concurrency())
    {
      Rcout << "Requested " << threads_to_use << " threads but only " << std::thread::hardware_concurrency() << " available" << std::endl;
    }
    // Create local thread count to not overwrite nthreads,
    // would get reported wrongly by get_parameters()
    unsigned int current_threads = threads_to_use;
    for (int n = 0; n < n_trees; n += current_threads)
    {
      if (n >= (n_trees - current_threads + 1))
      {
        current_threads = n_trees % current_threads;
      }

      std::vector<std::thread> threads(current_threads);
      for (int t = 0; t < current_threads; ++t)
      {
        // Rcout << "Dispatching thread " << (n + t + 1) << "/" << n_trees << std::endl;
        threads[t] = std::thread(&ClassificationRPF::create_tree_family, this, std::ref(initial_leaves), n + t);
      }
      for (auto &t : threads)
      {
        if (t.joinable())
          t.join();
      }
    }
  }
  else
  {
    for (int n = 0; n < n_trees; ++n)
    {
      create_tree_family(initial_leaves, n);
    }
  }

  // optionally purify tree
  if (purify_forest)
  {
    this->purify_3();
  }
  else
  {
    purified = false;
  }
}

/*  retrospectively change parameters of existing class object,
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data',
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(StringVector keys, NumericVector values)
{
  if (keys.size() != values.size())
  {
    Rcout << "Size of input vectors is not the same. " << std::endl;
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    if (keys[i] == "deterministic")
    {
      this->deterministic = values[i];
    }
    else if (keys[i] == "nthreads")
    {
      this->nthreads = values[i];
    }
    else if (keys[i] == "purify")
    {
      this->purify_forest = values[i];
    }
    else if (keys[i] == "n_trees")
    {
      this->n_trees = values[i];
    }
    else if (keys[i] == "n_splits")
    {
      this->n_splits = values[i];
    }
    else if (keys[i] == "t_try")
    {
      this->t_try = values[i];
    }
    else if (keys[i] == "split_try")
    {
      this->split_try = values[i];
    }
    else if (keys[i] == "max_interaction")
    {
      this->max_interaction = values[i];
    }
    else if (keys[i] == "cv")
    {
      this->cross_validate = values[i];
    }
    else if (keys[i] == "loss")
    {
      if (keys[i] == "L1")
      {
        this->loss = LossType::L1;
        this->calcLoss = &ClassificationRPF::L1_loss;
      }
      else if (keys[i] == "L2")
      {
        this->loss = LossType::L2;
        this->calcLoss = &ClassificationRPF::L2_loss;
      }
      else if (keys[i] == "median")
      {
        this->loss = LossType::median;
        this->calcLoss = &ClassificationRPF::median_loss;
      }
      else if (keys[i] == "logit")
      {
        this->loss = LossType::logit;
        this->calcLoss = &ClassificationRPF::logit_loss;
      }
      else if (keys[i] == "exponential")
      {
        this->loss = LossType::exponential;
        this->calcLoss = &ClassificationRPF::exponential_loss;
      }
      else
      {
        Rcout << "Unkown loss function." << std::endl;
      }
    }
    else if (keys[i] == "delta")
    {
      this->delta = values[i];
    }
    else if (keys[i] == "epsilon")
    {
      this->epsilon = values[i];
    }
    else if (keys[i] == "split_decay_rate") 
    {
      this->split_decay_rate_ = values[i];
    }
    else if (keys[i] == "max_candidates") 
    {
      this->max_candidates_ = static_cast<size_t>(values[i]);
    }
    else
    {
      Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}
