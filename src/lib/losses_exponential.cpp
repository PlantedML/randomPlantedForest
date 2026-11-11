// Classification losses: the exponential loss's part of DriftFree (cpf.hpp), binary and SAMME.
#include "cpf.hpp"
#include <cmath>
#include <limits>

// binary: w for y = 1, w for y = -1, F; SAMME: w in the column of the row's class (the one coded 1; the last if none
// is), the K scores
void DriftFree::add_exp(const double *y, const double *F, double w, double *out) const
{
  if (kind == exp_binary) {
    out[y[0] > 0 ? 0 : 1] += w;
    out[2] += F[0];
    return;
  }
  size_t c = 0;
  while (c + 1 < K && !(y[c] > 0))
    ++c;
  out[c] += w;
  for (size_t q = 0; q < K; ++q)
    out[K + q] += F[q];
}

// The penalised loss after the region's step (cpf.hpp): binary 2 sqrt(V1 V0), SAMME K exp(mean_k log V_k)
double DriftFree::value_exp(const double *t, double n, double *step, double *scale, Scratch &sc) const
{
  if (kind == exp_binary) {
    const double h = a / 2, m = t[2] / n;
    const double V1 = t[0] + h * std::exp(-0.5 * m), V0 = t[1] + h * std::exp(0.5 * m);
    if (scale)
      *scale = V1 + V0;
    const double l1 = std::log(V1), l0 = std::log(V0);
    if (step)
      step[0] = l1 - l0;
    return 2 * std::exp(0.5 * (l1 + l0));
  }
  std::vector<double> &lv = sc.m;
  lv.resize(K);
  const double ak = a / (double)K;
  double mean = 0, tot = 0;
  for (size_t k = 0; k < K; ++k) {
    const double V = t[k] + ak * std::exp(-t[K + k] / n);
    tot += V;
    lv[k] = std::log(V);
    mean += lv[k];
  }
  mean /= (double)K;
  if (scale)
    *scale = tot;
  if (step)
    for (size_t k = 0; k < K; ++k)
      step[k] = lv[k] - mean;
  return (double)K * std::exp(mean);
}
