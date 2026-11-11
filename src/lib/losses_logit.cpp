// Classification losses: the drift-free losses' common arithmetic (DriftFree, cpf.hpp) and the logit loss's part.
#include "cpf.hpp"
#include <cmath>
#include <limits>

namespace {

// softmax(0, f) of K scores (the reference class at 0) into pi; returns log(1 + sum_p e^f_p)
double softmax_ref(const double *f, size_t K, double *pi)
{
  double mx = 0;
  for (size_t p = 0; p < K; ++p)
    mx = std::max(mx, f[p]);
  double den = std::exp(-mx);
  for (size_t p = 0; p < K; ++p) {
    pi[p] = std::exp(f[p] - mx);
    den += pi[p];
  }
  for (size_t p = 0; p < K; ++p)
    pi[p] /= den;
  return mx + std::log(den);
}

// The Cholesky factor L (row-major K x K, lower) of the symmetric matrix whose lower triangle by rows is h; false if
// it is not positive definite
bool cholesky(const double *h, size_t K, double *L)
{
  for (size_t p = 0; p < K; ++p) {
    double *Lp = L + p * K;
    for (size_t q = 0; q <= p; ++q) {
      double v = *h++;
      const double *Lq = L + q * K;
      for (size_t k = 0; k < q; ++k)
        v -= Lp[k] * Lq[k];
      if (q < p)
        Lp[q] = v / Lq[q];
      else if (!(v > 0))
        return false;
      else
        Lp[p] = std::sqrt(v);
    }
  }
  return true;
}

// x = (L L')^-1 g, with z = L^-1 g on the way; returns z'z
double cholesky_solve(const double *L, size_t K, const double *g, double *z, double *x)
{
  double zz = 0;
  for (size_t p = 0; p < K; ++p) {
    double v = g[p];
    for (size_t k = 0; k < p; ++k)
      v -= L[p * K + k] * z[k];
    z[p] = v / L[p * K + p];
    zz += z[p] * z[p];
  }
  if (x)
    for (size_t p = K; p-- > 0;) {
      double v = z[p];
      for (size_t k = p + 1; k < K; ++k)
        v -= L[k * K + p] * x[k];
      x[p] = v / L[p * K + p];
    }
  return zz;
}

} // namespace

DriftFree::DriftFree(Kind kind_, size_t K_, double a_) : kind(kind_), K(K_), a(a_)
{
  switch (kind) {
  case logit:
    T = 2 * K + h_size();
    break;
  case exp_binary:
    T = 3;
    break;
  case samme:
    T = 2 * K;
    break;
  }
}

double DriftFree::value(const double *t, double n, double *step, double *scale, Scratch &sc) const
{
  return kind == logit ? value_logit(t, n, step, scale, sc) : value_exp(t, n, step, scale, sc);
}

double DriftFree::score(const double *s, const double *b, double ns, double nb, double vp, Scratch &sc) const
{
  double scs = 0, scb = 0;
  const double v = value(s, ns, nullptr, &scs, sc) + value(b, nb, nullptr, &scb, sc) - vp;
  if (std::isfinite(v) && v < -1e-12 * (scs + scb))
    return v;
  return std::numeric_limits<double>::infinity();
}

void DriftFree::add_logit(const double *y, const double *F, const double *p, double *out) const
{
  double *h = out + K, *ms = out + K + h_size();
  for (size_t q = 0; q < K; ++q) {
    const double pq = p[q];
    out[q] += y[q] - pq;
    for (size_t r = 0; r < q; ++r)
      *h++ -= pq * p[r];
    *h++ += pq * (1 - pq);
    ms[q] += F[q];
  }
}

// -gain + L_ps(m) of one Newton step from the region's sums (cpf.hpp)
double DriftFree::value_logit(const double *t, double n, double *step, double *scale, Scratch &sc) const
{
  std::vector<double> &m = sc.m, &pi = sc.pi, &G = sc.G, &Hl = sc.H, &L = sc.L, &z = sc.z;
  const size_t Hn = h_size();
  m.resize(K);
  pi.resize(K);
  G.resize(K);
  Hl.resize(Hn);
  L.assign(K * K, 0.0);
  z.resize(K);
  if (scale)
    *scale = n;
  const double K1 = (double)(K + 1);
  double msum = 0;
  for (size_t p = 0; p < K; ++p) {
    m[p] = t[K + Hn + p] / n;
    msum += m[p];
  }
  const double lse = softmax_ref(m.data(), K, pi.data());
  for (size_t p = 0; p < K; ++p)
    G[p] = t[p] + a * (1.0 / K1 - pi[p]);
  const double *h = t + K;
  size_t j = 0;
  for (size_t p = 0; p < K; ++p)
    for (size_t q = 0; q <= p; ++q, ++j)
      Hl[j] = h[j] + (a / 2) * ((p == q ? 1.0 : 0.0) - 1.0 / K1);
  if (!cholesky(Hl.data(), K, L.data()))
    return std::numeric_limits<double>::quiet_NaN();
  const double gain = cholesky_solve(L.data(), K, G.data(), z.data(), step) / 2;
  const double pseudo = -a / K1 * msum + a * lse;
  return -gain + pseudo;
}

// The exact minimiser of the penalised loss sum_rows l(y_i, F_i + c) + L_ps(m + c): Newton from c = 0 with the true
// Hessian (the pseudo-rows' bound if it is not positive definite) and step halving until the loss falls by a quarter
// of the predicted decrease
void DriftFree::logit_exact_step(
    const std::vector<int> &rows,
    const std::vector<std::vector<double>> &F,
    const std::vector<std::vector<double>> &Y,
    std::vector<double> &c
) const
{
  const size_t Hn = h_size();
  const double K1 = (double)(K + 1);
  c.assign(K, 0.0);
  if (rows.empty())
    return;
  std::vector<double> m(K, 0.0), f(K), pi(K), grad(K), H(Hn), L(K * K), z(K), d(K), trial(K);
  for (int i : rows)
    for (size_t p = 0; p < K; ++p)
      m[p] += F[i][p];
  for (double &v : m)
    v /= (double)rows.size();
  // the penalised loss at step x
  auto objective = [&](const std::vector<double> &x) {
    double L = 0;
    for (int i : rows) {
      double yf = 0;
      for (size_t p = 0; p < K; ++p) {
        f[p] = F[i][p] + x[p];
        yf += Y[i][p] * f[p];
      }
      L += softmax_ref(f.data(), K, pi.data()) - yf;
    }
    double s = 0;
    for (size_t p = 0; p < K; ++p) {
      f[p] = m[p] + x[p];
      s += f[p];
    }
    return L - a / K1 * s + a * softmax_ref(f.data(), K, pi.data());
  };
  double obj = objective(c);
  for (int it = 0; it < 60; ++it) {
    // gradient and Hessian at c
    std::fill(grad.begin(), grad.end(), 0.0);
    std::fill(H.begin(), H.end(), 0.0);
    auto add = [&](const double *y, double w) {
      size_t j = 0;
      for (size_t p = 0; p < K; ++p) {
        grad[p] += w * (pi[p] - (y ? y[p] : 1.0 / K1));
        for (size_t q = 0; q <= p; ++q, ++j)
          H[j] += w * ((p == q ? pi[p] : 0.0) - pi[p] * pi[q]);
      }
    };
    for (int i : rows) {
      for (size_t p = 0; p < K; ++p)
        f[p] = F[i][p] + c[p];
      softmax_ref(f.data(), K, pi.data());
      add(Y[i].data(), 1.0);
    }
    for (size_t p = 0; p < K; ++p)
      f[p] = m[p] + c[p];
    softmax_ref(f.data(), K, pi.data());
    std::vector<double> H_rows = H;
    add(nullptr, a);
    if (!cholesky(H.data(), K, L.data())) {
      size_t j = 0;
      for (size_t p = 0; p < K; ++p)
        for (size_t q = 0; q <= p; ++q, ++j)
          H_rows[j] += (a / 2) * ((p == q ? 1.0 : 0.0) - 1.0 / K1);
      if (!cholesky(H_rows.data(), K, L.data()))
        break;
    }
    // the Newton decrement grad' H^-1 grad and the direction d = H^-1 grad
    const double dec = cholesky_solve(L.data(), K, grad.data(), z.data(), d.data());
    if (!std::isfinite(dec))
      break;
    // near the minimiser (the decrement below the objective's rounding scale times 1e6), full Newton steps, which
    // converge quadratically there; stop once the step is at the level of rounding
    if (dec <= 1e-10 * (1 + std::fabs(obj))) {
      double dmax = 0, cmax = 0;
      for (size_t p = 0; p < K; ++p) {
        c[p] -= d[p];
        dmax = std::max(dmax, std::fabs(d[p]));
        cmax = std::max(cmax, std::fabs(c[p]));
      }
      if (dmax <= 1e-14 * (1 + cmax))
        break;
      continue;
    }
    double s = 1, next = obj;
    for (int bt = 0; bt < 60; ++bt) {
      for (size_t p = 0; p < K; ++p)
        trial[p] = c[p] - s * d[p];
      next = objective(trial);
      if (next <= obj - 0.25 * s * dec)
        break;
      s *= 0.5;
    }
    if (!(next < obj))
      break;
    c = trial;
    obj = next;
  }
}
