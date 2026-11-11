#include "cpf.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include "internal_utils.hpp"
using namespace rpf_utils;
#include <limits>
#include <unordered_set>

// ----------------- rpf subclass for classification -----------------

/**
 * \brief Create a prediction model based on Random Forests for classification data sets.
 */

// Maps a loss-name string to `loss` and, for L1 / L2 / median, `calcLoss`
void ClassificationRPF::set_loss_function(const std::string &loss)
{
  calcLoss = nullptr;
  if (loss == "L1") {
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  } else if (loss == "L2") {
    this->loss = LossType::L2;
    this->calcLoss = &ClassificationRPF::L2_loss;
  } else if (loss == "median") {
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  } else if (loss == "logit") {
    this->loss = LossType::logit;
  } else if (loss == "exponential") {
    this->loss = LossType::exponential;
  } else {
    throw std::invalid_argument("Unknown loss function '" + loss + "'.");
  }
}

ClassificationRPF::ClassificationRPF(
    const std::string &loss, const RPFParams &params, bool exact_step, double prior_rows
)
    : RandomPlantedForest(params), exact_step_(exact_step), prior_rows_(prior_rows)
{
  set_loss_function(loss);
}

void ClassificationRPF::set_drift_free()
{
  const size_t K = value_size > 0 ? (size_t)value_size : 1;
  if (loss == LossType::logit)
    df_ = DriftFree(DriftFree::logit, K, prior_rows_);
  else if (loss == LossType::exponential)
    df_ = DriftFree(K == 1 ? DriftFree::exp_binary : DriftFree::samme, K, prior_rows_);
  else
    df_ = DriftFree();
}

// The fit of one leaf holding the rows of a family that made no split (fit_null_leaf()): the leaf update of
// create_tree_family() for that leaf, from the initial scores 0 (and weights 1). L1 and L2: the leaf's mean of the
// class indicators. logit: with exact_step the exact minimiser of the penalised loss, log((n_k + a / K1) / (n_ref +
// a / K1)), otherwise its Newton step; exponential its step (cpf.hpp, DriftFree), binary log((n_1 + a / 2) / (n_0 +
// a / 2)), SAMME log(n_k + a / K) - mean_j log(n_j + a / K).
std::vector<double>
ClassificationRPF::null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows)
{
  if (loss == LossType::median)
    return calcMedian(Y, rows);
  if (loss == LossType::L1 || loss == LossType::L2)
    return calcMean(Y, rows);
  const size_t K = value_size;
  std::vector<double> value(K, 0.0);
  if (loss == LossType::logit && exact_step_) {
    std::vector<double> n(K, 0.0);
    double n_ref = 0;
    for (int i : rows) {
      double s = 0;
      for (size_t p = 0; p < K; ++p) {
        n[p] += Y[i][p];
        s += Y[i][p];
      }
      n_ref += 1 - s;
    }
    const double ak = df_.a / (double)(K + 1);
    for (size_t p = 0; p < K; ++p)
      value[p] = std::log((n[p] + ak) / (n_ref + ak));
    return value;
  }
  std::vector<double> t(df_.T, 0.0), zero(K, 0.0), p0(K, 1.0 / (double)(K + 1));
  for (int i : rows) {
    if (loss == LossType::logit)
      df_.add_logit(Y[i].data(), zero.data(), p0.data(), t.data());
    else
      df_.add_exp(Y[i].data(), zero.data(), 1.0, t.data());
  }
  DriftFree::Scratch sc;
  df_.value(t.data(), (double)rows.size(), value.data(), nullptr, sc);
  for (double &v : value)
    if (!std::isfinite(v))
      v = 0;
  return value;
}

// ----------------- the loss of a tree family -----------------
// The split searches (splits_*.cpp) and create_tree_family() (rpf.cpp) are regression's: the same candidate pool, the
// same sampling of candidates, leaves and thresholds (and therefore the same random draws), the same aging. L1 / L2 /
// median score each try with calcLoss from the individuals on both sides (I_s / I_b, in the order of the leaf's
// individuals) and the sums of Y the search passes, and after a split update the working rows' Y; for L2, calcLoss is
// L2_loss (rpf.cpp), regression's expression in regression's order, so loss "L2" on a 0/1 response is the regression
// fit, bit for bit, ties included. logit and exponential (cpf.hpp, DriftFree) keep each row's scores (and the
// exponential its weight) and score a try from the sums of the rows' terms over its two sides: by the running sums
// the searches build (scan_terms() etc.), or with set_parameters("scan", 0) by a pass over the leaf's rows (score()).
// The winner's steps always come from a pass over its rows (finalize()), so the two give the same splits and leaf
// values but for exact ties between candidates (the sums are added in another order).

class ClassificationRPF::Loss : public FamilyLoss {
public:
  Loss(ClassificationRPF &forest, std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx)
      : f(forest), Y(samples_Y), df(forest.df_)
  {
    const size_t n = (size_t)f.sample_size, K = f.value_size;
    if (f.keep_weights_)
      kept_rows = boot_idx;
    drift = f.loss == LossType::logit || f.loss == LossType::exponential;
    if (!drift)
      return;
    T = df.T;
    scores.assign(n, std::vector<double>(K, 0.0));
    if (f.loss == LossType::exponential) {
      wt.assign(n, 1.0);
    } else {
      keep_terms = f.logit_terms_ < 0 ? K + 1 <= keep_terms_max : f.logit_terms_ != 0;
      if (keep_terms) {
        terms.assign(n * T, 0.0);
        std::vector<double> p0(K, 1.0 / (double)(K + 1));
        for (size_t i = 0; i < n; ++i)
          df.add_logit(Y[i].data(), scores[i].data(), p0.data(), &terms[i * T]);
      } else {
        probs.assign(n, std::vector<double>(K, 1.0 / (double)(K + 1)));
      }
    }
    side_s.resize(T);
    side_b.resize(T);
  }

  double score(Split &trial, const Leaf &leaf, int k, double sp, const std::vector<std::vector<double>> &X) override
  {
    if (drift) { // the pass: the sums of the terms over the two sides, in the leaf's row order
      std::fill(side_s.begin(), side_s.end(), 0.0);
      std::fill(side_b.begin(), side_b.end(), 0.0);
      size_t ns = 0, nb = 0;
      for (int ind : leaf.individuals) {
        if (X[ind][k] < sp) {
          add_terms(ind, side_s.data());
          ++ns;
        } else {
          add_terms(ind, side_b.data());
          ++nb;
        }
      }
      trial.min_sum = INF;
      if (ns == 0 || nb == 0)
        return INF;
      parent.resize(T);
      for (size_t q = 0; q < T; ++q)
        parent[q] = side_s[q] + side_b[q];
      const double vp = df.value(parent.data(), (double)(ns + nb), nullptr, nullptr, dsc);
      trial.min_sum = df.score(side_s.data(), side_b.data(), (double)ns, (double)nb, vp, dsc);
      return trial.min_sum;
    }
    trial.Y = &Y;
    trial.I_s.clear();
    trial.I_b.clear();
    for (int ind : leaf.individuals) {
      if (X[ind][k] < sp)
        trial.I_s.push_back(ind);
      else
        trial.I_b.push_back(ind);
    }
    (f.*f.calcLoss)(trial);
    return trial.min_sum;
  }

  // The winner's sides (I_s / I_b) and steps (M_s / M_b), by a pass over its rows; its score stays the search's
  void finalize(Split &winner, const std::vector<std::vector<double>> &X) override
  {
    if (std::isinf(winner.min_sum) || winner.leaf_index == nullptr)
      return;
    if (!drift) {
      score(winner, *winner.leaf_index, winner.split_coordinate - 1, winner.split_point, X);
      return;
    }
    const int k = winner.split_coordinate - 1;
    winner.I_s.clear();
    winner.I_b.clear();
    std::fill(side_s.begin(), side_s.end(), 0.0);
    std::fill(side_b.begin(), side_b.end(), 0.0);
    for (int ind : winner.leaf_index->individuals) {
      if (X[ind][k] < winner.split_point) {
        winner.I_s.push_back(ind);
        add_terms(ind, side_s.data());
      } else {
        winner.I_b.push_back(ind);
        add_terms(ind, side_b.data());
      }
    }
    const size_t K = f.value_size;
    winner.M_s.assign(K, 0.0);
    winner.M_b.assign(K, 0.0);
    if (f.loss == LossType::logit && f.exact_step_) {
      df.logit_exact_step(winner.I_s, scores, Y, winner.M_s);
      df.logit_exact_step(winner.I_b, scores, Y, winner.M_b);
    } else {
      df.value(side_s.data(), (double)winner.I_s.size(), winner.M_s.data(), nullptr, dsc);
      df.value(side_b.data(), (double)winner.I_b.size(), winner.M_b.data(), nullptr, dsc);
    }
    // a step that is not finite (not seen; a Hessian that is not positive definite) leaves the rows where they are
    for (std::vector<double> *M : {&winner.M_s, &winner.M_b})
      for (double &v : *M)
        if (!std::isfinite(v))
          v = 0;
  }

  void apply(
      const Split &split,
      const std::vector<std::vector<double>> &X,
      std::vector<std::vector<double>> &samples_Y,
      std::vector<double> &update_s,
      std::vector<double> &update_b
  ) override;

  void finish(size_t n) override
  {
    // (tests) the rows' state as the family's growth leaves it (cpf.hpp)
    if (!f.keep_weights_ || n >= f.kept_weights_.size())
      return;
    KeptWeights kw;
    kw.rows = std::move(kept_rows);
    kw.scores = std::move(scores);
    if (!wt.empty()) {
      kw.weights.resize(wt.size());
      for (size_t i = 0; i < wt.size(); ++i)
        kw.weights[i].assign(1, wt[i]);
    }
    f.kept_weights_[n] = std::move(kw);
  }

  // ---- the running sums (rpf.hpp) ----
  size_t scan_terms() const override { return drift && f.scan_ ? T : 0; }

  void add_terms(int ind, double *out) const override
  {
    if (keep_terms) {
      const double *t = &terms[(size_t)ind * T];
      for (size_t q = 0; q < T; ++q)
        out[q] += t[q];
    } else if (f.loss == LossType::logit) {
      df.add_logit(Y[ind].data(), scores[ind].data(), probs[ind].data(), out);
    } else {
      df.add_exp(Y[ind].data(), scores[ind].data(), wt[ind], out);
    }
  }

  void prefix_terms(const Leaf &leaf, const size_t *order, size_t m, double *out) const override
  {
    if (m == 0)
      return;
    std::fill(out, out + T, 0.0);
    add_terms(leaf.individuals[order[0]], out);
    for (size_t i = 1; i < m; ++i) {
      double *o = out + i * T;
      if (keep_terms) {
        const double *t = &terms[(size_t)leaf.individuals[order[i]] * T];
        for (size_t q = 0; q < T; ++q)
          o[q] = o[q - T] + t[q];
      } else {
        std::copy(o - T, o, o);
        add_terms(leaf.individuals[order[i]], o);
      }
    }
  }

  double region_value(const double *t, size_t n) const override
  {
    return df.value(t, (double)n, nullptr, nullptr, dsc);
  }

  double scan_score(const double *s, const double *b, size_t ns, size_t nb, double vp) const override
  {
    if (ns == 0 || nb == 0)
      return INF;
    return df.score(s, b, (double)ns, (double)nb, vp, dsc);
  }

private:
  ClassificationRPF &f;
  std::vector<std::vector<double>> &Y; // the family's working rows' Y
  const DriftFree &df;
  bool drift = false;      // logit or exponential
  size_t T = 0;            // their terms per row
  bool keep_terms = false; // logit: each row's terms are kept (terms), otherwise its probabilities (probs)
  std::vector<std::vector<double>> scores, probs;
  std::vector<double> wt, terms;
  std::vector<double> side_s, side_b, parent, refresh_p; // scratch
  mutable DriftFree::Scratch dsc;                        // scratch of the loss's arithmetic
  std::vector<int> kept_rows;                            // (keep_weights_, tests) each working row's row in the data
  // logit: row i's probabilities, or its terms, after its scores changed
  void refresh(size_t i)
  {
    const size_t K = f.value_size;
    double *pr;
    if (keep_terms) {
      refresh_p.resize(K);
      pr = refresh_p.data();
    } else {
      pr = probs[i].data();
    }
    const double *F = scores[i].data();
    double mx = 0;
    for (size_t q = 0; q < K; ++q)
      mx = std::max(mx, F[q]);
    double den = std::exp(-mx);
    for (size_t q = 0; q < K; ++q) {
      pr[q] = std::exp(F[q] - mx);
      den += pr[q];
    }
    for (size_t q = 0; q < K; ++q)
      pr[q] /= den;
    if (keep_terms) {
      double *t = &terms[i * T];
      std::fill(t, t + T, 0.0);
      df.add_logit(Y[i].data(), F, pr, t);
    }
  }
};

std::unique_ptr<FamilyLoss>
ClassificationRPF::family_loss(std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx)
{
  return std::unique_ptr<FamilyLoss>(new Loss(*this, samples_Y, boot_idx));
}

// The leaves' increments, the steps of the split's two sides (M_s, M_b), and the update of the split leaf's rows: L1 /
// L2 / median their Y; logit their scores (and probabilities or terms); exponential their scores and weights
void ClassificationRPF::Loss::apply(
    const Split &curr_split,
    const std::vector<std::vector<double>> &samples_X,
    std::vector<std::vector<double>> &samples_Y,
    std::vector<double> &update_s,
    std::vector<double> &update_b
)
{
  const size_t K = f.value_size;
  update_s = curr_split.M_s;
  update_b = curr_split.M_b;
  const int k = curr_split.split_coordinate - 1;
  const double sp = curr_split.split_point;
  switch (f.loss) {
  case LossType::L1:
  case LossType::L2:
  case LossType::median: {
    for (int individual : curr_split.leaf_index->individuals) {
      if (samples_X[individual][k] < sp) {
        samples_Y[individual] -= update_s;
      } else {
        samples_Y[individual] -= update_b;
      }
    }
    break;
  }
  case LossType::logit: {
    for (int individual : curr_split.leaf_index->individuals) {
      const std::vector<double> &c = samples_X[individual][k] < sp ? update_s : update_b;
      std::vector<double> &F = scores[(size_t)individual];
      for (size_t q = 0; q < K; ++q)
        F[q] += c[q];
      refresh((size_t)individual);
    }
    break;
  }
  case LossType::exponential: {
    // binary: the weight times exp(-y c / 2); SAMME: times exp(-d_c), c the row's class
    std::vector<double> fac_s(K + 1), fac_b(K + 1);
    if (K == 1) {
      fac_s[0] = std::exp(-0.5 * update_s[0]);
      fac_s[1] = std::exp(0.5 * update_s[0]);
      fac_b[0] = std::exp(-0.5 * update_b[0]);
      fac_b[1] = std::exp(0.5 * update_b[0]);
    } else {
      for (size_t q = 0; q < K; ++q) {
        fac_s[q] = std::exp(-update_s[q]);
        fac_b[q] = std::exp(-update_b[q]);
      }
    }
    for (int individual : curr_split.leaf_index->individuals) {
      const bool small = samples_X[individual][k] < sp;
      const std::vector<double> &c = small ? update_s : update_b;
      const std::vector<double> &fac = small ? fac_s : fac_b;
      const std::vector<double> &y = samples_Y[individual];
      size_t cls;
      if (K == 1) {
        cls = y[0] > 0 ? 0 : 1;
      } else {
        cls = 0;
        while (cls + 1 < K && !(y[cls] > 0))
          ++cls;
      }
      wt[(size_t)individual] *= fac[cls];
      std::vector<double> &F = scores[(size_t)individual];
      for (size_t q = 0; q < K; ++q)
        F[q] += c[q];
    }
    break;
  }
  }
}

// fit forest to new data
void ClassificationRPF::fit()
{
  set_drift_free();
  kept_weights_.assign(keep_weights_ ? (size_t)n_trees : 0, KeptWeights());
  // Use the base class multithreaded trainer with RNG seeding identical to regression
  RandomPlantedForest::fit();
}

/*  retrospectively change parameters of existing class object,
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data',
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values)
{
  if (keys.size() != values.size()) {
    warn("Size of input vectors is not the same. ");
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i) {
    if (keys[i] == "deterministic") {
      this->deterministic = values[i];
    } else if (keys[i] == "nthreads") {
      this->nthreads = values[i];
    } else if (keys[i] == "purify") {
      this->purify_forest = values[i];
    } else if (keys[i] == "n_trees") {
      this->n_trees = values[i];
    } else if (keys[i] == "n_splits") {
      this->n_splits = values[i];
    } else if (keys[i] == "t_try") {
      this->t_try = values[i];
    } else if (keys[i] == "split_try") {
      this->split_try = values[i];
    } else if (keys[i] == "max_interaction") {
      this->max_interaction = values[i];
    } else if (keys[i] == "cv") {
      this->cross_validate = values[i];
    } else if (keys[i] == "split_decay_rate") {
      this->split_decay_rate_ = values[i];
    } else if (keys[i] == "max_candidates") {
      this->max_candidates_ = static_cast<size_t>(values[i]);
    } else if (keys[i] == "exact_step") {
      this->exact_step_ = values[i] != 0;
    } else if (keys[i] == "prior_rows") {
      this->prior_rows_ = values[i];
    } else if (keys[i] == "keep_weights") { // internal, for tests (cpf.hpp)
      this->keep_weights_ = values[i] != 0;
    } else if (keys[i] == "scan") { // internal, for tests (cpf.hpp)
      this->scan_ = values[i] != 0;
    } else if (keys[i] == "logit_terms") { // internal, for measurements (cpf.hpp)
      this->logit_terms_ = values[i] < 0 ? -1 : (values[i] != 0 ? 1 : 0);
    } else {
      warn("Unkown parameter key  '" + keys[i] + "' .");
    }
  }
  this->fit();
}

// The per-row criterion of a classification forest (rpf.hpp; the same as R's crit_rows()). L1 / L2 (and median): the
// squared error summed over the score columns (Y holds their 0/1 coding), H = 2I. logit and exponential: the log-loss
// of the probabilities predict() makes: binary p = 1 / (1 + exp(-g)); multiclass logit the softmax of the K - 1 scores
// and the reference level's score 0; multiclass exponential the softmax of the K scores. H = diag(p) - pp' over the
// free score columns.
double ClassificationRPF::oob_row_error(
    size_t i, const double *g, const double *d, double inv, const double *var, double coef
) const
{
  const size_t V = value_size;
  const std::vector<double> &y = Y[i];
  const bool logit_like = loss == LossType::logit;
  const bool expo = loss == LossType::exponential;
  if (!logit_like && !expo) {
    double e = 0, q = 0;
    for (size_t v = 0; v < V; ++v) {
      const double r = y[v] - g[v];
      e += r * r;
      if (d)
        q += d[v] * d[v] / inv;
      else if (var)
        q += var[v];
    }
    return e - coef * q;
  }
  // the probabilities of the free score columns, and the log-loss of the row's class
  std::vector<double> p(V);
  double e = 0;
  if (V == 1) {
    // exponential codes the classes as -1 / 1, logit as 0 / 1
    const bool second = expo ? y[0] > 0 : y[0] > 0.5;
    p[0] = 1 / (1 + std::exp(-g[0]));
    // -log p(y): log(1 + exp(-g)) for the second class, log(1 + exp(g)) for the first, stably
    const double z = second ? -g[0] : g[0];
    e = z > 0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
  } else {
    // logit: the reference score 0 joins the softmax; exponential: the K = V scores alone
    double mx = logit_like ? 0.0 : -std::numeric_limits<double>::infinity();
    for (size_t v = 0; v < V; ++v)
      mx = std::max(mx, g[v]);
    double den = logit_like ? std::exp(-mx) : 0.0;
    for (size_t v = 0; v < V; ++v)
      den += std::exp(g[v] - mx);
    const double lse = mx + std::log(den);
    int cls = -1; // the column of the row's class (-1: the reference level of logit)
    for (size_t v = 0; v < V; ++v)
      if (y[v] > (expo ? 0.0 : 0.5))
        cls = (int)v;
    for (size_t v = 0; v < V; ++v)
      p[v] = std::exp(g[v] - lse);
    e = cls < 0 ? lse : lse - g[cls];
  }
  double q = 0;
  if (d) {
    double a = 0, b = 0;
    for (size_t v = 0; v < V; ++v) {
      a += p[v] * d[v] * d[v];
      b += p[v] * d[v];
    }
    q = 0.5 * (a - b * b) / inv;
  } else if (var)
    for (size_t v = 0; v < V; ++v)
      q += 0.5 * p[v] * (1 - p[v]) * var[v];
  return e - coef * q;
}

int ClassificationRPF::score_centring() const
{
  if (value_size < 2)
    return 0;
  if (loss == LossType::logit)
    return 1;
  if (loss == LossType::exponential)
    return 2;
  return 0;
}
