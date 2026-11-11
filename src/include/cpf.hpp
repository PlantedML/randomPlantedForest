#ifndef CPF_H
#define CPF_H

#include "rpf.hpp"

// The arithmetic of the classification losses "logit" and "exponential" on a region's sums (losses_logit.cpp,
// losses_exponential.cpp). A region R (the rows of a leaf, or of one side of a candidate split) moves every row's score
// by one step c. Its prior is a pseudo-rows (a / K1 per class, K1 classes) at the region's mean score m, so a region's
// step minimises the penalised loss: its rows' loss after the step plus the pseudo-rows' loss after it. The pseudo-rows
// sit at the level the region has now, so they pull towards "no change": a region at its optimum is not moved again
// (no drift), a pure region converges to the level log((n + a / K1) / (a / K1)) instead of growing without bound, and
// no clamp is needed. A split is scored by the penalised losses after the two sides' steps minus the parent's after its
// own step, so the parent's own step is not credited to a split, and every new region pays for its pseudo-rows: a pure
// leaf, or a repeated split, does not gain.
//
// Everything a region's criterion and step need is a sum of per-row terms (T doubles per row) and the row count n:
// - logit (K score columns, the reference class at score 0, y one-hot with the reference all 0, p = softmax(0, F)):
//   the terms y - p (K), the lower triangle by rows of diag(p) - pp' (K (K + 1) / 2), and the scores F (K). With m the
//   mean score and pi = softmax(0, m):
//     g = sum (y - p) + a (1 / K1 - pi),  H = sum (diag(p) - pp') + (a / 2)(I - 11' / K1)
//   (the pseudo-rows' Hessian at its upper bound, Boehning 1992; for K1 = 2 it is a / 4). The Newton step is
//   c = H^-1 g (Cholesky), its gain g' H^-1 g / 2. A region's value is -gain + L_ps(m), L_ps(f) = -(a / K1) sum_p f_p
//   + a log(1 + sum_p e^f_p) the pseudo-rows' loss at step 0; the rows' loss at step 0 cancels between the sides and
//   the parent. The step does not depend on which class is the reference.
// - exponential, binary (one score F, the log-odds, y = -1 / 1, the row's weight w = exp(-y F / 2)): the terms w for
//   y = 1, w for y = -1, and F. With V1 = W1 + (a / 2) e^(-m / 2) and V0 = W0 + (a / 2) e^(m / 2) the step is
//   c = log(V1 / V0) and the value, the penalised loss after it, 2 sqrt(V1 V0).
// - exponential, K >= 3 classes (SAMME: K centred scores g, the row of class c has the weight w = exp(-g_c)): the
//   terms w in its class's column (K) and the scores (K). With V_k = W_k + (a / K) e^(-m_k) the step is
//   d_k = log V_k - mean_j log V_j and the value K exp(mean_j log V_j).
// The exponential's step is the exact minimiser of the penalised loss. Logit's is one Newton step; with exact_step the
// two steps of the chosen split are the exact minimisers (Newton with step halving on the winner's rows).
struct DriftFree {
  enum Kind {
    logit,
    exp_binary,
    samme
  };
  Kind kind = logit;
  size_t K = 1; // score columns
  size_t T = 0; // terms per row
  double a = 1; // pseudo-rows per region
  // The buffers of value(), owned by its caller (one per thread)
  struct Scratch {
    std::vector<double> m, pi, G, H, L, z;
  };
  DriftFree() = default;
  DriftFree(Kind kind_, size_t K_, double a_);
  size_t h_size() const { return K * (K + 1) / 2; }
  // out[0 .. T) += a row's terms. logit: y (K), the scores F (K), the probabilities p (K); exponential: y (binary 1,
  // SAMME K), F, the weight w
  void add_logit(const double *y, const double *F, const double *p, double *out) const;
  void add_exp(const double *y, const double *F, double w, double *out) const;
  // The value of a region with the term sums t and n rows (see above), and its step if `step` is not null. scale: the
  // size of the region's loss (logit n, exponential the penalised loss at step 0), for the no-gain threshold. NaN if
  // the step is not defined.
  double value(const double *t, double n, double *step, double *scale, Scratch &sc) const;
  double value_logit(const double *t, double n, double *step, double *scale, Scratch &sc) const; // losses_logit.cpp
  double value_exp(const double *t, double n, double *step, double *scale, Scratch &sc) const; // losses_exponential.cpp
  // The criterion of a split into the sides s and b (term sums, rows), vp the parent region's value: the sides' values
  // minus the parent's, INF unless it is a gain (below -1e-12 times the loss's size)
  double score(const double *s, const double *b, double ns, double nb, double vp, Scratch &sc) const;
  // logit: the exact minimiser of the penalised loss of the region `rows` (scores F, outcomes Y), into c
  void logit_exact_step(
      const std::vector<int> &rows,
      const std::vector<std::vector<double>> &F,
      const std::vector<std::vector<double>> &Y,
      std::vector<double> &c
  ) const;
};

class ClassificationRPF : public RandomPlantedForest {

public:
  // exact_step: logit's leaves get the exact minimisers of the penalised loss (otherwise one Newton step); prior_rows:
  // the pseudo-rows a of the logit and exponential losses (cpf.hpp, DriftFree)
  ClassificationRPF(const std::string &loss, const RPFParams &params, bool exact_step = true, double prior_rows = 1);
  void fit() override;
  void set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values) override;
  ~ClassificationRPF() {};
  // The per-row out-of-bag criterion of the loss (rpf.hpp): L1 / L2 (and median) the squared error summed over the
  // score columns, H = 2I; logit and exponential the log-loss of the probabilities predict() makes, H = diag(p) - pp'
  // over the free score columns
  double
  oob_row_error(size_t i, const double *g, const double *d, double inv, const double *var, double coef) const override;
  // The centring of the score columns (rpf.hpp): multiclass logit 1, multiclass exponential 2, otherwise 0
  int score_centring() const override;
  // Internal, for tests: with set_parameters("keep_weights", 1) (which refits), each family keeps its sample rows'
  // state as it is when its growth ends: each sample row's row in the data (0-based), its scores (logit and
  // exponential) and its loss weight (exponential). Empty otherwise.
  struct KeptWeights {
    std::vector<int> rows;
    std::vector<std::vector<double>> weights, scores;
  };
  const std::vector<KeptWeights> &kept_weights() const { return kept_weights_; }
  // The loss's arithmetic (logit, exponential), for the test hook rpf_classif_check()
  const DriftFree &drift_free() const { return df_; }
  bool exact_step() const { return exact_step_; }

protected:
  // Maps a loss name to `loss`; throws on unknown names.
  void set_loss_function(const std::string &loss);

private:
  enum LossType {
    L1,
    L2,
    median,
    logit,
    exponential
  };
  LossType loss;
  bool exact_step_ = true;
  double prior_rows_ = 1;
  DriftFree df_;
  // Internal, for tests and measurements (set_parameters, refits): "scan" 0 scores every try of logit and
  // exponential by a pass over its leaf's rows instead of the running sums; "logit_terms" -1 (the default: kept for
  // up to keep_terms_max classes), 0 or 1 whether logit keeps each row's T terms (otherwise its probabilities)
  bool scan_ = true;
  int logit_terms_ = -1;
  // logit keeps each row's terms (T = 2K + K (K + 1) / 2 doubles: 3, 7, 18 and 63 for 2, 3, 5 and 10 classes) for up
  // to this many classes, otherwise its K probabilities, from which the terms are computed when the sums are built.
  // Kept terms made fits 14% (2 and 3 classes), 20% (5) and 19-24% (10) faster at 10,000 rows; above 10 classes they
  // would take over 0.5 KB per sample row and thread, growing with the square of the classes.
  static constexpr size_t keep_terms_max = 10;
  // Sets df_ for the loss and the data's value_size (fit)
  void set_drift_free();
  std::vector<double> null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows) override;
  // The loss of one tree family (FamilyLoss, rpf.hpp): logit and exponential its rows' scores and weights, the running
  // sums and the steps; L1 / L2 / median its working rows' Y and the criterion from the sides' rows (cpf.cpp)
  class Loss;
  std::unique_ptr<FamilyLoss>
  family_loss(std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx) override;
  void L1_loss(Split &split);
  void median_loss(Split &split);
  bool keep_weights_ = false;
  std::vector<KeptWeights> kept_weights_;
};

#endif
