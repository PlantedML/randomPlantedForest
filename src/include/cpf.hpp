#ifndef CPF_H
#define CPF_H

#include "rpf.hpp"

// The logit loss's terms that do not depend on the candidate split: for every individual i of a tree family's sample,
// with e_p = exp(W[i][p]) and den = 1 + sum_p e_p, the K values Y[i][p] * log(e_p / den) and then
// (1 - sum_p Y[i][p]) * log(1 / den), row-major n x (K + 1) in `terms`. create_tree_family() fills them once and
// refreshes the split leaf's rows after each logit update, so a candidate reads them instead of recomputing them
// (losses_logit.cpp). It also holds one candidate's scratch.
struct LogitTerms {
  const std::vector<std::vector<double>> *W = nullptr;
  const std::vector<std::vector<double>> *Y = nullptr;
  size_t K = 0;
  std::vector<double> terms;
  // y_class[i]: the one column of row i whose y is not 0 (K: the reference class), -1 if there is none or several.
  // irregular[i]: row i is not one-hot, or a cached term outside its class is not +-0 (a NaN); n_irregular: their count
  std::vector<int> y_class;
  std::vector<char> irregular;
  long n_irregular = 0;
  // scratch of one candidate
  std::vector<double> e, W_new, W_new_den;
  std::vector<int> inds, by_class;
  std::vector<size_t> class_start, class_next;
  // Points at W and Y (K columns) and computes every row.
  void reset(const std::vector<std::vector<double>> &W_, const std::vector<std::vector<double>> &Y_, size_t K_);
  // Recomputes row i (after W[i] changed).
  void update_row(size_t i);
};
// The current tree family's LogitTerms on this thread, or nullptr. A pointer only: the object lives on
// create_tree_family()'s stack.
extern thread_local LogitTerms *tls_logit_terms_ptr;

class ClassificationRPF : public RandomPlantedForest {

public:
  using RandomPlantedForest::calcOptimalSplit;
  ClassificationRPF(const std::string &loss, const RPFParams &params, double delta = 0.1, double epsilon = 0);
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
  // Internal, for tests: with set_parameters("keep_weights", 1) (which refits), each family keeps the loss weights of
  // its sample rows as they are when its growth ends (one vector per sample row), and each sample row's row in the
  // data (0-based). Empty otherwise.
  struct KeptWeights {
    std::vector<int> rows;
    std::vector<std::vector<double>> weights;
  };
  const std::vector<KeptWeights> &kept_weights() const { return kept_weights_; }

protected:
  // Maps a loss name to `loss` and `calcLoss`; throws on unknown names.
  void set_loss_function(const std::string &loss);

private:
  double delta;
  double epsilon;
  enum LossType {
    L1,
    L2,
    median,
    logit,
    logit_2,
    logit_3,
    logit_4,
    exponential,
    exponential_2,
    exponential_3
  };
  LossType loss;
  void (ClassificationRPF::*calcLoss)(Split &);
  void create_tree_family(std::vector<Leaf> initial_leaves, size_t n) override;
  std::vector<double> null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows) override;
  Split calcOptimalSplit(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  // Mode-specific split calculators (classification versions using calcLoss and weights)
  Split calcOptimalSplit_leaves(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  Split calcOptimalSplit_curTrees1(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  Split calcOptimalSplit_curTrees2(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  // Mode 4: histogram-binned (classification variant)
  Split calcOptimalSplit_hist(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  Split calcOptimalSplit_resTrees(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<RandomPlantedForest::ResultingTreeCandidate> &possible_trees,
      TreeFamily &curr_family,
      std::vector<std::vector<double>> &weights
  );
  // Helpers shared by the split modes above. They mirror the regression split modes
  // (splits_*.cpp): same candidate pool, same sampling (and RNG draws), same aging;
  // only the loss is evaluated with calcLoss.
  std::vector<size_t> sampleCandidates(const std::vector<double> &ages, bool clamp_to_positive = false);
  void evaluateSplit(
      Split &curr_split, const Leaf &leaf, int k, double split_point, const std::vector<std::vector<double>> &X
  );
  void ensureThresholdWeights(
      const std::shared_ptr<DecisionTree> &tree, int kdim, const std::vector<std::vector<double>> &X
  );
  void updateThresholdWeight(
      const std::shared_ptr<DecisionTree> &tree, size_t leaf_idx, bool append, const std::vector<std::vector<double>> &X
  );
  void L1_loss(Split &split);
  void median_loss(Split &split);
  void logit_loss(Split &split);
  // logit_loss() without LogitTerms: every term computed here
  void logit_loss_direct(Split &split);
  void logit_loss_2(Split &split);
  void logit_loss_3(Split &split);
  void logit_loss_4(Split &split);
  void exponential_loss(Split &split);
  void exponential_loss_2(Split &split);
  void exponential_loss_3(Split &split);
  // The exponential loss for K >= 3 classes (value_size = K score columns, the classes coded 1 / -1): SAMME (Zhu, Zou,
  // Rosset and Hastie 2009), the loss exp(-y.f / K) with y_k = 1 for the row's class and -1 / (K - 1) otherwise, and
  // sum_k f_k = 0. The forest stores g = f / (K - 1), so predict()'s softmax of g is the loss's population link. Then
  // y.f / K = g_c (c the row's class), and every row has one weight w = exp(-g_c), kept in each of its K columns. A
  // leaf with the class weight sums W_k and shares p_k = W_k / W is moved by the minimiser of sum_k W_k exp(-a_k) over
  // sum_k a_k = 0: a_k = log p_k - mean_j log p_j. samme_loss() is the split criterion: per side, the loss change
  // sum_k W_k exp(-a_k) - W at the shares clamped to [delta, 1 - delta]; it leaves the shares in M_s / M_b. The update
  // (create_tree_family) takes the step from the shares clamped to [epsilon, 1 - epsilon] and multiplies each row's
  // weight by exp(-a_c). Binary data keep the code of exponential_loss(), which is this loss for K = 2.
  void samme_loss(Split &split);
  // a_k = log c(p_k) - mean_j log c(p_j), with c the clamp to [clamp, 1 - clamp]
  void samme_step(const std::vector<double> &share, double clamp, std::vector<double> &a) const;
  bool keep_weights_ = false;
  std::vector<KeptWeights> kept_weights_;
};

#endif
