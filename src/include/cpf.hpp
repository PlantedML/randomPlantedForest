#ifndef CPF_H
#define CPF_H

#include "rpf.hpp"

class ClassificationRPF : public RandomPlantedForest {

public:
  using RandomPlantedForest::calcOptimalSplit;
  ClassificationRPF(const std::string &loss, const RPFParams &params, double delta = 0.1, double epsilon = 0);
  void fit() override;
  void set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values) override;
  ~ClassificationRPF() {};

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
  void logit_loss_2(Split &split);
  void logit_loss_3(Split &split);
  void logit_loss_4(Split &split);
  void exponential_loss(Split &split);
  void exponential_loss_2(Split &split);
  void exponential_loss_3(Split &split);
};

#endif
