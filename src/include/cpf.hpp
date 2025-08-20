#ifndef CPF_H
#define CPF_H

#include "rpf.hpp"

class ClassificationRPF : public RandomPlantedForest
{

public:
  using RandomPlantedForest::calcOptimalSplit;
  ClassificationRPF(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                    const String loss = "L2", const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0, 0.1, 0, 0.1, 50,1});
  void set_parameters(StringVector keys, NumericVector values);
  ~ClassificationRPF(){};

private:
  double delta;
  double epsilon;
  enum LossType
  {
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
  void fit() override;
  Split calcOptimalSplit(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<SplitCandidate>& possible_splits,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights) ;
  // Mode-specific split calculators (classification versions using calcLoss and weights)
  Split calcOptimalSplit_leaves(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<SplitCandidate>& possible_splits,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights);
  Split calcOptimalSplit_curTrees1(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<SplitCandidate>& possible_splits,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights);
  Split calcOptimalSplit_curTrees2(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<SplitCandidate>& possible_splits,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights);
  // Mode 4: histogram-binned (classification variant)
  Split calcOptimalSplit_hist(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<SplitCandidate>& possible_splits,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights);
  Split calcOptimalSplit_resTrees(
      const std::vector<std::vector<double>>& Y,
      const std::vector<std::vector<double>>& X,
      std::vector<RandomPlantedForest::ResultingTreeCandidate>& possible_trees,
      TreeFamily& curr_family,
      std::vector<std::vector<double>>& weights);
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
