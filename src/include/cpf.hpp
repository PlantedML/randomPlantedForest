#ifndef CPF_H
#define CPF_H

#include <vector>
#include "rpf.hpp"
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
struct CPFParams
{
  int max_interaction;
  int n_trees;
  int n_splits;
  int split_try;
  double t_try;
  bool purify_forest;
  bool deterministic;
  int nthreads;
  bool cross_validate;
  LossType loss;
};

class ClassificationRPF : public RandomPlantedForest
{

public:
  using RandomPlantedForest::calcOptimalSplit;
  ClassificationRPF(const std::vector<std::vector<double>> &samples_Y, const std::vector<std::vector<double>> &samples_X,
                    const std::string loss = "L2", const std::vector<double> parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0, 0, 0.1});
  ~ClassificationRPF() {};
  CPFParams get_parameters();

private:
  double delta;
  double epsilon;
  LossType loss;
  void (ClassificationRPF::*calcLoss)(Split &);
  void create_tree_family(std::vector<Leaf> initial_leaves, size_t n) override;
  void fit() override;
  Split calcOptimalSplit(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                         std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family,
                         std::vector<std::vector<double>> &weights);
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
