#ifndef RPF_H
#define RPF_H

#include <vector>
#include <stdexcept>
#include <Rcpp.h>
#include "trees.hpp"
using namespace Rcpp;
class RandomPlantedForest
{

public:
  RandomPlantedForest(
      const std::vector<std::vector<double>> &samples_Y,
      const std::vector<std::vector<double>> &samples_X,
      const std::vector<double> parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0});
  RandomPlantedForest() {};
  void set_data(const std::vector<std::vector<double>> &samples_Y, const std::vector<std::vector<double>> &samples_X);
  std::vector<std::vector<double>> predict_matrix(const std::vector<std::vector<double>> &X, const std::vector<double> components = {0});
  std::vector<std::vector<double>> predict_vector(const std::vector<double> &X, const std::vector<double> components = {0});
  void purify_1();
  void purify_2();
  void purify_3();
  void print();
  void cross_validation(int n_sets = 4, std::vector<int> splits = {5, 50}, std::vector<double> t_tries = {0.2, 0.5, 0.7, 0.9}, std::vector<int> split_tries = {1, 2, 5, 10});
  double MSE(const std::vector<std::vector<double>> &Y_predicted, const std::vector<std::vector<double>> &Y_true);
  void get_parameters();
  void set_parameters(std::vector<std::string> keys, std::vector<double> values);
  List get_model();
  virtual ~RandomPlantedForest() {};
  bool is_purified();

protected:
  double MSE_vec(const std::vector<double> &Y_predicted, const std::vector<double> &Y_true);
  std::vector<std::vector<double>> X; /**< Nested vector feature samples of size (sample_size x feature_size) */
  std::vector<std::vector<double>> Y; /**< Corresponding values for the feature samples */
  int max_interaction;                /**< Maximum level of interaction determining maximum number of split dimensions for a tree */
  int n_trees;                        /**< Number of trees generated per family */
  int n_splits;                       /**< Number of performed splits for each tree family */
  std::vector<int> n_leaves;          /**< */
  double t_try = 0.4;                 /**< */
  int split_try = 10;                 /**< */
  size_t value_size = 1;
  int feature_size = 0;       /**< Number of feature dimension in X */
  int sample_size = 0;        /**< Number of samples of X */
  bool purify_forest = 0;     /**< Whether the forest should be purified */
  bool purified = false;      /**< Track if forest is currently purified */
  bool deterministic = false; /**< Choose whether approach deterministic or random */
  // bool parallelize = false;                   /**< Perform algorithm in parallel or serialized */
  int nthreads = 1;            /**< Number threads used for parallelisation */
  bool cross_validate = false; /**< Determines if cross validation is performed */
  std::vector<double> upper_bounds;
  std::vector<double> lower_bounds;
  std::vector<TreeFamily> tree_families; /**<  random planted forest containing result */
  std::vector<double> predict_single(const std::vector<double> &X, std::set<int> component_index);
  void L2_loss(Split &split);
  virtual void fit();
  virtual void create_tree_family(std::vector<Leaf> initial_leaves, size_t n);
  virtual Split calcOptimalSplit(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                 std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family);
};

#endif // RPF_HPP
