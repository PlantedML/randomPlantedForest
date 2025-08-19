#ifndef RPF_H
#define RPF_H

#include "trees.hpp"

using namespace Rcpp;

class RandomPlantedForest
{

public:
  RandomPlantedForest(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                      const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0, 0.1, 50, 1, 3});
  RandomPlantedForest(){};
  void set_data(const NumericMatrix &samples_Y, const NumericMatrix &samples_X);
  NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
  NumericMatrix predict_vector(const NumericVector &X, const NumericVector components = {0});
  void purify_1();
  void purify_2();
  void purify_3();
  void print();
  void cross_validation(int n_sets = 4, IntegerVector splits = {5, 50}, NumericVector t_tries = {0.2, 0.5, 0.7, 0.9}, IntegerVector split_tries = {1, 2, 5, 10});
  double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true);
  void get_parameters();
  void set_parameters(StringVector keys, NumericVector values);
  List get_model();
  virtual ~RandomPlantedForest(){};
  bool is_purified();
  
protected:
  double MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true);
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
  struct SplitCandidate;
  // overload possibleExists for your vector of SplitCandidate
  static bool possibleExists(
    int dim,
    const std::vector<SplitCandidate>& possible_splits,
    const std::set<int>& resulting_dims
  );
  // helpers for different split-structure modes
  Split calcOptimalSplit_leaves(const std::vector<std::vector<double>> &Y,
                                const std::vector<std::vector<double>> &X,
                                std::vector<SplitCandidate> &possible_splits,
                                TreeFamily &curr_family);
  Split calcOptimalSplit_curTrees2(const std::vector<std::vector<double>> &Y,
                                   const std::vector<std::vector<double>> &X,
                                   std::vector<SplitCandidate> &possible_splits,
                                   TreeFamily &curr_family);
  Split calcOptimalSplit_curTrees1(const std::vector<std::vector<double>> &Y,
                                   const std::vector<std::vector<double>> &X,
                                   std::vector<SplitCandidate> &possible_splits,
                                   TreeFamily &curr_family);
  struct ResultingTreeCandidate { std::shared_ptr<DecisionTree> tree; double age = 0.0; ResultingTreeCandidate() = default; explicit ResultingTreeCandidate(std::shared_ptr<DecisionTree> t):tree(std::move(t)){} };
  bool resultingTreeExists(const std::vector<ResultingTreeCandidate>& pool, const std::set<int>& dims);
  Split calcOptimalSplit_resTrees(const std::vector<std::vector<double>> &Y,
                                  const std::vector<std::vector<double>> &X,
                                  std::vector<ResultingTreeCandidate> &possible_trees,
                                  TreeFamily &curr_family);
  virtual Split calcOptimalSplit(const std::vector<std::vector<double>> &Y,
                                 const std::vector<std::vector<double>> &X,
                                 std::vector<SplitCandidate> &possible_splits,
                                 TreeFamily &curr_family);
  // exponential‐decay rate for split age
  double split_decay_rate_;
  size_t max_candidates_;
  // track each split candidate and how long it’s sat unchosen
  struct SplitCandidate {
    int dim;
    std::shared_ptr<DecisionTree> tree;
    size_t leaf_idx;
    double age = 0.0;
    // legacy ctor without leaf index (defaults to 0) — keep but prefer the 4-arg form from callers
    explicit SplitCandidate(int d, std::shared_ptr<DecisionTree> t, double a=0.0)
      : dim(d), tree(std::move(t)), leaf_idx(0), age(a) {}
    SplitCandidate(int d, std::shared_ptr<DecisionTree> t, size_t li, double a=0.0)
      : dim(d), tree(std::move(t)), leaf_idx(li), age(a) {}
  };
  // Which split structure to use (0=res_trees, 1=cur_trees_2, 2=cur_trees_1, 3=leaves)
  int split_structure_mode_ = 3;
  
  bool leafCandidateExists(const std::vector<SplitCandidate>&,
                           const std::shared_ptr<DecisionTree>&,
                           size_t leaf_idx, int dim);
  bool delete_leaves;
};

#endif // RPF_HPP
