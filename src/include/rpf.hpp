// Public API for the Random Planted Forest (regression base). This header
// declares the externally visible training, prediction, and model-introspection
// methods used from R via the Rcpp module in `src/randomPlantedForest.cpp`.
//
// Key entry points:
// - ctor(Y, X, parameters): construct and fit a model (calls set_data + fit)
// - set_data(Y, X): load data (no training) and initialize bounds
// - fit(): build tree families according to split_structure_mode_
// - predict_matrix/predict_vector(): batch/single predictions
// - purify_1/2/3(): optional post-processing to orthogonalize components
// - cross_validation(): coarse k-fold search over a few parameters (legacy)
// - get_parameters()/set_parameters(): inspect or update configuration
// - get_model(): export current forest (for R printing/plotting)
// - is_purified(): flag indicating whether purify_* was applied last
//
// Implementation notes:
// - Training orchestrated in `lib/training.cpp`
// - Prediction logic in `lib/predict.cpp`
// - Split calculators in `lib/splits_*.cpp`
// - Utilities (RNG, sampling, caching) in `lib/internal_utils.cpp`
#ifndef RPF_H
#define RPF_H

#include "trees.hpp"

using namespace Rcpp;

class RandomPlantedForest
{

public:
  // Construct and fit a random planted forest on Y ~ X with configuration in
  // `parameters` (see R docs for positional mapping; last value selects
  // split-structure mode). Calls set_data() then fit().
  RandomPlantedForest(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                      const NumericVector parameters = {1, 50, 30, 10, 0.4, 0, 0, 0, 0, 0.1, 50, 1, 3});
  RandomPlantedForest(){};
  // Load or replace data without fitting; computes bounds and resets state.
  void set_data(const NumericMatrix &samples_Y, const NumericMatrix &samples_X);
  // Predict for a matrix or a single vector. `components = {0}` means the full
  // model; otherwise a set of component indices to evaluate (expert mode).
  NumericMatrix predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
  NumericMatrix predict_vector(const NumericVector &X, const NumericVector components = {0});
  // Optional post-processing to redistribute effects across component orders.
  void purify_1();
  void purify_2();
  void purify_3();
  // Human-readable dump of forest structure to R console.
  void print();
  // Legacy coarse CV over a few parameters; mainly for internal experiments.
  void cross_validation(int n_sets = 4, IntegerVector splits = {5, 50}, NumericVector t_tries = {0.2, 0.5, 0.7, 0.9}, IntegerVector split_tries = {1, 2, 5, 10});
  // Mean-squared error helper for matrix outputs.
  double MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true);
  // Inspect/update configuration; `set_parameters` may trigger a refit.
  void get_parameters();
  void set_parameters(StringVector keys, NumericVector values);
  // Export a list representation of the current forest for printing/plotting.
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
  // Seeds generated on the main thread from R's RNG, one per tree family
  std::vector<unsigned long long> tree_seeds_;
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
  // LRU cap for per-leaf per-feature caches
  size_t leaf_feature_cache_cap_ = 64;
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
  // Which split structure to use (0=res_trees, 1=cur_trees_2, 2=cur_trees_1, 3=leaves, 4=hist)
  int split_structure_mode_ = 3;
  
  // Histogram mode buffers
  size_t num_bins_ = 64; // total number of global bins per feature (smaller default for speed)
  // For each feature k in [0, feature_size), store K-1 cut points (ascending)
  std::vector<std::vector<double>> feature_cut_points_;
  // For each feature k, per-sample bin id in [0, K-1]
  std::vector<std::vector<int>> sample_bin_id_;
  // For the current bootstrapped working set (per-family), cache per-feature bin ids
  // Moved to thread-local storage in implementation to avoid races under multithreading
  // std::vector<std::vector<int>> working_bin_id_;
  
  bool leafCandidateExists(const std::vector<SplitCandidate>&,
                           const std::shared_ptr<DecisionTree>&,
                           size_t leaf_idx, int dim);
  bool delete_leaves;

  // Mode 4: histogram-binned split evaluation
  Split calcOptimalSplit_hist(const std::vector<std::vector<double>> &Y,
                              const std::vector<std::vector<double>> &X,
                              std::vector<SplitCandidate> &possible_splits,
                              TreeFamily &curr_family);
};

#endif // RPF_HPP
