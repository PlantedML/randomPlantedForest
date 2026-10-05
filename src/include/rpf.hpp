// Public API for the Random Planted Forest (regression base). The core uses
// only standard C++ types; language bindings (the Rcpp module in
// `src/rcpp_interface.h`) convert their own types at the boundary.
//
// Key entry points:
// - ctor(params): configure; set_data(Y, X) loads data; fit() trains
// - predict_matrix/predict_vector(): batch/single predictions
// - purify(): optional post-processing to orthogonalize components
// - get_model()/set_model(), get_grid_leaves()/set_grid_leaves(): export and
//   restore the forest for serialization
// - verbose_out / seed_source: injected by the binding (console, RNG)
//
// Implementation notes:
// - Training orchestrated in `lib/training.cpp`
// - Prediction logic in `lib/predict.cpp`
// - Split calculators in `lib/splits_*.cpp`
// - Utilities (RNG, sampling, caching) in `lib/internal_utils.cpp`
#ifndef RPF_H
#define RPF_H

#include "trees.hpp"
#include <cstdint>
#include <functional>

typedef std::vector<std::vector<double>> Matrix2D; /**< row-major: rows x columns */

struct RPFParams
{
  int max_interaction = 1;
  int n_trees = 50;
  int n_splits = 30;
  int split_try = 10;
  double t_try = 0.4;
  bool purify_forest = false;
  bool deterministic = false;
  int nthreads = 1;
  bool cross_validate = false;
  double split_decay_rate = 0.1;
  size_t max_candidates = 50;
  bool delete_leaves = true;
  int split_structure_mode = 3; /**< 0=res_trees, 1=cur_trees_2, 2=cur_trees_1, 3=leaves, 4=hist */
};

// Forest export for serialization: one entry per tree family.
struct TreeExport
{
  std::set<int> variables;
  std::vector<std::vector<double>> values;      /**< per leaf: value_size values */
  std::vector<std::vector<Interval>> intervals; /**< per leaf: feature_size intervals */
};
typedef std::vector<TreeExport> FamilyExport;

// Purified grid export: one entry per tree family.
struct GridTreeExport
{
  std::set<int> variables;
  std::vector<int> dims;
  Matrix2D values; /**< grid cells (column-major over dims) x value_size */
};
struct GridFamilyExport
{
  Matrix2D lim_list; /**< cell limits per feature, shared by the family's trees */
  std::vector<GridTreeExport> trees;
};

class RandomPlantedForest
{

public:
  RandomPlantedForest(){};
  explicit RandomPlantedForest(const RPFParams &params);
  virtual ~RandomPlantedForest(){};

  // Load or replace training data without fitting; computes bounds.
  void set_data(const Matrix2D &samples_Y, const Matrix2D &samples_X);
  // Train tree families on the loaded data.
  virtual void fit();
  // Restore shape metadata without training data (serialization path).
  void set_shape(int feature_size_in, int value_size_in, int sample_size_in,
                 const std::vector<double> &lower, const std::vector<double> &upper);
  // Store training data without fitting and without recomputing bounds.
  void set_training_data(const Matrix2D &samples_Y, const Matrix2D &samples_X);
  const Matrix2D &get_X() const { return X; }
  const Matrix2D &get_Y() const { return Y; }
  const std::vector<double> &get_lower_bounds() const { return lower_bounds; }
  const std::vector<double> &get_upper_bounds() const { return upper_bounds; }
  int get_feature_size() const { return feature_size; }
  int get_value_size() const { return (int)value_size; }
  int get_sample_size() const { return sample_size; }
  // Export and restore tree structure (serialization path).
  std::vector<FamilyExport> get_model() const;
  void set_model(const std::vector<FamilyExport> &model);
  // Export and restore per-tree purified grids; restoring sets purified = true.
  std::vector<GridFamilyExport> get_grid_leaves() const;
  void set_grid_leaves(const std::vector<GridFamilyExport> &grid);
  // Predict n rows of column-major X (n x p). `components = {0}` means the full
  // model, `{-1}` the intercept; otherwise a set of component indices with X
  // holding only those columns. `nthreads = 0` uses the forest's own setting.
  // Returns column-major n x value_size.
  std::vector<double> predict_matrix(const double *X, int n, int p, const std::set<int> &components, int nthreads = 0);
  Matrix2D predict_vector(const std::vector<double> &X, const std::set<int> &components);
  // Optional post-processing to redistribute effects across component orders.
  void purify_1();
  void purify_2();
  // Unified purifier: mode 1 = grid path, mode 2 = fast exact (KD-tree)
  void purify(int maxp_interaction, int nthreads, int mode);
  // Unified entry with explicit threading control
  void purify_fast_exact(int maxp_interaction, int nthreads);
  // Human-readable dump of forest structure.
  void print(std::ostream &out);
  // Legacy coarse CV over a few parameters; currently a no-op.
  void cross_validation(int n_sets = 4, const std::vector<int> &splits = {5, 50},
                        const std::vector<double> &t_tries = {0.2, 0.5, 0.7, 0.9},
                        const std::vector<int> &split_tries = {1, 2, 5, 10});
  // Mean-squared error over all entries.
  double MSE(const Matrix2D &Y_predicted, const Matrix2D &Y_true);
  // Inspect/update configuration; `set_parameters` refits.
  void get_parameters(std::ostream &out);
  virtual void set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values);
  bool is_purified();

  // Destination for warnings; nullptr silences them. The core never writes to
  // std::cout itself, as R forbids that in packages.
  std::ostream *verbose_out = nullptr;
  // Source of per-tree seeds, called once per tree on the calling thread when
  // fitting. Defaults to std::random_device.
  std::function<std::uint64_t()> seed_source;

protected:
  // Write a line to verbose_out, if set.
  void warn(const std::string &msg);
  // Internal per-family worker (grid-based mode 1)
  void purify_3_family(TreeFamily &curr_family, int maxp_interaction);
  // Internal per-family worker for fast exact purifier (mode 2)
  void purify_fast_exact_family(TreeFamily &curr_family, int maxp_interaction);
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
  // Per-tree seeds drawn from seed_source on the main thread, one per tree family
  std::vector<unsigned long long> tree_seeds_;
  std::vector<double> predict_single(const std::vector<double> &X, std::set<int> component_index);
  void L2_loss(Split &split);
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
