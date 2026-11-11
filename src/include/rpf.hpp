// Public API for the Random Planted Forest (regression base). The core uses
// only standard C++ types; language bindings (the Rcpp module in
// `src/rcpp_interface.h`) convert their own types at the boundary.
//
// Key entry points:
// - ctor(params): configure; set_data(Y, X) loads data; fit() trains
// - predict_matrix/predict_vector(): batch/single predictions;
//   predict_components(): many components in one pass
// - purify(): optional post-processing to orthogonalize components;
//   purify_oob(): the same, plus every family's components at its
//   out-of-bag training rows (from the in-bag counts of the fit)
// - purify_columns(), store_components(): the purified components at the
//   training rows, and the pieces of chosen components, without the purified
//   store of all components (the pruning; `lib/prune_columns.cpp`)
// - set_oob_path() / oob_path(): the out-of-bag error along the splits, from
//   one fit (`lib/oob_path.cpp`)
// - set_adaptive(): the next fit() chooses the number of splits out of bag
//   while growing, one family at a time (`lib/oob_path.cpp`)
// - coarsen_select(), coarsen_oob(), coarsen_eval(): the panel form of
//   purified components (coarsen.hpp, lib/coarsen.cpp)
// - get_model()/set_model(), get_grid_leaves()/set_grid_leaves(),
//   get_purified_rects()/set_purified_rects(): export and restore the forest
//   for serialization
// - verbose_out / seed_source: injected by the binding (console, RNG)
//
// Implementation notes:
// - Training orchestrated in `lib/training.cpp`
// - Prediction logic in `lib/predict.cpp`
// - Split calculators in `lib/splits_*.cpp`, for every loss (FamilyLoss; classification's in `lib/cpf.cpp`)
// - Utilities (RNG, sampling, caching) in `lib/internal_utils.cpp`
#ifndef RPF_H
#define RPF_H

#include "trees.hpp"
#include "kdtree.hpp"
#include "rects.hpp"
#include "coarsen.hpp"
#include <cstdint>
#include <functional>
#include <limits>

typedef std::vector<std::vector<double>> Matrix2D; /**< row-major: rows x columns */

struct LeavesScratch; // internal_utils.hpp

// What a loss other than regression's squared error (classification, cpf.hpp) does while a tree family grows. The
// split searches (splits_*.cpp) and create_tree_family() are the same for every loss; given a FamilyLoss they score
// the tries with it and let it update the family's rows, otherwise (nullptr: regression) they score by the sums of Y
// and subtract the leaf means from Y. One object per tree family, owned by create_tree_family().
struct FamilyLoss {
  virtual ~FamilyLoss() = default;
  // The criterion of splitting `leaf` at sp in feature k (0-based): smaller is better, and INF or NaN never wins. The
  // search has set trial.sum_s / trial.sum_b (value_size each) to the sums of the family's Y over the leaf's rows
  // below sp and not below, each accumulated in the order in which regression accumulates them; score() may fill
  // the other fields of `trial`.
  virtual double score(Split &trial, const Leaf &leaf, int k, double sp, const std::vector<std::vector<double>> &X) = 0;
  // The search's winner (min_sum finite, tree_index, leaf_index, split_coordinate, split_point, sum_s and sum_b set):
  // its I_s / I_b and the rest, as score() computes them for it.
  virtual void finalize(Split &winner, const std::vector<std::vector<double>> &X) = 0;
  // The split is done: update the family's rows (Y or the loss's own state) and set the two new leaves' increments.
  virtual void apply(
      const Split &split,
      const std::vector<std::vector<double>> &X,
      std::vector<std::vector<double>> &Y,
      std::vector<double> &inc_s,
      std::vector<double> &inc_b
  ) = 0;
  // The family (index n) has finished growing.
  virtual void finish(size_t n) = 0;
  // The running-sum scan. A loss whose criterion is a function of the sums of per-row terms over a split's two sides
  // (classification's logit and exponential) returns their number T > 0. The searches then build those sums as they
  // build the sums of Y (prefix sums along a leaf's sorted rows, sums per bin, or one pass over the leaf's rows per
  // try) and score a threshold with scan_score(), at a cost that does not grow with the leaf's rows where the sums are
  // prefix sums. 0 (the default): score() for every try.
  virtual size_t scan_terms() const { return 0; }
  // out[0 .. T) += row ind's terms
  virtual void add_terms(int ind, double *out) const
  {
    (void)ind;
    (void)out;
  }
  // out[i T .. (i + 1) T) = the sums of the terms of the rows leaf.individuals[order[0 .. i]], i < m
  virtual void prefix_terms(const Leaf &leaf, const size_t *order, size_t m, double *out) const
  {
    (void)leaf;
    (void)order;
    (void)m;
    (void)out;
  }
  // The parent region's part of the criterion, from its term sums t and n rows (once per leaf and feature)
  virtual double region_value(const double *t, size_t n) const
  {
    (void)t;
    (void)n;
    return 0;
  }
  // The criterion of the split whose sides have the term sums s (ns rows below the threshold) and b (nb rows), vp the
  // parent's region_value(): smaller is better, INF if the split gains nothing
  virtual double scan_score(const double *s, const double *b, size_t ns, size_t nb, double vp) const
  {
    (void)s;
    (void)b;
    (void)ns;
    (void)nb;
    (void)vp;
    return std::numeric_limits<double>::infinity();
  }
};

struct RPFParams {
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
  int min_leaf_rows = 1;        /**< rows on each side of every split (1: no limit), every split structure */
  bool exact_draw = true;       /**< leaves: the exact draw of candidates (false: the fast draw, LeavesSampler) */
  /** Allowed feature sets (1-based): a tree on the features T is grown only if T is a subset of one of them.
   *  Empty = unrestricted (only max_interaction applies). */
  std::vector<std::set<int>> interaction_constraints;
};

// Forest export for serialization: one entry per tree family.
struct TreeExport {
  std::set<int> variables;
  std::vector<std::vector<double>> values;      /**< per leaf: value_size values */
  std::vector<std::vector<Interval>> intervals; /**< per leaf: feature_size intervals */
};
typedef std::vector<TreeExport> FamilyExport;

// Purified grid export: one entry per tree family.
struct GridTreeExport {
  std::set<int> variables;
  std::vector<int> dims;
  Matrix2D values; /**< grid cells (column-major over dims) x value_size */
};
struct GridFamilyExport {
  Matrix2D lim_list; /**< cell limits per feature, shared by the family's trees */
  std::vector<GridTreeExport> trees;
};

// Sizes of the rectangle store of purify(mode = 2), for diagnostics.
struct PurifyStats {
  bool rectangles = false; /**< whether the forest holds a rectangle store */
  double n_components = 0, n_rects = 0, n_groups = 0, n_index_nodes = 0, n_table_cells = 0;
  double bytes = 0; /**< bytes allocated for the rectangles, tables and indices */
};

// The out-of-bag (OOB) components of purify_oob(): before the families are merged, every family's components are
// evaluated at its OOB training rows (in-bag count 0) and summed per row and component, separately over the families
// of even and odd index (0-based), in family order. All matrices are column-major.
struct OobComponents {
  std::vector<std::vector<int>> sets;         /**< the components (1-based features, ascending) */
  std::vector<double> sum0, sum1;             /**< n x (sets.size() * value_size): the two sums */
  std::vector<double> intercept0, intercept1; /**< n x value_size: the sums of the families' intercepts */
  std::vector<int> count;                     /**< n x 2: the number of OOB families of each half */
};

// The purified components at the training rows (purify_columns()), from one pass over the families without the
// purified store: each family is purified as purify(mode = 2) purifies it, its components are evaluated at all
// training rows and added up, and the family's pieces are freed. All matrices are column-major.
struct ComponentColumns {
  OobComponents oob;             /**< the sets; with `oob`, also the OOB sums and counts (otherwise empty) */
  std::vector<double> own;       /**< n x (sets.size() * value_size): the forest's components at the training rows */
  std::vector<double> intercept; /**< value_size: the forest's purified intercept, summed as purify() sums it */
};

// The pieces of chosen purified components (store_components()): their rectangles merged over the families exactly as
// purify(mode = 2) merges and stores them, without the other components. A component's value is the sum over the
// families: divide by n_trees.
struct ComponentStore {
  int p = 0, value_size = 1, n_trees = 1;
  std::vector<std::vector<double>> cuts; /**< the cells of every feature (0-based), as the purified store's */
  rpf_rect::Accumulator pieces;          /**< per component (1-based S): its constant and its (j, rectangles) groups */
};

// A ComponentStore built for queries: the groups' tables and k-d indices, as the purified store builds them
class BuiltStore {
public:
  // Takes the pieces out of `store`; throws std::invalid_argument if the store is inconsistent
  BuiltStore(ComponentStore &store, unsigned build_threads);
  // The components `sets` (1-based) at the n rows of column-major X (n x p, all features), as predict_components()
  // of the purified forest gives them, bitwise: column-major n x (sets.size() * value_size). Throws if a set is not
  // in the store. nthreads: at least 1.
  std::vector<double>
  eval(const double *X, int n, int p, const std::vector<std::vector<int>> &sets, int nthreads) const;
  int value_size() const { return vs; }

private:
  rpf_rect::PurifiedForest pf;
  int vs = 1, n_trees = 1;
};

// The out-of-bag (OOB) path of a fit with checkpoints (set_oob_path()): every family's prediction at each of its OOB
// training rows (in-bag count 0) after each checkpoint's number of splits, summed per row over the families of even
// and of odd index (0-based), in family order. All matrices are column-major.
struct OobPath {
  std::vector<int> checkpoints;   /**< the split counts, ascending */
  std::vector<double> sum0, sum1; /**< n x (checkpoints.size() * value_size), column c * value_size + v */
  std::vector<int> count;         /**< n x 2: the number of OOB families of each half */
};

// The adaptive fit (set_adaptive()): after every family's prediction at checkpoint c (0-based) is known, the stopper
// gets the OOB path at that checkpoint alone (as OobPath, one checkpoint: value_size columns) and returns 0 to grow on,
// or the 1-based index (1..c+1) of the chosen checkpoint to stop there. It is called on the calling thread, checkpoint
// by checkpoint.
typedef std::function<int(size_t c, const OobPath &path)> AdaptiveStopper;

// What the last adaptive fit did (adaptive_stats()), for diagnostics
struct AdaptiveStats {
  int pilot_bound = 0;          /**< the bound (checkpoints) the first `nthreads` families set together */
  double growths = 0;           /**< families grown (to a bound, or further when the stop came later) */
  double regrowths = 0;         /**< of these, grown again from their seed */
  double splits_grown = 0;      /**< splits grown in all, regrowths included */
  double pauses = 0;            /**< families paused at the bound */
  double resumes = 0;           /**< paused families that grew on */
  double max_paused_gb = 0;     /**< the largest estimated size of the paused families */
  double seconds_grow = 0;      /**< growing (the calling thread's wall time) */
  double seconds_stopper = 0;   /**< the stopper's calls */
  double seconds_rebuild = 0;   /**< rebuilding the chosen forest from the records */
  std::vector<int> bound_trace; /**< the bound after each estimate */
};

namespace rpf_cnt {
class SubsetCounter;
}

// Runs fn(family index) for n_families families, at most threads_to_use at a time, in batches of increasing index;
// after_batch(start, batch) runs on the calling thread after each batch. An exception in fn is rethrown on the calling
// thread. (lib/purify.cpp)
void run_families(
    size_t n_families,
    unsigned threads_to_use,
    const std::function<void(size_t)> &fn,
    const std::function<void(size_t, size_t)> &after_batch
);

// The components comps[k] (nullptr: zero) of the rectangle store pf at the n rows of column-major X (all features of
// pf), divided by n_trees, into out (column-major n x (comps.size() * vs), zero on entry), on up to `threads` threads:
// per chunk of rows the cells of all features once, then the components. (lib/predict.cpp)
void rects_components(
    const rpf_rect::PurifiedForest &pf,
    const std::vector<const rpf_rect::Component *> &comps,
    const double *X,
    int n,
    int vs,
    int n_trees,
    unsigned threads,
    double *out
);

class RandomPlantedForest {

public:
  RandomPlantedForest() {};
  explicit RandomPlantedForest(const RPFParams &params);
  virtual ~RandomPlantedForest() {};

  // Load or replace training data without fitting; computes bounds.
  void set_data(const Matrix2D &samples_Y, const Matrix2D &samples_X);
  // Train tree families on the loaded data.
  virtual void fit();
  // Restore shape metadata without training data (serialization path).
  void set_shape(
      int feature_size_in,
      int value_size_in,
      int sample_size_in,
      const std::vector<double> &lower,
      const std::vector<double> &upper
  );
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
  // Export and restore per-tree purified grids (mode 1); restoring sets purified = true.
  std::vector<GridFamilyExport> get_grid_leaves() const;
  void set_grid_leaves(const std::vector<GridFamilyExport> &grid);
  // Export and restore the rectangle components of purify(mode = 2) (nullptr if there are none); restoring sets
  // purified = true.
  const rpf_rect::PurifiedForest *get_purified_rects() const { return purified ? rect_forest_.get() : nullptr; }
  void set_purified_rects(rpf_rect::PurifiedForest pf);
  // The purified grids of one tree family (0-based), as one element of get_grid_leaves(): the stored grids of mode 1,
  // or the components of mode 2 evaluated on the family's grid (this needs the training data).
  GridFamilyExport get_grid_family(int family) const;
  // Sizes of the rectangle store of purify(mode = 2), for diagnostics; `rectangles` is false without one.
  PurifyStats purify_stats() const;
  // Predict n rows of column-major X (n x p). `components = {0}` means the full
  // model, `{-1}` the intercept; otherwise a set of component indices with X
  // holding only those columns. `nthreads = 0` uses the forest's own setting.
  // Returns column-major n x value_size.
  std::vector<double> predict_matrix(const double *X, int n, int p, const std::set<int> &components, int nthreads = 0);
  Matrix2D predict_vector(const std::vector<double> &X, const std::set<int> &components);
  // The components with the feature sets `sets` (1-based) at the n rows of column-major X (n x p, all features), in
  // one pass. Returns column-major n x (sets.size() * value_size), column k * value_size + c for set k and class c;
  // each column equals predict_matrix() of that component. The components a purified forest does not store are
  // zero and not evaluated.
  std::vector<double>
  predict_components(const double *X, int n, int p, const std::vector<std::vector<int>> &sets, int nthreads = 0);
  // The feature sets (1-based) of the components a purified forest stores; every other component is zero.
  std::vector<std::vector<int>> purified_components() const;
  // Optional post-processing to redistribute effects across component orders: mode 1 = dense grids per tree,
  // mode 2 = exact KD-tree counts, the components kept as weighted rectangles (rects.hpp)
  void purify(int maxp_interaction, int nthreads, int mode);
  // Human-readable dump of forest structure.
  void print(std::ostream &out);
  // Legacy coarse CV over a few parameters; currently a no-op.
  void cross_validation(
      int n_sets = 4,
      const std::vector<int> &splits = {5, 50},
      const std::vector<double> &t_tries = {0.2, 0.5, 0.7, 0.9},
      const std::vector<int> &split_tries = {1, 2, 5, 10}
  );
  // Mean-squared error over all entries.
  double MSE(const Matrix2D &Y_predicted, const Matrix2D &Y_true);
  // Inspect/update configuration; `set_parameters` refits.
  void get_parameters(std::ostream &out);
  virtual void set_parameters(const std::vector<std::string> &keys, const std::vector<double> &values);
  bool is_purified();
  // Whether the last purify() computed the components only up to an order below max_interaction.
  bool is_purify_capped() const;
  // Restore that state after set_grid_leaves() (serialization path). The order is not known: it counts as 1.
  void set_purify_capped(bool capped);
  // The order the last purify() was capped at, or 0 if it computed all orders.
  int get_purify_cap() const;
  // Restore it after set_grid_leaves() (serialization path); 0 is a purification of all orders.
  void set_purify_cap(int cap);
  // In-bag counts: inbag[f][i] = how often training row i is in the bootstrap sample of family f (saturating at 255;
  // all 1 for a deterministic forest). Empty before a fit, or after set_model() until set_inbag() (serialization).
  const std::vector<std::vector<uint8_t>> &get_inbag() const { return inbag_; }
  void set_inbag(std::vector<std::vector<uint8_t>> inbag);
  // purify(maxp_interaction, nthreads, 2), and the out-of-bag components of the training rows (see OobComponents).
  // The purified forest is the same as after purify(). Needs the training data and the in-bag counts.
  OobComponents purify_oob(int maxp_interaction, int nthreads);
  // The purified components at the training rows without the purified store (see ComponentColumns), for the pruning:
  // each family is purified as purify_oob() purifies it (the same range counts and weights), its components are
  // evaluated at all training rows tree by tree (no rectangles) and added to the columns on the calling thread in
  // family order, so nothing depends on the number of threads. The values equal those of purify_oob() and
  // predict_components() of the purified forest up to rounding; the intercept is the store's, bitwise. With `oob`, the
  // OOB sums as purify_oob() gives them (needs the in-bag counts). The forest is not changed (an earlier purification
  // stays as it is). Needs the training data.
  ComponentColumns purify_columns(int maxp_interaction, int nthreads, bool oob);
  // The pieces of the components `sets` (1-based; see ComponentStore), from one pass over the families restricted to
  // them. The forest is not changed. Needs the training data.
  ComponentStore store_components(const std::vector<std::vector<int>> &sets, int nthreads) const;
  // Track the OOB path in the next fit(): `checkpoints` are increasing split counts in 1..splits; empty turns it off
  // (the default: a plain fit does none of this). The fit itself is unchanged. A classification forest tracks its
  // value_size score columns (the leaves' values, on the link scale for logit and exponential).
  void set_oob_path(const std::vector<int> &checkpoints);
  // Whether the last fit tracked an OOB path; oob_path() returns it (see OobPath).
  bool has_oob_path() const;
  OobPath oob_path() const;
  // The next fit() is adaptive (not deterministic): the families grow one at a time, `nthreads` together,
  // and `stopper` decides at each of the `checkpoints` (increasing split counts, the last one splits) whether growing
  // goes on, from the corrected OOB error (see AdaptiveStopper). The forest at the chosen checkpoint is rebuilt from
  // the records of the splits, and n_splits becomes its number of splits: the forest that a fit with that many splits
  // and the same seed gives. `patience` (checkpoints), `tolerance` (paired standard errors) and `floor` (splits) are
  // the stopper's; with them the fit estimates how far to grow each family. `pause_gb`: the memory (GB) that families
  // waiting at that estimate may hold. One fit only: afterwards fits are plain again.
  void set_adaptive(
      const std::vector<int> &checkpoints,
      AdaptiveStopper stopper,
      int patience,
      double tolerance,
      double floor,
      double pause_gb
  );
  AdaptiveStats adaptive_stats() const;
  // The corrected out-of-bag error of training row i whose OOB score is g (value_size values), as the adaptive fit
  // estimates its bound: the per-row criterion of the forest's loss (regression, and classification with L1 / L2: the
  // squared error summed over the score columns; logit and exponential: the log-loss of the probabilities predict()
  // makes of g) minus coef times its second-order term, with H the Hessian of the per-row criterion in the score: given
  // d (the difference of the even- and odd-family OOB means), d'Hd / (2 inv); given var (the per-column variances
  // across families), sum_v H_vv var_v / 2; neither: none. For the squared error H = 2I, the regression correction.
  virtual double
  oob_row_error(size_t i, const double *g, const double *d, double inv, const double *var, double coef) const;
  // oob_row_error() of every training row, for checks: G and D column-major n x value_size (a row of D with a NaN: no
  // d), inv and coef per row
  std::vector<double> oob_row_errors(
      const std::vector<double> &G,
      const std::vector<double> &D,
      const std::vector<double> &inv,
      const std::vector<double> &coef
  ) const;
  // How the score columns are centred across the classes before their variances are summed (the size of a
  // component, the variance a coarsened form keeps): 0 not (regression, one column, L1 / L2), 1 with the reference
  // score 0 joined (multiclass logit), 2 across the columns (multiclass exponential)
  virtual int score_centring() const { return 0; }
  // The panel form of purified components (coarsen.hpp, lib/coarsen.cpp). The forest must hold its training data; the
  // components' rectangles come from passes over the families restricted to them. With value_size V > 1 (multiclass)
  // every component has V score columns, and one form (C, q, cuts) serves them all.
  // coarsen_select(): for each set S (1-based) of order >= 3, the panels the forest's f_S needs (forest_panels) and,
  // if more than P, every candidate (coarse features C, maximal q vector, equal-mass cuts) scored by the variance kept
  // on the training rows (of the columns centred by score_centring(), summed); the best one, its training-row values
  // and its stored forms (one per column). If at most P, its panels: the |S| - 2 features with the fewest own cut
  // points and these cut points (one set for all columns). fvals: column-major n_rows x n_cols, the forest's f_S at
  // the training rows, V columns per set (column k * V + v, as purify_columns() gives them).
  std::vector<rpf_coarse::Selected> coarsen_select(
      const std::vector<std::vector<int>> &sets, int P, int nthreads, const double *fvals, int n_rows, int n_cols
  );
  // coarsen_oob(): each family's coarsened components at its out-of-bag rows, summed as purify_oob() sums (needs the
  // in-bag counts); V columns per spec.
  rpf_coarse::OobSums coarsen_oob(const std::vector<rpf_coarse::Spec> &specs, int nthreads);
  // coarsen_eval(): the forest's f_S coarsened on any C (possibly empty: f_S itself) at the n rows of column-major X
  // (n x p, all features). The forest must be purified in mode 2 and have one score column.
  std::vector<double> coarsen_eval(const rpf_coarse::Spec &spec, const double *X, int n, int p);

  // Destination for warnings; nullptr silences them. The core never writes to
  // std::cout itself, as R forbids that in packages.
  std::ostream *verbose_out = nullptr;
  // Source of per-tree seeds, called once per tree on the calling thread when
  // fitting. Defaults to std::random_device.
  std::function<std::uint64_t()> seed_source;
  // Called now and then on the calling thread during a long fit (the adaptive fit); it may throw to stop the fit (a
  // user interrupt). Empty: never called.
  std::function<void()> check_interrupt;

protected:
  // Write a line to verbose_out, if set.
  void warn(const std::string &msg);
  // Internal per-family worker (grid-based mode 1)
  void purify_3_family(TreeFamily &curr_family, int maxp_interaction);
  // Internal per-family worker for fast exact purifier (mode 2); adds the family's components to `acc`. The range
  // counts come from `counter` if given (the same integers), from `kdt` otherwise. With `only`, only those components
  // (1-based feature sets) are produced.
  void purify_fast_exact_family(
      const TreeFamily &curr_family,
      int maxp_interaction,
      const rpf_kd::KDTree &kdt,
      const std::vector<std::vector<double>> &cuts,
      rpf_rect::Accumulator &acc,
      const rpf_cnt::SubsetCounter *counter = nullptr,
      const std::set<std::vector<int>> *only = nullptr
  ) const;
  // The components `only` as the families write them (purify_fast_exact_family() restricted to them, in the cells of
  // `cuts` = forest_cuts()), merged over the families in index order as purify_mode2() merges them: the store's
  // rectangles before its index
  rpf_rect::Accumulator merged_components(
      const std::set<std::vector<int>> &only, int nthreads, const std::vector<std::vector<double>> &cuts
  ) const;
  // The cells of the purified store: per feature (0-based) the sorted distinct interval end points of all leaves
  std::vector<std::vector<double>> forest_cuts() const;
  // One family's purified components at n rows (`cells`: n x stride, row-major, the rows' cells in `cuts`) added to
  // out[S] (n x value_size, row-major; created for every S the family's purification writes), and its intercept to
  // `intercept` (value_size). The purification of purify_fast_exact_family() (the same range counts from `counter`
  // and weights), evaluated tree by tree without rectangles (lib/prune_columns.cpp).
  void family_columns(
      const TreeFamily &curr_family,
      int maxp_interaction,
      const std::vector<std::vector<double>> &cuts,
      const int32_t *cells,
      size_t stride,
      size_t n,
      const rpf_cnt::SubsetCounter &counter,
      std::map<std::vector<int>, std::vector<double>> &out,
      std::vector<double> &intercept
  ) const;
  // Components of the last mode-2 purification as weighted rectangles (null after mode 1 or before purify)
  std::shared_ptr<rpf_rect::PurifiedForest> rect_forest_;
  // Threads for purify(): nthreads_param if > 0, otherwise min(the forest's nthreads, the hardware concurrency)
  unsigned purify_threads(int nthreads_param);
  // The OOB sums of purify_oob(), filled family by family (see there)
  struct OobSums;
  // Mode 2 for all families on threads_to_use threads. With oob, also the OOB sums, and the range counts from
  // per-subset k-d trees (subset_count.hpp: the same counts, so the same purification).
  void purify_mode2(int maxp_interaction, unsigned threads_to_use, OobSums *oob);
  // inbag_[f][i]: how often training row i is in the bootstrap sample of family f (saturating at 255); filled by
  // create_tree_family
  std::vector<std::vector<uint8_t>> inbag_;
  // Start family n's in-bag counts (deterministic: every row once), and count one draw of row i
  void init_inbag(size_t n) { inbag_[n].assign((size_t)sample_size, deterministic ? 1 : 0); }
  void add_inbag(size_t n, int i)
  {
    uint8_t &c = inbag_[n][(size_t)i];
    if (c < 255)
      ++c;
  }
  // The OOB path (oob_path.cpp). While family n grows, its tracker (path_tracker.hpp, owned by the family's thread)
  // follows the family's prediction at its OOB rows split by split; after the family is grown, its snapshots are added
  // to the sums on the calling thread in family order (path_merge), so the sums do not depend on the threads.
  struct PathTracker;
  std::vector<int> path_checkpoints_;
  std::vector<double> path_sum_[2]; /**< as OobPath::sum0 / sum1 */
  std::vector<std::shared_ptr<PathTracker>> path_family_;
  // Clear the sums for a new fit (fit()); the checkpoints must lie in 1..n_splits
  void path_reset();
  // A tracker for family n (growing as `fam`) after its bootstrap (nullptr if the path is off). In an adaptive fit
  // it also keeps the fit of the family's null leaf on its working outcome Y (null_leaf_value()), for the rebuild of
  // a family that has no split at the chosen checkpoint.
  std::shared_ptr<PathTracker>
  path_start(size_t n, TreeFamily *fam, const std::vector<std::vector<double>> *Y = nullptr);
  // A split found by the search, before the family's leaves change: move the tracked rows of the split leaf to the two
  // new leaves and add M_s / M_b to their predictions (the change of the family's prediction); in an adaptive fit,
  // also record the split. A classification forest passes its leaves' increments inc_s / inc_b instead (for logit and
  // exponential not M_s / M_b).
  void path_split(
      PathTracker &trk,
      const Split &split,
      const std::vector<double> *inc_s = nullptr,
      const std::vector<double> *inc_b = nullptr
  ) const;
  // `k` splits are done: snapshot the predictions if k is the next checkpoint. Returns whether the family grows on
  // (always, except in an adaptive fit)
  bool path_step(PathTracker &trk, int k) const;
  // Family n is grown: its last snapshot, then add its snapshots to the sums and free its tracker
  void path_merge(size_t n);
  // The adaptive fit (oob_path.cpp; set_adaptive()): its state, while it runs and (the records) until the rebuild
  struct AloneFit;
  std::shared_ptr<AloneFit> alone_;
  AdaptiveStats alone_stats_;
  // The fit one family at a time (instead of fit()'s loop), then the forest at the chosen checkpoint
  void fit_alone(const std::vector<Leaf> &initial_leaves);
  // path_step() of the adaptive fit, `k` splits done: at a checkpoint, pause, finish or grow on (see oob_path.cpp)
  bool alone_step(PathTracker &trk, size_t k) const;
  // Grow families on threads until every one has reached the bound
  void alone_pool(const std::vector<Leaf> &initial_leaves);
  // Family n has stopped growing: keep its record and predictions, free the rest
  void alone_finish(size_t n);
  // The end of the growing: the paused families finish, all threads are joined
  void alone_end();
  // The forest at checkpoint c (0-based), rebuilt from the records
  void alone_rebuild(size_t c);
  // Free the per-leaf caches of a growing family (pure memos of its rows, recomputed on demand: the fit is the same)
  static void free_leaf_caches(TreeFamily &family);
  GridFamilyExport family_grid_from_rects(size_t f) const;
  std::vector<std::vector<double>> X; /**< Nested vector feature samples of size (sample_size x feature_size) */
  std::vector<std::vector<double>> Y; /**< Corresponding values for the feature samples */
  int max_interaction; /**< Maximum level of interaction determining maximum number of split dimensions for a tree */
  std::vector<std::set<int>> interaction_constraints_; /**< allowed feature sets (1-based); empty = unrestricted */
  // Whether a tree with these split dims (0 = the null tree, ignored) may be grown: at most
  // max_interaction dims and, if constraints are given, a subset of an allowed set.
  bool dims_allowed(const std::set<int> &dims) const;
  int n_trees;               /**< Number of trees generated per family */
  int n_splits;              /**< Number of performed splits for each tree family */
  std::vector<int> n_leaves; /**< */
  double t_try = 0.4;        /**< */
  int split_try = 10;        /**< */
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
  // purify() computed the components only up to an order below max_interaction, so the grids sum to a truncated
  // prediction: predict() of the whole forest then scans the leaves, which purification leaves unchanged
  bool purify_capped = false;
  // The order of that cap (0: all orders); a later purify() above it purifies again
  int purify_cap = 0;
  std::vector<double> predict_single(const std::vector<double> &X, std::set<int> component_index);
  std::vector<double> predict_single_rect(const std::vector<double> &X, const std::set<int> &component_index) const;
  void L2_loss(Split &split);
  // Grow tree family n from the initial leaves (the same code for every loss: see FamilyLoss)
  void create_tree_family(std::vector<Leaf> initial_leaves, size_t n);
  // The loss of a family that grows on the working rows samples_Y (sample_size rows; boot_idx: each one's row in the
  // data): nullptr for regression (the squared error); ClassificationRPF returns its loss's.
  virtual std::unique_ptr<FamilyLoss>
  family_loss(std::vector<std::vector<double>> &samples_Y, const std::vector<int> &boot_idx);
  // After a split of cur_trees_1 or cur_trees_2: the cached sampling weight of the leaf idx of `tree` in every feature
  // whose weights are cached, replaced (append = false) or appended
  void update_leaf_weight(
      TreeFamily &family, DecisionTree &tree, size_t idx, bool append, const std::vector<std::vector<double>> &X
  );
  // A family has finished growing: drop its empty trees and everything only the growing needs (the leaves' rows and
  // caches, the trees' weight caches, cur_trees_1's ranks), and release spare capacity
  static void finish_family(TreeFamily &family);
  // A family that made no split: its null leaf gets null_leaf_value() of its rows (rpf.cpp)
  void fit_null_leaf(TreeFamily &family, const std::vector<std::vector<double>> &Y);
  virtual std::vector<double> null_leaf_value(const std::vector<std::vector<double>> &Y, const std::vector<int> &rows);
  struct SplitCandidate;
  // overload possibleExists for your vector of SplitCandidate
  static bool
  possibleExists(int dim, const std::vector<SplitCandidate> &possible_splits, const std::set<int> &resulting_dims);
  // helpers for different split-structure modes
  Split calcOptimalSplit_leaves(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
  // leaves (both outcomes): this round's n_candidates entries of the pool, by the deterministic, exact or fast draw,
  // leaving out entries without a threshold that keeps min_leaf_rows_ rows on each side (splits_leaves.cpp). sc holds
  // the exact draw's keys and the fast draw's sampler, which lasts one tree family.
  std::vector<size_t> draw_leaves_candidates(
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      size_t n_candidates,
      LeavesScratch &sc
  );
  // leaves: after the drawn entries were aged, their new weights into the fast draw's sampler
  void leaves_sampler_after_aging(
      std::vector<SplitCandidate> &possible_splits, const std::vector<size_t> &sample_idxs, LeavesScratch &sc
  );
  Split calcOptimalSplit_curTrees2(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
  Split calcOptimalSplit_curTrees1(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
  struct ResultingTreeCandidate {
    std::shared_ptr<DecisionTree> tree;
    double age = 0.0;
    ResultingTreeCandidate() = default;
    explicit ResultingTreeCandidate(std::shared_ptr<DecisionTree> t) : tree(std::move(t)) {}
  };
  bool resultingTreeExists(const std::vector<ResultingTreeCandidate> &pool, const std::set<int> &dims);
  Split calcOptimalSplit_resTrees(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<ResultingTreeCandidate> &possible_trees,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
  Split calcOptimalSplit(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
  // exponential‐decay rate for split age
  double split_decay_rate_;
  size_t max_candidates_;
  // leaves: draw the candidates exactly (one random number per pool entry), or by the fast draw
  // (LeavesSampler): the same distribution, other random numbers, so other fits
  bool exact_draw_ = true;
  // LRU cap for per-leaf per-feature caches
  size_t leaf_feature_cache_cap_ = 64;
  // track each split candidate and how long it’s sat unchosen
  struct SplitCandidate {
    int dim;
    std::shared_ptr<DecisionTree> tree;
    size_t leaf_idx;
    double age = 0.0;
    // leaves: exp(-split_decay_rate_ * age) as last computed, for the age w_age (NaN: not yet). The same
    // expression as before, recomputed only when the age changes, so the draw's weights are bitwise the same.
    double w = 0.0;
    double w_age = std::numeric_limits<double>::quiet_NaN();
    // leaves and hist with min_leaf_rows_ > 1: whether the leaf has a threshold on dim with min_leaf_rows_ rows on
    // each side (-1: not yet known). Entries without one are not drawn.
    signed char rows_ok = -1;
    // legacy ctor without leaf index (defaults to 0) — keep but prefer the 4-arg form from callers
    explicit SplitCandidate(int d, std::shared_ptr<DecisionTree> t, double a = 0.0)
        : dim(d), tree(std::move(t)), leaf_idx(0), age(a)
    {
    }
    SplitCandidate(int d, std::shared_ptr<DecisionTree> t, size_t li, double a = 0.0)
        : dim(d), tree(std::move(t)), leaf_idx(li), age(a)
    {
    }
  };
  // Which split structure to use (0=res_trees, 1=cur_trees_2, 2=cur_trees_1, 3=leaves, 4=hist)
  int split_structure_mode_ = 3;
  // At least this many rows of the bootstrap sample on each side of every split (1: no limit), in every split
  // structure: leaves and hist leave out pool entries without such a threshold (rows_ok); cur_trees_1 and res_trees
  // weigh a leaf by its number of such thresholds, so a leaf without one is not drawn; cur_trees_2 skips such a leaf.
  // Each search draws its thresholds only from those that keep the rows.
  int min_leaf_rows_ = 1;

  // Histogram mode buffers
  size_t num_bins_ = 64; // total number of global bins per feature (smaller default for speed)
  // For each feature k in [0, feature_size), store K-1 cut points (ascending)
  std::vector<std::vector<double>> feature_cut_points_;
  // For each feature k, per-sample bin id in [0, K-1]
  std::vector<std::vector<int>> sample_bin_id_;
  // For the current bootstrapped working set (per-family), cache per-feature bin ids
  // Moved to thread-local storage in implementation to avoid races under multithreading
  // std::vector<std::vector<int>> working_bin_id_;

  bool delete_leaves;

  // Mode 4: histogram-binned split evaluation
  Split calcOptimalSplit_hist(
      const std::vector<std::vector<double>> &Y,
      const std::vector<std::vector<double>> &X,
      std::vector<SplitCandidate> &possible_splits,
      TreeFamily &curr_family,
      FamilyLoss *loss = nullptr
  );
};

#endif // RPF_HPP
