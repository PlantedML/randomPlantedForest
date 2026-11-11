#ifndef TREES_H
#define TREES_H

#include "helper.hpp"
#include "grid.hpp"
#include <memory>
#include <unordered_map>

using namespace utils;
using namespace grid;

namespace rpf_ct1 {
struct FamilyRanks; // cur_trees_1.hpp
}

class DecisionTree;
// The split search's memos for one leaf, for each feature dimension: a stable order of indices into `individuals`
// sorted by the feature value (reusable across evaluations), the sorted feature values along it (for lower_bound), the
// unique sorted values (threshold sampling), the unique count (to skip leaves with too few thresholds) and, with
// min_leaf_rows, the range of thresholds that keep enough rows on each side. They are
// functions of the leaf's rows and are needed only while the family grows.
struct LeafCache {
  std::unordered_map<int, std::vector<size_t>> order_cache;
  std::unordered_map<int, std::vector<double>> sorted_vals_cache;
  std::unordered_map<int, std::vector<double>> unique_vals_cache;
  std::unordered_map<int, size_t> unique_count_cache;
  // tree_then_leaf (cur_trees_1) with min_leaf_rows: the range [lo, hi_excl) of distinct-value indices whose
  // thresholds keep min_leaf_rows rows on each side (rpf_ct1::rows_range)
  std::unordered_map<int, std::pair<int, int>> rows_range_cache;
};

struct Leaf {
  std::vector<int> individuals;    /**< considered samples for each leaf */
  std::vector<double> value;       /**< residual */
  std::vector<Interval> intervals; /**< min/max for each feature of the interval */

  // The caches sit behind one pointer: allocated on first use, freed by free_cache() when the family has finished
  // growing, so a finished leaf holds one pointer instead of four hash maps. Copying a leaf copies its caches.
  LeafCache &cache()
  {
    if (!cache_)
      cache_.reset(new LeafCache());
    return *cache_;
  }
  void free_cache() { cache_.reset(); }
  bool has_cache() const { return static_cast<bool>(cache_); }

  Leaf() = default;
  Leaf(const Leaf &o)
      : individuals(o.individuals), value(o.value), intervals(o.intervals),
        cache_(o.cache_ ? new LeafCache(*o.cache_) : nullptr)
  {
  }
  Leaf &operator=(const Leaf &o)
  {
    if (this != &o) {
      individuals = o.individuals;
      value = o.value;
      intervals = o.intervals;
      if (o.cache_ && cache_)
        *cache_ = *o.cache_; // assigns the maps, as the four members did
      else
        cache_.reset(o.cache_ ? new LeafCache(*o.cache_) : nullptr);
    }
    return *this;
  }
  Leaf(Leaf &&) noexcept = default;
  Leaf &operator=(Leaf &&) noexcept = default;

private:
  std::unique_ptr<LeafCache> cache_;
};

/**
 * \brief A split performed with a score at a leaf_index in tree_index.
 *
 * Remembers data for the two news leaves.
 */
struct Split {
  double min_sum;                           /**< minimal achievable sum of squared residuals */
  std::shared_ptr<DecisionTree> tree_index; /**< pointer to tree */
  Leaf *leaf_index;                         /**< pointer to leaf containing interval */
  int split_coordinate;                     /**< coordinate for splitting */
  double split_point;                       /**< splitpoint */
  double M_sp;
  double M_bp;
  std::vector<double> sum_s;
  std::vector<double> sum_b;
  std::vector<int> I_s;    /**< individuals smaller than splitpoint */
  std::vector<int> I_b;    /**< individuals bigger than splitpoint */
  std::vector<double> M_s; /**< mean or median of individuals smaller than splitpoin */
  std::vector<double> M_b; /**< mean or median of individuals bigger than splitpoint */
  const std::vector<std::vector<double>> *W;
  const std::vector<std::vector<double>> *Y;
  Split()
      : min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), M_sp(0),
        M_bp(0) {};
};

/**
 * \brief Decision trees contain split data as leaves for respective splitting
 * dimensions.
 */
class DecisionTree {
public:
  DecisionTree() {};
  DecisionTree(std::set<int> dims, std::vector<Leaf> &first_leaves) : split_dims(dims), leaves(first_leaves) {};
  DecisionTree(std::set<int> dims) : split_dims(dims) {};
  std::set<int> get_split_dims() const;
  std::vector<Leaf> get_leaves() const;

private:
  friend class RandomPlantedForest;
  friend class ClassificationRPF;
  std::set<int> split_dims; /**<  dimensions of the performed splits */
  std::vector<Leaf> leaves; /**<  leaves of tree containing intervals and approximating value */
  LeafGrid GridLeaves;
  // Cached per-dimension weighted sampling over leaves for cur_trees_1
  // epoch that increments whenever leaves are structurally changed
  int weights_epoch = 0;
  // For each feature dimension k, remember which epoch the cache corresponds to
  // vector-backed caches to avoid unordered_map overhead
  std::vector<int> weights_epoch_by_dim_v; // length == feature_size (lazy-sized)
  // For each feature dimension k, Fenwick tree (1-based) of per-leaf weights (width of valid thresholds)
  std::vector<std::vector<double>> fenwick_by_dim_v; // length == feature_size (lazy-sized)
  // For each feature dimension k, raw per-leaf weights array (0-based index over leaves)
  std::vector<std::vector<double>> leaf_weights_by_dim_v; // length == feature_size (lazy-sized)
  // For each feature dimension k, total weight across all leaves
  std::vector<double> weights_total_by_dim_v; // length == feature_size (lazy-sized)
  // cur_trees_1: the family's ranks of its working rows (cur_trees_1.hpp), held by the null tree while the family
  // grows (training only)
  std::shared_ptr<rpf_ct1::FamilyRanks> ct1_ranks;
};

typedef std::map<std::set<int>, std::shared_ptr<DecisionTree>, setComp> TreeFamily;

std::shared_ptr<DecisionTree> treeExists(const std::set<int> &split_dims, TreeFamily &tree_family);
// Legacy overload kept in trees.cpp for backward compatibility in R-facing helpers.
bool possibleExists(
    const int dim,
    const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits,
    const std::set<int> &resulting_dims
);
bool leafExists(std::vector<Interval> &intervals, const std::shared_ptr<DecisionTree> tree);

#endif // TREES_H
