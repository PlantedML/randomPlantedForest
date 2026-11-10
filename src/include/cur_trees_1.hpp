// Helpers of split_structure "cur_trees_1" (lib/splits_cur_trees_1.cpp), for regression and classification alike.
//
// A tree family keeps the dense rank of each of its working rows in each feature (FamilyRanks, held by the family's
// null tree while the family grows). Ranks order rows exactly as their values do, ties included (-0 and +0 share a
// rank), so
//   - the number of distinct values of a leaf in feature k is its number of distinct ranks: one pass over its rows,
//     with no sort and nothing stored per row;
//   - a threshold test x < sp is rank(x) < rank(sp).
// The ranks are a memo of the working rows: they do not change the fit.
#ifndef CUR_TREES_1_HPP
#define CUR_TREES_1_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "trees.hpp"

namespace rpf_ct1 {

struct FamilyRanks {
  const std::vector<std::vector<double>> *X = nullptr; // the working rows the ranks were built from
  size_t n = 0;
  size_t p = 0;
  std::vector<std::vector<int32_t>> rank; // rank[k][row]: the row's rank among the distinct values of feature k
  std::vector<std::vector<double>> value; // value[k][rank]: the distinct values of feature k, ascending
  std::vector<uint32_t> stamp;            // scratch: stamp[rank] == epoch marks a rank as seen in the current pass
  uint32_t epoch = 0;
  std::vector<int32_t> buf; // scratch: a leaf's distinct ranks
  uint32_t next_epoch();
};

// The family's ranks: those in `slot` (the null tree's) if they were built from X, else built now and stored there.
FamilyRanks &family_ranks(std::shared_ptr<FamilyRanks> &slot, const std::vector<std::vector<double>> &X);

// The number of distinct values of feature k among the leaf's rows (cached in the leaf's unique_count_cache).
size_t unique_count(FamilyRanks &fr, Leaf &leaf, int k);

// The leaf's sampling weight in feature k, its number of thresholds the search may draw, without sorting the leaf:
// with min_rows <= 1, max(0, distinct values - 2 * leaf_size + 1), as rpf_utils::leaf_threshold_weight() gives it;
// with min_leaf_rows (min_rows > 1), max(0, hi_excl - lo) of rows_range(). A leaf without a threshold weighs 0, so it
// is never drawn.
double leaf_weight(FamilyRanks &fr, Leaf &leaf, int k, int leaf_size, int min_rows = 1);

// min_leaf_rows (min_rows > 1): the distinct-value indices [lo, hi_excl) of the leaf in feature k (indices into
// unique_values()) whose thresholds leave at least leaf_size distinct values and min_rows rows on each side, as
// rpf_utils::min_rows_range() gives them from the sorted values. From the ranks in O(rows): the min_rows-th smallest
// and largest rank, then the number of distinct ranks up to each. Cached in the leaf's rows_range_cache.
void rows_range(FamilyRanks &fr, Leaf &leaf, int k, int leaf_size, int min_rows, int &lo, int &hi_excl);

// The leaf's distinct values of feature k, ascending, as rpf_utils::compute_unique_sorted_values() gives them from the
// leaf's stably sorted values: a run of equal values is represented by its first row in the leaf's order (this
// matters only for -0 and +0). Cached in the leaf's unique_vals_cache, and the count in the leaf's unique_count_cache.
const std::vector<double> &unique_values(FamilyRanks &fr, const std::vector<std::vector<double>> &X, Leaf &leaf, int k);

// The rank of value v, a value of feature k among the working rows.
int32_t rank_of(const FamilyRanks &fr, int k, double v);

// One threshold of one leaf in feature k (a "lane"): the leaf's rows, the ranks of feature k and the threshold's
// rank. lane_sums() fills in the number of rows with rank < thr and the sums of Y[, 0] over the rows below and not
// below, each accumulated in the leaf's row order, as a row-by-row pass with a branch would. Y has one column.
struct Lane {
  const std::vector<int> *rows;
  const int32_t *rk;
  int32_t thr;
  size_t ns;
  double s, b;
};
void lane_sums(std::vector<Lane> &lanes, const std::vector<std::vector<double>> &Y);

} // namespace rpf_ct1

#endif // CUR_TREES_1_HPP
