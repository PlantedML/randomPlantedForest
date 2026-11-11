// A growing family's tracker for the out-of-bag (OOB) path (oob_path.cpp)
#ifndef RPF_PATH_TRACKER_HPP
#define RPF_PATH_TRACKER_HPP

#include "rpf.hpp"
#include <cstdint>
#include <map>
#include <set>
#include <vector>

// One split of a family in the adaptive fit, enough to rebuild its leaves (40 bytes): the split tree (an index into
// the family's list of trees' dims), the split leaf's index in that tree, the 1-based feature, the split point and the
// two halves' values M_s, M_b (their first columns; with value_size > 1 all columns are in the family's rec_vals)
struct RpfSplitRecord {
  int32_t dims;
  int32_t leaf;
  int32_t coord;
  double point;
  double ms;
  double mb;
};

struct RandomPlantedForest::PathTracker {
  std::vector<int> rows;                                          /**< the family's OOB training rows */
  std::vector<double> pred;                                       /**< [p * value_size + v]: its current prediction */
  std::map<std::set<int>, std::vector<std::vector<int>>> in_leaf; /**< per tree (its dims), per leaf: the rows (p) */
  std::vector<double> snaps; /**< [(c * rows.size() + p) * value_size + v]: the predictions at each checkpoint */
  size_t next = 0;           /**< the next checkpoint */
  // the adaptive fit
  size_t family = 0;
  TreeFamily *fam = nullptr;                /**< the growing family (its caches are freed when it pauses) */
  bool stopped = false;                     /**< it stopped before the last checkpoint */
  std::vector<std::vector<double>> csnaps;  /**< per checkpoint reached: the predictions at its OOB rows */
  std::vector<RpfSplitRecord> rec;          /**< its splits in order */
  std::vector<double> rec_vals;             /**< value_size > 1: per record the 2 * value_size values (M_s, then M_b) */
  std::vector<uint32_t> rec_at;             /**< per checkpoint reached: the number of records up to it */
  std::map<std::set<int>, int32_t> dims_id; /**< the split trees' dims, numbered as they first occur */
  std::vector<std::set<int>> dims_list;
  std::vector<double> null_value; /**< the fit of the null leaf if the family makes no split (path_start()) */
};

#endif // RPF_PATH_TRACKER_HPP
