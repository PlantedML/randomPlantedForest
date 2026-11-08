// A growing family's tracker for the out-of-bag (OOB) path (oob_path.cpp)
#ifndef RPF_PATH_TRACKER_HPP
#define RPF_PATH_TRACKER_HPP

#include "rpf.hpp"
#include <map>
#include <set>
#include <vector>

struct RandomPlantedForest::PathTracker {
  std::vector<int> rows;                                          /**< the family's OOB training rows */
  std::vector<double> pred;                                       /**< [p * value_size + v]: its current prediction */
  std::map<std::set<int>, std::vector<std::vector<int>>> in_leaf; /**< per tree (its dims), per leaf: the rows (p) */
  std::vector<double> snaps; /**< [(c * rows.size() + p) * value_size + v]: the predictions at each checkpoint */
  size_t next = 0;           /**< the next checkpoint */
};

#endif // RPF_PATH_TRACKER_HPP
