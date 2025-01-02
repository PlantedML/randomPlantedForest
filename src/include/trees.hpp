#ifndef TREES_H
#define TREES_H

#include "helper.hpp"
#include "grid.hpp"
#include <vector>

using namespace utils;
using namespace grid;

class DecisionTree;
struct Leaf
{
  std::vector<int> individuals;    /**< considered samples for each leaf */
  std::vector<double> value;       /**< residual */
  std::vector<Interval> intervals; /**< min/max for each feature of the interval */
};

/**
 * \brief A split performed with a score at a leaf_index in tree_index.
 *
 * Remembers data for the two news leaves.
 */
struct Split
{
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
  Split() : min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), M_sp(0), M_bp(0){};
};

/**
 * \brief Decision trees contain split data as leaves for respective splitting
 * dimensions.
 */
class DecisionTree
{
public:
  DecisionTree(){};
  DecisionTree(std::set<int> dims, std::vector<Leaf> &first_leaves) : split_dims(dims), leaves(first_leaves){};
  DecisionTree(std::set<int> dims) : split_dims(dims){};
  std::set<int> get_split_dims() const;
  std::vector<Leaf> get_leaves() const;

private:
  friend class RandomPlantedForest;
  friend class ClassificationRPF;
  std::set<int> split_dims; /**<  dimensions of the performed splits */
  std::vector<Leaf> leaves; /**<  leaves of tree containing intervals and approximating value */
  LeafGrid GridLeaves;

};

typedef std::map<std::set<int>, std::shared_ptr<DecisionTree>, setComp> TreeFamily;

std::shared_ptr<DecisionTree> treeExists(const std::set<int> &split_dims, TreeFamily &tree_family);
bool possibleExists(const int dim, const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, const std::set<int> &resulting_dims);
bool leafExists(std::vector<Interval> &intervals, const std::shared_ptr<DecisionTree> tree);

#endif // TREES_H
