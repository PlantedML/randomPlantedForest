
#include "trees.hpp"

std::set<int> DecisionTree::get_split_dims() const
{
  return split_dims;
}

std::vector<Leaf> DecisionTree::get_leaves() const
{
  return leaves;
}

// ----------------- helper functions -----------------

/**
 * \brief Check whether a tree with specified split_dims already exists in tree_family
 *
 * \param split_dims defining the tree to be searched for.
 * \param tree_family the family to be tested whether containing the tree.
 */
std::shared_ptr<DecisionTree> treeExists(const std::set<int> &split_dims, TreeFamily &tree_family)
{
  if (tree_family.find(split_dims) != tree_family.end())
    return tree_family[split_dims];
  return nullptr;
}

/**
 * \brief Check whether a tree with resulting_dims for a split_coordinate is already in possible_splits
 *
 * \param dim defining the dimension of the split.
 * \param possible_splits containing all possible splits.
 * \param resulting_dims as union set of split dimension and dimensions of tree which is splitted.
 */
bool possibleExists(const int dim, const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, const std::set<int> &resulting_dims)
{
  for (auto &elem : possible_splits)
  {
    if (elem.first == dim && elem.second->get_split_dims() == resulting_dims)
      return 1;
  }
  return 0;
}

/**
 * \brief Check whether a tree has a leaf with specific interval.
 *
 * \param interval to be compared with.
 * \param tree to be searched for leaf with interval.
 */
bool leafExists(std::vector<Interval> &intervals, const std::shared_ptr<DecisionTree> tree)
{
  bool exists = false;
  for (auto &leaf : tree->get_leaves())
  {

    bool same_intervals = true;
    for (unsigned int i = 0; i < intervals.size(); ++i)
    {
      if (leaf.intervals[i] != intervals[i])
      {
        same_intervals = false;
        break;
      }
    }

    if (same_intervals)
      exists = true;
  }
  return exists;
}
