#ifndef RANDOM_PLANTED_FOREST_H
#define RANDOM_PLANTED_FOREST_H

#include <decisionTree.h>
#include <optional>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <memory>


typedef std::vector<std::shared_ptr<DecisionTree>> TreeFamily;

const double INF = std::numeric_limits<double>::infinity();

struct Split {
        double min_sum;			// minimal achievable sum of squared residuals
        std::shared_ptr<DecisionTree> tree_index;       // pointer to tree // todo: check if needed or only access to interval
        std::shared_ptr<Leaf> leaf_index;		// pointer to leaf containing interval
        int split_coordinate;           // coordinate for splitting
        double split_point;		// splitpoint
        std::set<int> I_s;              // individuals smaller than splitpoint
        std::set<int> I_b;              // individuals bigger than splitpoint
        double I_s_mean;                // mean of individuals smaller than splitpoin
        double I_b_mean;                // mean of individuals bigger than splitpoint
        Split(): min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), I_s_mean(0.0), I_b_mean(0.0) {};
};

std::shared_ptr<DecisionTree> treeExists(const std::set<int> split_dims, TreeFamily &tree_family);

bool possibleExists(const Split &curr_split, const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, const std::set<int> &resulting_dims);

class RandomPlantedForest {
	
    public:
        RandomPlantedForest(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                            int max_interaction=2, int n_trees=50, int n_splits=30, std::vector<int> n_leaves=std::vector<int>(),
                            int split_try=10, double t_try=0.4, std::vector<std::vector<int>> variables=std::vector<std::vector<int>>(), bool purify_forest=false);
        void fit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X);
        double predict(const std::vector<double> &X);
        std::vector<double> predict(const std::vector<std::vector<double>> &X);
        void purify();
        // todo: getter/setter
        std::vector<TreeFamily> get_forest();

    private:
        int max_interaction;                        //
        int n_trees;                                //
        int n_splits;                               // number of performed splits for each tree family
        std::vector<int> n_leaves;                  //
        double t_try;                               //
        int split_try;                              //
        int feature_size;                           // number of feature dimension in X
        int sample_size;                            // number of samples of X
        bool purify_forest;                         // whether the forest should be purified
        bool purified = false;                      // track if forest is currently purified
        std::vector<std::vector<int>> variables;    // split dimensions for initial trees
        std::vector<double> upper_bounds;           //
        std::vector<double> lower_bounds;           //
        std::vector<TreeFamily> tree_families;      // random planted forest conatining result
        Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                               const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family);
};


#endif // RANDOM_PLANTED_FOREST_H
