#ifndef RANDOM_PLANTED_FOREST_H
#define RANDOM_PLANTED_FOREST_H

#include <decisionTree.h>
#include <optional>
#include <set>
#include <map>
#include <limits>
#include <cmath>


typedef std::vector<DecisionTree> TreeFamily;

const double INF = std::numeric_limits<double>::infinity();

struct Split {
        double min_sum;			// minimal achievable sum of squared residuals
        int tree_index;			// index of the tree
        int interval_index;		// index of the interval
        int split_coordinate;           // coordinate for splitting
        double split_point;		// splitpoint
        Split(): min_sum(INF), tree_index(0), interval_index(0), split_coordinate(0), split_point(0) {};
};

class RandomPlantedForest {
	
	public:
                RandomPlantedForest(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                    int max_interaction=2, int n_trees=50, int n_splits=30, std::optional<std::vector<int>> n_leaves=std::nullopt,
                                    int split_try=10, double t_try=0.4, std::optional<std::vector<std::vector<int>>> variables=std::nullopt);
		void fit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X);
		double predict(const std::vector<double> &X);
                std::vector<double> predict(const std::vector<std::vector<double>> &X);
		// todo: getter/setter
		
	private:
		int max_interaction;
		int n_trees;
		int n_splits;
                std::vector<int> n_leaves;
		double t_try;
		int split_try;
		int feature_size;
		int sample_size;
		std::vector<std::vector<int>> variables;
		std::vector<double> upper_bounds;
		std::vector<double> lower_bounds;
		std::vector<TreeFamily> tree_families;
                Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                       const std::multimap<int, DecisionTree*> &possible_splits, TreeFamily &curr_family);
};

#endif // RANDOM_PLANTED_FOREST_H
