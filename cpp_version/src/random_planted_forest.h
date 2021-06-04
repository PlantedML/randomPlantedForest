#ifndef RANDOM_PLANTED_FOREST_H
#define RANDOM_PLANTED_FOREST_H

#include <decision_tree.h>
#include <optional>
#include <exception>
#include <stdexcept>

typedef std::vector<Tree> tree_family;

class RandomPlantedForest {
	
	public:
		RandomPlantedForest(const std::vector<double> &Y, const std::vector<std::vector<double>> &X, 
							int max_interaction=2, int n_trees=50, int n_splits=30, std::optional<int> n_leaves=std::nullopt, 
							int split_try=10, double t_try=0.4, std::optional<std::vector<int>> variables=std::nullopt);
		void fit();
		vector<double> predict();
		
	private:
		int max_interaction;
		int n_trees;
		int n_splits;
		int n_leaves;
		double t_try;
		int split_try;
		int x_dim;
		std::vector<int> variables;
		std::vector<double> upper_bounds;
		std::vector<double> lower_bounds;
		std::vector<std::vector<double>> X;
		std::vector<double> Y;
		std::vector<tree_family> trees;
};


#endif // RANDOM_PLANTED_FOREST_H
