#include <random_planted_forest.h>
#include <iostream>

RandomPlantedForest::RandomPlantedForest(const std::vector<double> &Y, const std::vector<std::vector<double>> &X, 
										 int max_interaction, int n_trees, int n_splits, std::optional<int> n_leaves, 
										 int split_try, double t_try, std::optional<std::vector<int>> variables)
{
    this->Y = Y;
    this->X = X;

    // Check if all vector in x are same length
    for(const auto &vec:X){
        this->x_dim = vec.size();
        if(vec.size() != x_dim){
            // throw std::invalid_argument("Dimensions of X mismatch.");
			std::cout << "Dimensions of X mismatch." << std::endl;
        }
    }
	if(Y.size() != X.size()){
		std::cout << "Dimensions of X and Y mismatch." << std::endl;
	}

    this->max_interaction = max_interaction;
    this->n_trees = n_trees;
    this->n_splits = n_splits;
    this->split_try = split_try;
    this->t_try = t_try;

    // if arguments not specified set to default
    if(n_leaves.has_value()){
        this->n_leaves = n_leaves.value();
    }else{
        this->n_leaves = x_dim;
    }
    if(variables.has_value()){
        this->variables = variables.value();
    }else{
        this->variables = {};
    }


	for(int n=0; n<n_trees; ++n){
		tree_family curr_fam;
		this->trees.push_back(curr_fam);
	}
}
