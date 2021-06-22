#include <randomPlantedForest.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>

RandomPlantedForest::RandomPlantedForest(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                         int max_interaction, int n_trees, int n_splits, std::optional<std::vector<int>> n_leaves,
                                         int split_try, double t_try, std::optional<std::vector<std::vector<int>>> variables){
    // Check if all vector in x are same length
    for(const auto &vec:X){
        this->feature_size = vec.size();
        if(vec.size() != feature_size) std::cout << "Dimensions of X mismatch." << std::endl;
    }

    // initilize class members
    this->max_interaction = max_interaction;
    this->n_trees = n_trees;
    this->n_splits = n_splits;
    this->split_try = split_try;
    this->t_try = t_try;

    // if arguments not specified set to default
    if(n_leaves.has_value()){
        if(n_leaves.value().size() != feature_size) std::cout << "Number of nodes for leafes has wrong dimension." << std::endl;
        this->n_leaves = n_leaves.value();
    }else{
        this->n_leaves = std::vector<int>(feature_size, 1);
    }
    if(variables.has_value()){
        this->variables = variables.value();
        // todo: check for valid values in feature range
    }else{
        this->variables = std::vector<std::vector<int>>(feature_size);
        for(int i = 0; i<feature_size; ++i) this->variables[i] = std::vector<int> {i+1};
    }

    // construct tree families
    this->fit(Y, X);
}


Split RandomPlantedForest::calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                            const std::multimap<int, DecisionTree*> &possible_splits, TreeFamily &curr_family){
    Split curr_split;
    std::set<int> tree_dims;
    int k;
    double curr_sum, leaf_size, sample_point, Y_s_mean, Y_b_mean, Y_s_sum, Y_b_sum;

    // setup random device
    std::random_device rd;
    std::mt19937 gen(rd());

    // sample possible splits
    int n_candidates = ceil(t_try*possible_splits.size()); // number of candidates that will be considered
    std::vector<int> split_candidates(possible_splits.size());
    std::iota(split_candidates.begin(), split_candidates.end(), 0); // consecutive indices of possible candidates
    std::shuffle(split_candidates.begin(), split_candidates.end(), gen); // shuffle for random order

    // go through all trees in current family
    for(DecisionTree curr_tree:curr_family){

        // consider a fraction of possible splits
        auto candidate = possible_splits.begin();
        for(size_t n=0; n<n_candidates; ++n){
            std::advance(candidate, split_candidates[n]); // get random split candidate without replacement

            k = candidate->first - 1; // split dim of current candidate, converted to index starting at 0
            leaf_size = n_leaves[k-1];

            // test if coordinate for splitting of current candidate is in current tree and results in a valid tree
            tree_dims = curr_tree.split_dims;
            tree_dims.insert(k+1); // add candidate's coordinate to current tree
            if(tree_dims != candidate->second->split_dims) continue; // if split dimensions do not match consider next candidate

            // go through all leaves of current tree
            for(Leaf leaf:curr_tree.leaves){
                std::set<int> curr_individuals = leaf.individuals; // consider individuals of current leaf

                // extract sample points according to individuals from X and Y
                std::map<double, double> unique_samples;
                for(int individual: curr_individuals){
                    unique_samples[X[individual][k]] = Y[individual];
                }

                // check if number of sample points is within limit
                if(unique_samples.size() < 2*leaf_size) continue;

                // consider split_try-number of random samples
                for(int t=0; t<split_try; ++t){

                    // get samplepoint
                    auto sample_pos = unique_samples.begin();
                    std::uniform_int_distribution<> distrib(leaf_size, unique_samples.size() + 1 - 2*leaf_size);
                    std::advance(sample_pos, distrib(gen)); // consider only sample points with offset
                    sample_point = sample_pos->first;

                    std::vector<double> Y_s; // samples smaller than the samplepoint
                    std::vector<double> Y_b; // samples bigger than the samplepoint

                    // sum of respective samples
                    Y_s_sum = 0;
                    Y_b_sum = 0;

                    // check which samples are greater/smaller than samplepoint
                    for(int individual: curr_individuals){
                        if(X[individual][k]<sample_point){
                            Y_s.push_back(Y[individual]);
                            Y_s_sum += Y[individual];
                        }else{
                            Y_b.push_back(Y[individual]);
                            Y_b_sum += Y[individual];
                        }
                    }

                    // get mean
                    Y_s_mean = Y_s_sum / Y_s.size();
                    Y_b_mean = Y_b_sum / Y_b.size();

                    // accumulate squared mean
                    curr_sum = 0;
                    std::for_each(Y_s.begin(), Y_s.end(), [&Y_s_mean, &curr_sum](double val){ curr_sum += pow(val - Y_s_mean, 2) + pow(val, 2); });
                    std::for_each(Y_b.begin(), Y_b.end(), [&Y_b_mean, &curr_sum](double val){ curr_sum += pow(val - Y_b_mean, 2) + pow(val, 2); });

                    // update split if squared sum is smaller
                    if(curr_sum < curr_split.min_sum){
                        curr_split.min_sum = curr_sum;
                        // todo: change type to pointer
                        //curr_split.tree_index = curr_tree;
                        //curr_split.interval_index = leaf;
                        curr_split.split_coordinate = k+1;
                        curr_split.split_point = sample_point;
                    }
                }
            }

        }
    }

    return curr_split;
}


void RandomPlantedForest::fit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X){

    // Check for correct input dimensions
    for(const auto &vec:X){
        this->feature_size = vec.size();
        if(vec.size() != feature_size){
			std::cout << "Dimensions of X mismatch." << std::endl;
        }
    }
    if(Y.size() != X.size()){
            std::cout << "Dimensions of X and Y mismatch." << std::endl;
    }

    this->sample_size = X.size();
    this->upper_bounds = std::vector<double>(feature_size);
    this->lower_bounds = std::vector<double>(feature_size);

    // get upper/lower bounds
    double minVal, maxVal, currVal;
    for(size_t i=0; i<feature_size; ++i){
            minVal = maxVal = X[0][i];
            for(size_t j=0; j<sample_size; ++j){
                    currVal = X[j][i];
                    if(currVal<minVal) minVal = currVal;
                    if(currVal>maxVal) maxVal = currVal;
            }
            this->upper_bounds[i] = maxVal;
            this->lower_bounds[i] = minVal;
    }

    // setup initial set of individuals
    std::set<int> initial_individuals;
    auto pos = initial_individuals.begin();
    for(int i = 0; i < sample_size-1; ++i) pos = initial_individuals.insert(pos, i);

    // initilize intervals with lower and upper bounds
    std::vector<Interval> initial_intervals(feature_size);
    for(size_t i = 0; i<feature_size; ++i) initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};

    // set properties of first leaf
    Leaf initial_leaf;
    {
        initial_leaf.value = 0;
        initial_leaf.individuals = initial_individuals;
        initial_leaf.intervals = initial_intervals;
    }

    // store possible splits in map with splitting variable as key and pointer to resulting tree
    std::multimap<int, DecisionTree*> possible_splits;

    // initial trees according to variables
    DecisionTree initial_tree;
    TreeFamily initial_family;
    for(size_t i=0; i<variables.size(); ++i){
        initial_tree.split_dims = std::set<int> (variables[i].begin(),variables[i].end());
        initial_tree.leaves.push_back(initial_leaf); // add leaf
        initial_family.push_back(initial_tree); // save tree with one leaf in the beginning
    }

    // initilize tree families
    this->tree_families = std::vector<TreeFamily>(1, initial_family); // todo: change to (n_trees, initial_family);

    // construct some helper objects
    int split_count;
    Split curr_split;
    std::vector<std::vector<double>> samples_X = std::vector<std::vector<double>>(sample_size);
    std::vector<double> samples_Y = std::vector<double>(sample_size); // question: what used for

    // iterate over families of trees and modify
    for(TreeFamily curr_family:tree_families){

        // reset possible splits
        possible_splits.clear();
        for(int i=0; i<variables.size(); ++i){
            for(int j=0; j<variables[i].size(); ++j){
                possible_splits.insert(std::pair<int, DecisionTree*>(variables[i][j], &curr_family[i])); // add pointer to resulting tree with split dimension as key
            }
        }

        // sample data points with replacement
        for(size_t i = 0; i<sample_size; ++i){
            // todo: switch to 'std::uniform_int_distribution' for equally-likely numbers
            int sample_index = std::rand() % sample_size;
            samples_Y[i] = Y[sample_index];
            samples_X[i] = X[sample_index];
        }

        // modify existing or add new trees through splitting
        split_count = 0; // reset number of splits
        while(split_count<n_splits){ // question: while or for loop
            ++split_count; // keep track of number of performed splits
            curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family); // find optimal split

            // continue only if we get a significant result
            if(!std::isinf(curr_split.min_sum)){
                std::cout << curr_split.min_sum << std::endl;
            }

        }

        // todo: clear individuals of each tree
    }
}




