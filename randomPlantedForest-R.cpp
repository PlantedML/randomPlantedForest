#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>
#include <optional>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <Rcpp.h>

using namespace Rcpp;


//  ----------------- functions for converting R and Cpp types ----------------- 

Rcpp::IntegerVector from_std_vec(std::vector<int> v) {
    return Rcpp::IntegerVector(v.begin(), v.end());
}

Rcpp::NumericVector from_std_vec(std::vector<double> v) {
    return Rcpp::NumericVector(v.begin(), v.end());
}

std::vector<int> to_std_vec(Rcpp::IntegerVector rv) {
    return std::vector<int>(rv.begin(), rv.end());
}

std::vector<double> to_std_vec(Rcpp::NumericVector rv) {
    return std::vector<double>(rv.begin(), rv.end());
}

std::vector<std::vector<double>> to_std_vec(Rcpp::NumericMatrix rv) {
    std::vector<std::vector<double>> X;
    for(int i=0; i<rv.rows(); i++) X.push_back(to_std_vec(rv(i, _ )));
    return X;
}

std::set<int> to_std_set(Rcpp::NumericVector rv) {
    return std::set<int>(rv.begin(), rv.end());
}


// ----------------- custom data types ----------------- 

typedef std::pair<double, double> Interval;

struct Leaf{
    std::set<int> individuals;          // considered samples for each leaf
    double value;                       // residual
    std::vector<Interval> intervals;    // min/max for each feature of the interval
};

class DecisionTree {
    
    friend class RandomPlantedForest;
    
    public:
        DecisionTree() {};
        DecisionTree(std::set<int> dims, std::vector<Leaf> first_leaves):
            split_dims(dims), leaves(first_leaves) {};
        DecisionTree(std::set<int> dims): split_dims(dims) {};
        std::set<int> get_split_dims() const;
        std::vector<Leaf> get_leaves() const;
        
    private:
        std::set<int> split_dims;       // dimensions of the performed splits
        std::vector<Leaf> leaves;       // leaves of tree containing intervals and approximating value
        // idea: save intervals as interval-tree with nodes and corresponding values
};

std::set<int> DecisionTree::get_split_dims() const{
    return split_dims;
}

std::vector<Leaf> DecisionTree::get_leaves() const{
    return leaves;
}

struct TreeFamily {
    std::vector<std::shared_ptr<DecisionTree>> trees;
    double constant = 0;
};

const double INF = std::numeric_limits<double>::infinity();

struct Split {
    double min_sum;                 // minimal achievable sum of squared residuals
    std::shared_ptr<DecisionTree> tree_index;   // pointer to tree
    Leaf* leaf_index;               // pointer to leaf containing interval
    int split_coordinate;           // coordinate for splitting
    double split_point;             // splitpoint
    std::set<int> I_s;              // individuals smaller than splitpoint
    std::set<int> I_b;              // individuals bigger than splitpoint
    double I_s_mean;                // mean of individuals smaller than splitpoin
    double I_b_mean;                // mean of individuals bigger than splitpoint
    Split(): min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), I_s_mean(0.0), I_b_mean(0.0) {};
};


// ----------------- helper functions -----------------

// helper function to check whether a tree with specified split_dims already exists in tree_family
std::shared_ptr<DecisionTree> treeExists(const std::set<int> split_dims, TreeFamily &tree_family){
    for(auto& tree: tree_family.trees){
        if(tree->get_split_dims() == split_dims){
            return tree; // if found, return pointer to tree, otherwise nullptr
        }
    }
    return nullptr;
}

// helper function to check whether a tree with resulting_dims for a split_coordinate is already in possible_splits
bool possibleExists(const int dim, const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, const std::set<int> &resulting_dims){
    for(auto& elem:possible_splits){
        if(elem.first == dim && elem.second->get_split_dims() == resulting_dims) return 1;
    }
    return 0;
}

bool leafExists(std::vector<Interval>& intervals, const std::shared_ptr<DecisionTree> tree){
    bool exists = false;
    for(auto& leaf: tree->get_leaves()){
        
        bool same_intervals = true;
        for(auto& curr_interval: leaf.intervals){
            for(auto& interval: intervals){
                if(curr_interval != interval){
                    same_intervals = false;
                    break;
                }
            }
        }
        
        if(same_intervals){
            exists = true;
        }
    }
    return exists;
}

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> vec){
    if(vec.size() == 0) return vec;
    std::vector<std::vector<double>> transposed(vec[0].size(),std::vector<double>());
    
    for(int i=0; i<vec.size(); ++i){
        for(int j=0; j<vec[i].size(); ++j){
            transposed[j].push_back(vec[i][j]);
        }
    }
    
    for(auto col: transposed){
        for(auto el: col){
            std::cout << el << ", ";
        }
        std::cout << std::endl;
    }
    
    return transposed;
}

// ----------------- main rpf class -----------------

class RandomPlantedForest {

    public:
        RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                            int max_interaction=2, int n_trees=50, int n_splits=30, double t_try=0.4);
        void fit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X);
        Rcpp::NumericVector predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
        Rcpp::NumericVector predict_vector(const NumericVector &X, const NumericVector components = {0});
        void purify();
        // todo: getter/setter
        void print();

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
        bool deterministic = true;                  // choose whether approach deterministic or random
        std::vector<std::vector<int>> variables;    // split dimensions for initial trees
        std::vector<double> upper_bounds;           //
        std::vector<double> lower_bounds;           //
        std::vector<TreeFamily> tree_families;      // random planted forest conatining result
        double predict_single(const std::vector<double> &X, std::set<int> component_index);
        Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                               const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family);
};

// constructor
RandomPlantedForest::RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                                         int max_interaction, int n_trees, int n_splits, double t_try){

    std::vector<double> Y = to_std_vec(samples_Y);
    std::vector<std::vector<double>> X = to_std_vec(samples_X);

    // Check if all vector in x are same length
    for(const auto &vec:X){
        this->feature_size = vec.size();
        if(vec.size() != feature_size) std::cout << "Dimensions of X mismatch." << std::endl;
    }

    // initialize class members
    this->max_interaction = max_interaction;
    this->n_trees = n_trees;
    this->n_splits = n_splits;
    this->split_try = 10;
    this->t_try = t_try;
    this->purify_forest = false;
    this->n_leaves = std::vector<int>(feature_size, 1);
    this->variables = std::vector<std::vector<int>>(feature_size);
    for(int i = 0; i<feature_size; ++i) this->variables[i] = std::vector<int> {i+1};

    // construct tree families
    this->fit(Y, X);
}

// determine optimal split
Split RandomPlantedForest::calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                            const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family){
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
    
    if(!deterministic){
        std::shuffle(split_candidates.begin(), split_candidates.end(), gen); // shuffle for random order
    }
    
    if(true){
        std::cout << "Current candidates: (" << n_candidates << ") ";
        for(size_t n=0; n<n_candidates; ++n){
            auto candidate = possible_splits.begin();
            std::advance(candidate, split_candidates[n]);
            std::cout << candidate->first << "- ";
            for(auto dim: candidate->second->split_dims) std::cout << dim << ",";
            std::cout << "; ";
        }
        std::cout << std::endl;
    }

    // go through all trees in current family
    for(auto& curr_tree:curr_family.trees){
        
        // skip if tree has no leaves
        if(curr_tree->leaves.size() == 0){ continue; }

        // consider a fraction of possible splits
        for(size_t n=0; n<n_candidates; ++n){
            
            auto candidate = possible_splits.begin();
            std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
            
            k = candidate->first - 1; // split dim of current candidate, converted to index starting at 0
            leaf_size = n_leaves[k];
            
            // Test if splitting in the current tree w.r.t. the coordinate "k" is an element of candidate tree
            tree_dims = curr_tree->split_dims;
            tree_dims.insert(k+1);
            if(tree_dims != candidate->second->split_dims) continue; // if split dimensions do not match consider next candidate
            
            // go through all leaves of current tree
            for(auto& leaf: curr_tree->leaves){
                std::set<int> curr_individuals = leaf.individuals; // consider individuals of current leaf
                
                // extract sample points according to individuals from X and Y
                std::map<double, double> unique_samples;
                for(int individual: curr_individuals){
                    unique_samples[X[individual][k]] = Y[individual];
                }
                
                // check if number of sample points is within limit
                if(unique_samples.size() < 2*leaf_size) continue;
                
                int start = 0;
                int end = split_try;
                if(deterministic){
                    start = 1;
                    end = unique_samples.size()-1;
                }
                
                // consider split_try-number of random samples
                for(int t = start; t<end; ++t){
                    
                    // get samplepoint
                    auto sample_pos = unique_samples.begin();
                    std::uniform_int_distribution<> distrib(leaf_size, unique_samples.size() - leaf_size + 1);
                    if(deterministic){
                        std::advance(sample_pos, t);
                    }else{
                        std::advance(sample_pos, distrib(gen)); // consider only sample points with offset
                    }
                    sample_point = sample_pos->first;
                    
                    std::set<int> I_s, I_b; // individuals smaller/bigger than the samplepoint
                    std::vector<double> Y_s, Y_b; // values of individuals smaller/bigger than the samplepoint
                    
                    // sum of respective samples
                    Y_s_sum = 0;
                    Y_b_sum = 0;
                    
                    // get samples greater/smaller than samplepoint
                    for(int individual: curr_individuals){
                        if(X[individual][k] < sample_point){
                            Y_s.push_back(Y[individual]);
                            I_s.insert(individual);
                            Y_s_sum += Y[individual];
                        }else{
                            Y_b.push_back(Y[individual]);
                            I_b.insert(individual);
                            Y_b_sum += Y[individual];
                        }
                    }
                    
                    // get mean
                    if(Y_s.size() != 0) Y_s_mean = Y_s_sum / Y_s.size();
                    if(Y_b.size() != 0) Y_b_mean = Y_b_sum / Y_b.size();
                    
                    // accumulate squared mean
                    curr_sum = 0;
                    std::for_each(Y_s.begin(), Y_s.end(), [&Y_s_mean, &curr_sum](double val){ curr_sum += pow(val - Y_s_mean, 2) - pow(val, 2); });
                    std::for_each(Y_b.begin(), Y_b.end(), [&Y_b_mean, &curr_sum](double val){ curr_sum += pow(val - Y_b_mean, 2) - pow(val, 2); });
                    
                    if(false){
                        std::cout << "Current Min=" << curr_sum << "; ";
                        std::cout << "Means=" <<  Y_s_mean << "/" <<  Y_b_mean << std::endl;
                    }
                    
                    // update split if squared sum is smaller
                    if(curr_sum < curr_split.min_sum){
                        curr_split.min_sum = curr_sum;
                        curr_split.tree_index = curr_tree;
                        curr_split.leaf_index =  &leaf;
                        curr_split.split_coordinate = k+1;
                        curr_split.split_point = sample_point;
                        curr_split.I_s_mean = Y_s_mean;
                        curr_split.I_b_mean = Y_b_mean;
                        curr_split.I_s = I_s;
                        curr_split.I_b = I_b;
                    }
                }
            }
        }
    }
    
    return curr_split;
}

// fit forest to new data
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
    
    if(false){
        std::cout << "Upper bounds: ";
        for(auto val: upper_bounds) std::cout << val << ", ";
        std::cout << "Lower bounds: ";
        for(auto val: lower_bounds) std::cout << val << ", ";
        std::cout << std::endl;
    }
    
    std::set<int> feature_dims;
    auto pos = feature_dims.begin();
    for(int i = 1; i <= feature_size; ++i) pos = feature_dims.insert(pos, i);
    
    // setup initial set of individuals
    std::set<int> initial_individuals;
    pos = initial_individuals.begin();
    for(int i = 0; i < sample_size; ++i) pos = initial_individuals.insert(pos, i);
    
    if(false){
        std::cout << "Initial individuals: (" << sample_size << ") ";
        for(auto val: initial_individuals) std::cout << val << ", ";
        std::cout << std::endl;
    }
    
    // initialize intervals with lower and upper bounds
    std::vector<Interval> initial_intervals(feature_size);
    for(size_t i = 0; i<feature_size; ++i) initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};
    
    if(false){
        std::cout << "Initial intervals: ";
        for(auto interval: initial_intervals) std::cout << interval.first << ", " << interval.second << "; ";
        std::cout << std::endl;
    }
    
    // set properties of first leaf
    Leaf initial_leaf;
    {
        initial_leaf.value = 0;
        initial_leaf.individuals = initial_individuals;
        initial_leaf.intervals = initial_intervals;
    }
    std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf
    
    // store possible splits in map with splitting variable as key and pointer to resulting tree
    std::multimap<int, std::shared_ptr<DecisionTree>> possible_splits;
    
    // construct some helper objects
    Split curr_split;
    std::vector<std::vector<double>> samples_X = std::vector<std::vector<double>>(sample_size);
    std::vector<double> samples_Y = std::vector<double>(sample_size);
    
    // initialize tree families
    this->tree_families = std::vector<TreeFamily>(n_trees); 
    
    // iterate over families of trees and modify
    for(size_t n=0; n<n_trees; ++n){
        
        DecisionTree initial_tree;
        TreeFamily curr_family;
        for(size_t i=0; i<variables.size(); ++i){
            initial_tree.split_dims = std::set<int> (variables[i].begin(),variables[i].end());
            initial_tree.leaves = initial_leaves;
            curr_family.trees.push_back(std::make_shared<DecisionTree>(initial_tree)); // save tree with one leaf in the beginning
        }
        
        if(true){
            std::cout << "Initial TreeFamily: (" << curr_family.trees.size() << ") ";
            for(auto tree: curr_family.trees){
                std::cout << "Dims = ";
                for(auto dim: tree->split_dims) std::cout << dim << ", ";
                std::cout << "; " << "Number of Leafs = " << tree->leaves.size();
                std::cout << " / ";
            }
            std::cout << std::endl;
        }
        
        // reset possible splits
        possible_splits.clear();
        for(int i=0; i<variables.size(); ++i){
            for(int j=0; j<variables[i].size(); ++j){
                // add pointer to resulting tree with split dimension as key
                possible_splits.insert(std::pair<int, std::shared_ptr<DecisionTree>>(variables[i][j], curr_family.trees[i]));
            }
        }
        
        if(true){
            std::cout << "Initial Possible Splits: ";
            for(auto split: possible_splits){
                std::cout << split.first << "-";
                for(auto dim: split.second->split_dims) std::cout << dim << ",";
                std::cout << "; ";
            }
            std::cout << std::endl;
        }
        
        // sample data points with replacement
        int sample_index;
        for(size_t i=0; i<sample_size; ++i){
            // todo: switch to 'std::uniform_int_distribution' for equally-likely numbers
            sample_index = std::rand() % sample_size;
            samples_Y[i] = Y[sample_index];
            samples_X[i] = X[sample_index];
        }
        
        // deterministic
        if(deterministic){
            samples_X = X;
            samples_Y = Y;
            this->t_try = 1;
        }
        
        // modify existing or add new trees through splitting
        for(size_t split_count=0; split_count<n_splits; ++split_count){
            
            // find optimal split
            curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family);
            
            // continue only if we get a significant result
            if(!std::isinf(curr_split.min_sum)){
                
                if(true){
                    std::cout << "Current Optimal Split: " << curr_split.min_sum << "; " << curr_split.split_coordinate << "- ";
                    for(auto dim: curr_split.tree_index->split_dims) std::cout << dim << ", ";
                    std::cout << "; " << curr_split.I_s.size() << "/" << curr_split.I_b.size() << "=" << curr_split.I_s.size()+curr_split.I_b.size() << "; " <<
                        curr_split.I_s_mean << "/" << curr_split.I_b_mean <<  std::endl;
                }

                // update possible splits
                if(curr_split.tree_index->split_dims.count(curr_split.split_coordinate) == 0 || curr_split.tree_index->leaves.size() == 1){

                    // consider all possible dimensions
                    for(int feature_dim = 1; feature_dim<=feature_size; ++feature_dim){
                        
                        // ignore dim if same as split coordinate or in dimensions of old tree
                        if(feature_dim == curr_split.split_coordinate || curr_split.tree_index->split_dims.count(feature_dim) > 0) continue;
                        
                        // create union of split coord, feature dim and dimensions of old tree
                        std::set<int> curr_dims = curr_split.tree_index->split_dims;
                        curr_dims.insert(curr_split.split_coordinate);
                        curr_dims.insert(feature_dim);
                        
                        if(curr_dims.size()>max_interaction) continue; 
                        
                        // check if resulting tree already exists in family
                        std::shared_ptr<DecisionTree> found_tree = treeExists(curr_dims, curr_family);
                        
                        // update possible_splits if not already existing
                        if(found_tree){ // if yes add pointer
                            possible_splits.insert(std::pair<int, std::shared_ptr<DecisionTree>>(feature_dim, found_tree));
                        }else{ // if not create new tree
                            curr_family.trees.push_back(std::make_shared<DecisionTree>(DecisionTree(curr_dims)));
                            possible_splits.insert(std::pair<int, std::shared_ptr<DecisionTree>>(feature_dim, curr_family.trees.back()));
                        }
                        
                        if(true){
                            std::cout << "Updated Possible Splits: " << std::endl;
                            for(auto split: possible_splits){
                                std::cout << split.first << "-";
                                for(auto dim: split.second->split_dims) std::cout << dim << ",";
                                std::cout << "; ";
                            }
                            std::cout << std::endl;
                        }
                    }
                }
                
                // update values of individuals of split interval with mean
                for(int individual: curr_split.leaf_index->individuals){ // todo: loop directly over I_s I_b
                    if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
                        samples_Y[individual] -= curr_split.I_s_mean;
                    }else{
                        samples_Y[individual] -= curr_split.I_b_mean;
                    }
                }
                
                // construct new leaves
                Leaf leaf_s, leaf_b;
                {
                    leaf_s.individuals = curr_split.I_s;
                    leaf_b.individuals = curr_split.I_b;
                    leaf_s.value = curr_split.I_s_mean;
                    leaf_b.value = curr_split.I_b_mean;
                    
                    // initialize interval with split interval
                    leaf_s.intervals = curr_split.leaf_index->intervals;
                    leaf_b.intervals = curr_split.leaf_index->intervals;
                    
                    // interval of leaf with smaller individuals has new upper bound in splitting dimension
                    leaf_s.intervals[curr_split.split_coordinate-1].second = curr_split.split_point;
                    // interval of leaf with bigger individuals has new lower bound in splitting dimension
                    leaf_b.intervals[curr_split.split_coordinate-1].first = curr_split.split_point;
                }
                
                if(false){
                    std::cout << "First leaf: intervals=" << std::endl;
                    for(auto interval: leaf_s.intervals) std::cout << interval.first << "," << interval.second << ";";
                    std::cout << "individuals=";
                    for(auto i: leaf_s.individuals) std::cout << i << ",";
                    std::cout << std::endl << "Second leaf: intervals=" << std::endl;
                    for(auto interval: leaf_b.intervals) std::cout << interval.first << "," << interval.second << ";";
                    std::cout << "individuals=";
                    for(auto i: leaf_b.individuals) std::cout << i << ",";
                    std::cout << std::endl;
                }
                
                // construct split_dims of resulting tree when splitting in split_coordninate
                std::set<int> resulting_dims = curr_split.tree_index->split_dims;
                resulting_dims.insert(curr_split.split_coordinate);
                
                // check if resulting tree already exists in family
                std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);
                
                // determine which tree is modified
                if(curr_split.tree_index->split_dims.count(curr_split.split_coordinate) ){ // if split variable is already in tree to be split
                    // change values
                    {
                        leaf_s.value += curr_split.leaf_index->value;
                        leaf_b.value += curr_split.leaf_index->value;
                    }
                    *curr_split.leaf_index = leaf_b; // replace old interval
                    curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
                } else{ // otherwise 
                    found_tree->leaves.push_back(leaf_s); //append new leaves
                    found_tree->leaves.push_back(leaf_b);
                }
                
                if(true){
                    std::cout << "Current TreeFamily: (" << curr_family.trees.size() << ") ";
                    for(auto tree: curr_family.trees){
                        std::cout << "Dims = ";
                        for(auto dim: tree->split_dims) std::cout << dim << ", ";
                        std::cout << "; " << "Number of Leafs = " << tree->leaves.size();
                        std::cout << " / ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }
        
        // remove empty trees
        curr_family.trees.erase(std::remove_if(curr_family.trees.begin(), curr_family.trees.end(), 
                                                [](std::shared_ptr<DecisionTree> tree) { return tree->leaves.size() == 0;}), curr_family.trees.end() );
        
        // todo: clear individuals of each tree
        
        // optional: purify tree
        if(purify_forest){
            this->purify();
        }else{
            purified = false;
        }
        
        tree_families[n] = curr_family;
    }
}

// predict single feature vector
double RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index){

    double total_res = 0;
    
    // consider all components
    if(component_index == std::set<int>{0}) {
        for(auto& tree_family: this->tree_families){
            for(auto& tree: tree_family.trees){
                for(auto& leaf: tree->leaves){
                    bool valid = true;
                    for(auto& dim: tree->split_dims){
                        if(!(leaf.intervals[dim-1].first <= X[dim-1]
                            && (leaf.intervals[dim-1].second > X[dim-1]
                            || leaf.intervals[dim-1].second == upper_bounds[dim-1]))){
                            valid = false;
                        }
                    }
                    if(valid) total_res += leaf.value;
                }
            }
            total_res += tree_family.constant;
        }
    }else{ // choose components for prediction
        for(auto& tree_family: this->tree_families){
            for(auto& tree: tree_family.trees){
                
                // only consider trees with same dimensions as component_index
                if(tree->split_dims != component_index) continue;
                
                std::vector<int> dims;
                for(auto dim: tree->split_dims) {dims.push_back(dim);}
                
                for(auto& leaf: tree->leaves){
                    bool valid = true;
                    for(int i = 0; i<dims.size(); ++i){

                        int dim = dims[i];
                        std::cout << X[i] << ", " << leaf.intervals[dim-1].first << std::endl;
                        
                        if(!(leaf.intervals[dim-1].first <= X[i]
                            && (leaf.intervals[dim-1].second > X[i]
                            || leaf.intervals[dim-1].second == upper_bounds[dim-1]))){
                            valid = false;
                        }
                    }
                    if(valid) total_res += leaf.value;
                }
            }
        }
    }
    
    double average_res = total_res / n_trees;
    
    return average_res;
}

// predict multiple feature vectors
Rcpp::NumericVector RandomPlantedForest::predict_matrix(const NumericMatrix &X, const NumericVector components){
    std::vector<std::vector<double>> feature_vec = to_std_vec(X);
    std::set<int> component_index = to_std_set(components);
    std::vector<double> predictions;
    Rcpp::NumericVector res;
    
    // todo: sanity check for X
    
    if(!this->purified) { std::cout << "Note: The estimator has not been purified." << std::endl;}
    if(feature_vec.empty()) { std::cout << "Feature vector is empty." << std::endl; return res; }
    
    if(component_index == std::set<int>{0} && feature_vec[0].size() != this->feature_size){ 
        std::cout << "Feature vector has wrong dimension." << std::endl; return res;
    }
    if(component_index != std::set<int>{0} && component_index.size() != feature_vec[0].size()){
        std::cout << "The input X has the wrong dimension in order to calculate f_i(x)" << std::endl; return res;
    }  
    
    for(auto& vec: feature_vec){
        predictions.push_back(predict_single(vec, component_index));
    }

    res = from_std_vec(predictions);
    return res;
}

Rcpp::NumericVector RandomPlantedForest::predict_vector(const NumericVector &X, const NumericVector components){
    std::vector<double> feature_vec = to_std_vec(X);
    std::set<int> component_index = to_std_set(components);
    std::vector<double> predictions;
    Rcpp::NumericVector res;
    
    // todo: sanity check for X
    
    if(!this->purified) { std::cout << "Note: The estimator has not been purified." << std::endl;}
    if(feature_vec.empty()) { std::cout << "Feature vector is empty." << std::endl; return res; }
    
    if(component_index == std::set<int>{0} && feature_vec.size() != this->feature_size){ 
        std::cout << "Feature vector has wrong dimension." << std::endl; return res;
    }
    
    if(component_index == std::set<int>{0}) {
        predictions.push_back( predict_single(feature_vec, component_index));
    }else{
        for(auto vec: feature_vec){
            predictions.push_back( predict_single(std::vector<double>{vec}, component_index));
        }
    }
    
    
    res = from_std_vec(predictions);
    return res;
}

void RandomPlantedForest::purify(){
    
    // go through all n_trees families 
    for(auto& curr_family: this->tree_families){
        
        // recap maximum number of dimensions of current family
        unsigned int curr_max = 0;
        for(auto& tree: curr_family.trees){
            if(tree->split_dims.size() > curr_max) curr_max = tree->split_dims.size();
        }

        while(curr_max >= 1){
            
            if(false){
                std::cout << curr_max << ", "; 
            }
            
            int numb_of_trees = curr_family.trees.size();
            // go through split dimensions of all trees
            for(int n=0; n<numb_of_trees; ++n){
                auto& curr_tree = curr_family.trees[n];
                
                std::set<int> curr_dims = curr_tree->split_dims;
                
                // check if number of dims same as current max_interaction
                if(curr_dims.size() == curr_max){
                    
                    if(curr_max == 1){
                        
                        // new leaf including intervals and value
                        Leaf new_leaf = curr_tree->leaves[0]; // initialize intervals with first leaf 
                        new_leaf.value = 0; // set value to zero
                        int curr_dim = *curr_dims.begin() - 1; // since maximal 1 split dimension, converted to index starting at 0
                        
                        // replace intervals with min/max at current split dimensions
                        new_leaf.intervals[curr_dim].first = lower_bounds[curr_dim];
                        new_leaf.intervals[curr_dim].second = upper_bounds[curr_dim]; 
                        
                        // go through all leaves
                        for(auto& leaf: curr_tree->leaves){
                            
                            //
                            double multiplier = (leaf.intervals[curr_dim].second - leaf.intervals[curr_dim].first) / (upper_bounds[curr_dim] - lower_bounds[curr_dim]);
                            
                            // update value of new leaf
                            new_leaf.value -= leaf.value * multiplier;
                            
                            // update constant of family
                            curr_family.constant += leaf.value * multiplier;
                        }
                        
                        // append new leaf
                        curr_tree->leaves.push_back(new_leaf);
                        
                    }else{
                        
                        // go through feature dims
                        for(int feature_dim=1; feature_dim<=feature_size; ++feature_dim){
                            
                            // continue only if dim in current tree
                            if(curr_tree->split_dims.count(feature_dim) != 0){
                                
                                //
                                std::set<int> tree_dims = curr_tree->split_dims;
                                tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree
                                
                                // check if tree with dimensions exists, if not create
                                std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
                                if(!tree) curr_family.trees.push_back(std::make_shared<DecisionTree>(DecisionTree(tree_dims)));
                                
                                // go through leafs of current tree
                                int n_leafs = curr_tree->leaves.size();
                                for(int l=0; l<n_leafs; ++l){
                                    auto& curr_leaf = curr_tree->leaves[l];
                                    
                                    double multiplier = (curr_leaf.intervals[feature_dim-1].second - curr_leaf.intervals[feature_dim-1].first) 
                                                        / (upper_bounds[feature_dim-1] - lower_bounds[feature_dim-1]);
                                    
                                    // new leaf including intervals and value
                                    Leaf new_leaf = curr_leaf; // initialize intervals with first leaf 
                                    new_leaf.intervals[feature_dim-1].first = lower_bounds[feature_dim-1];
                                    new_leaf.intervals[feature_dim-1].second = upper_bounds[feature_dim-1]; 
                                    new_leaf.value = -curr_leaf.value * multiplier; // update value of new leaf
                                    
                                    // append new leaf
                                    if(!leafExists(new_leaf.intervals, curr_tree)) curr_tree->leaves.push_back(new_leaf);
                                    new_leaf.value = curr_leaf.value * multiplier; // update value of new leaf
                                    if(!leafExists(new_leaf.intervals, tree)) tree->leaves.push_back(new_leaf);
                                }
                            }
                        }
                    }
                }
            }
            
            // update currently considered dimension size
            --curr_max;
        }
        
        std::cout  << std::endl;
    }
    
    purified = true;
}

void RandomPlantedForest::print(){
    for(int n=0; n<n_trees; ++n){
        TreeFamily family = tree_families[n];
        std::cout << n+1 << " TreeFamily: constant=" << family.constant << std::endl << std::endl;
        for(int m=0; m<family.trees.size(); ++m){
            DecisionTree tree = *(family.trees[m]);
            std::cout << m+1 << " Tree: ";
            std::cout << "Dims=";
            for(auto& dim: tree.split_dims) std::cout << dim << ",";
            std::cout << std::endl << "Leafs: (" << tree.leaves.size() << ")" << std::endl;
            for(auto& leaf: tree.leaves){
                std::cout << "Intervals=";
                for(auto& interval: leaf.intervals){
                    std::cout << interval.first << "," << interval.second << "/";
                }
                std::cout << " Value=" << leaf.value << std::endl;
            }  
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}


// ----------------- Rcpp include  -----------------

RCPP_MODULE(mod_rpf) {

    class_<RandomPlantedForest>("RandomPlantedForest")
    .constructor<const NumericVector, const NumericMatrix, int, int, int, double>()
    .method("predict_matrix", &RandomPlantedForest::predict_matrix)
    .method("predict_vector", &RandomPlantedForest::predict_vector)
    .method("purify", &RandomPlantedForest::purify)
    .method("print", &RandomPlantedForest::print)
    ;

}
