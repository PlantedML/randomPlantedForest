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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;


//  ----------------- functions for converting R and Cpp types ----------------- 

Rcpp::IntegerVector from_std_vec(std::vector<int> v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

Rcpp::NumericVector from_std_vec(std::vector<double> v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

Rcpp::NumericMatrix from_std_vec(std::vector<std::vector<double>> v) {
  if(v.empty()) return NumericMatrix();
  NumericMatrix m(v.size(), v[0].size());
  for(int row=0; row<v.size(); ++row){
    m(row, _ ) = NumericMatrix(1, v[row].size(), from_std_vec(v[row]).begin());
  }
  return m;
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

std::set<int> to_std_set(Rcpp::IntegerVector rv) {
  return std::set<int>(rv.begin(), rv.end());
}


// ----------------- custom data types ----------------- 

typedef std::pair<double, double> Interval;

struct Leaf{
  std::vector<int> individuals;          // considered samples for each leaf
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

namespace rpf{

struct Split {
  double min_sum;                 // minimal achievable sum of squared residuals
  std::shared_ptr<DecisionTree> tree_index;   // pointer to tree
  Leaf* leaf_index;               // pointer to leaf containing interval
  int split_coordinate;           // coordinate for splitting
  double split_point;             // splitpoint
  std::vector<int> I_s;              // individuals smaller than splitpoint
  std::vector<int> I_b;              // individuals bigger than splitpoint
  double I_s_mean;                // mean of individuals smaller than splitpoin
  double I_b_mean;                // mean of individuals bigger than splitpoint
  Split(): min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), I_s_mean(0.0), I_b_mean(0.0) {};
};

}


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


// ----------------- main rpf class -----------------

class RandomPlantedForest {
  
public:
  RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                      int max_interaction=2, int n_trees=50, int n_splits=30, double t_try=0.4);
  void set_data(const NumericVector &samples_Y, const NumericMatrix &samples_X);
  Rcpp::NumericVector predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
  Rcpp::NumericVector predict_vector(const NumericVector &X, const NumericVector components = {0});
  void purify();
  void print();
  void cross_validation(int n_sets, IntegerVector splits={5,50}, NumericVector t_tries={0.2,0.5,0.7,0.9}, IntegerVector split_tries={1,2,5,10});
  double MSE(const NumericVector &Y_predicted, const NumericVector &Y_true); 
  void set_deterministic(bool deterministic);
  void set_parallel(bool parallelize);
  void get_parameters();
  
private:
  std::vector<double> Y;
  std::vector<std::vector<double>> X;
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
  bool deterministic = false;                 // choose whether approach deterministic or random
  bool parallelize = true;                    // 
  std::vector<std::vector<int>> variables;    // split dimensions for initial trees
  std::vector<double> upper_bounds;           //
  std::vector<double> lower_bounds;           //
  std::vector<TreeFamily> tree_families;      // random planted forest conatining result
  void fit();
  void create_tree_family(std::vector<Leaf> initial_leaves, size_t n);
  struct CreateTreeFamilies : public Worker {
    std::vector<Leaf> initial_leaves;
    RandomPlantedForest* forest;
    CreateTreeFamilies(const std::vector<Leaf> initial_leaves, RandomPlantedForest* forest) 
      : initial_leaves(initial_leaves), forest(forest) {}
    void operator()(std::size_t begin, std::size_t end);
  };
  double predict_single(const std::vector<double> &X, std::set<int> component_index);
  rpf::Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                              const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family);
};

// constructor
RandomPlantedForest::RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                                         int max_interaction, int n_trees, int n_splits, double t_try){
  
  // initialize class members
  this->max_interaction = max_interaction;
  this->n_trees = n_trees;
  this->n_splits = n_splits;
  this->split_try = 10;
  this->t_try = t_try;
  this->purify_forest = false;
  
  // set data and data related members
  this->set_data(samples_Y, samples_X);
}

// determine optimal split
rpf::Split RandomPlantedForest::calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                                 const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family){
  rpf::Split curr_split;
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
        std::vector<int> curr_individuals = leaf.individuals; // consider individuals of current leaf
        
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
          
          std::vector<int> I_s, I_b; // individuals smaller/bigger than the samplepoint
          std::vector<double> Y_s, Y_b; // values of individuals smaller/bigger than the samplepoint
          
          // sum of respective samples
          Y_s_sum = 0;
          Y_b_sum = 0;
          
          // get samples greater/smaller than samplepoint
          for(int individual: curr_individuals){
            if(X[individual][k] < sample_point){
              Y_s.push_back(Y[individual]);
              I_s.push_back(individual);
              Y_s_sum += Y[individual];
            }else{
              Y_b.push_back(Y[individual]);
              I_b.push_back(individual);
              Y_b_sum += Y[individual];
            }
          }
          
          // ensure individuals are sorted and unique
          std::sort(I_s.begin(), I_s.end());
          std::sort(I_b.begin(), I_b.end());
          I_s.erase(std::unique(I_s.begin(), I_s.end()), I_s.end());
          I_b.erase(std::unique(I_b.begin(), I_b.end()), I_b.end());
          
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

void RandomPlantedForest::set_data(const NumericVector &samples_Y, const NumericMatrix &samples_X){
  
  this->Y = to_std_vec(samples_Y);
  this->X = to_std_vec(samples_X);
  
  // Check if all vector in x are same length
  for(const auto &vec:X){
    this->feature_size = vec.size();
    if(vec.size() != feature_size) std::cout << "Dimensions of X mismatch." << std::endl;
  }
  if(Y.size() != X.size()){
    std::cout << "Dimensions of X and Y mismatch." << std::endl;
  }
  
  this->n_leaves = std::vector<int>(feature_size, 1);
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
  
  this->variables = std::vector<std::vector<int>>(feature_size);
  for(int i = 0; i<feature_size; ++i) this->variables[i] = std::vector<int> {i+1};
  
  this->fit();
}

void RandomPlantedForest::create_tree_family(std::vector<Leaf> initial_leaves, size_t n){
  
  DecisionTree initial_tree;
  TreeFamily curr_family;
  for(size_t i=0; i<variables.size(); ++i){
    initial_tree.split_dims = std::set<int> (variables[i].begin(),variables[i].end());
    initial_tree.leaves = initial_leaves;
    curr_family.trees.push_back(std::make_shared<DecisionTree>(initial_tree)); // save tree with one leaf in the beginning
  }
  
  // store possible splits in map with splitting variable as key and pointer to resulting tree
  std::multimap<int, std::shared_ptr<DecisionTree>> possible_splits;
  for(int i=0; i<variables.size(); ++i){
    for(int j=0; j<variables[i].size(); ++j){
      // add pointer to resulting tree with split dimension as key
      possible_splits.insert(std::pair<int, std::shared_ptr<DecisionTree>>(variables[i][j], curr_family.trees[i]));
    }
  }
  
  if(false){
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
  std::vector<std::vector<double>> samples_X = std::vector<std::vector<double>>(sample_size);
  std::vector<double> samples_Y = std::vector<double>(sample_size);
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
  rpf::Split curr_split;
  for(size_t split_count=0; split_count<n_splits; ++split_count){
    
    // find optimal split
    curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family);
    
    // continue only if we get a significant result
    if(!std::isinf(curr_split.min_sum)){
      
      if(false){
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
          
          if(false){
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
      
      if(false){
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

void RandomPlantedForest::CreateTreeFamilies::operator()(std::size_t begin, std::size_t end){
  for(size_t n=begin; n<end; ++n){
    forest->create_tree_family(initial_leaves, n);
  }
}

// fit forest to new data
void RandomPlantedForest::fit(){
  
  // setup initial set of individuals
  std::vector<int> initial_individuals(sample_size);
  std::iota(initial_individuals.begin(), initial_individuals.end(), 0);
  
  // initialize intervals with lower and upper bounds
  std::vector<Interval> initial_intervals(feature_size);
  for(size_t i = 0; i<feature_size; ++i) initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};
  
  // set properties of first leaf
  Leaf initial_leaf;
  {
    initial_leaf.value = 0;
    initial_leaf.individuals = initial_individuals;
    initial_leaf.intervals = initial_intervals;
  }
  std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf
  
  // initialize tree families
  this->tree_families = std::vector<TreeFamily>(n_trees); 
  
  // iterate over families of trees and modify
  if(parallelize){
    std::cout << "Note: Using parallelized version." << std::endl;
    CreateTreeFamilies create_tree_families(initial_leaves, this);
    parallelFor(0, n_trees, create_tree_families);
  }else{
    for(size_t n=0; n<n_trees; ++n){
      create_tree_family(initial_leaves, n);
    }
  }
}

void RandomPlantedForest::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries){
  
  if(deterministic) {
    std::cout << "Note: Set model to non-deterministic. " << std::endl;
    deterministic = false; 
  }
  
  std::set<int> splits_vec = to_std_set(splits);
  std::vector<int> split_tries_vec = to_std_vec(split_tries);
  std::vector<double> t_tries_vec = to_std_vec(t_tries);
  
  if(splits_vec.size()!=2) {std::cout << "Min and max needed for number of splits." << std::endl; return;}
  
  // remember optimal parameter set and MSE
  double  MSE_sum = 0, curr_MSE = 0, MSE_min = INF, optimal_split = INF, optimal_t_try = INF, optimal_split_try = INF;
  int optimal_inter = 1;
  
  std::vector<int> order(sample_size);
  std::iota(order.begin(), order.end(), 0);
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(order.begin(), order.end(), g);
  double tmp = double(sample_size)/double(n_sets);
  int set_size = round(tmp);
  
  // remember original data samples
  NumericMatrix X_original = from_std_vec(X);
  NumericVector Y_original = from_std_vec(Y);
  
  // set level of interactions
  std::set<int> interactions{1};
  if(feature_size >= 2){
    interactions.insert(2);
    interactions.insert(feature_size);
  }
  
  // go through all parameter combinations
  for(int inter: interactions){
    this->max_interaction = inter;
    for(int splits=*splits_vec.begin(); splits<=*--splits_vec.end(); splits=ceil(splits*1.2)){
      this->n_splits = splits;
      for(auto t: t_tries){
        this->t_try = t;
        for(auto s: split_tries){	
          this->split_try = s;
          
          // k-fold cross-validation: go over all possible combinations as test set
          MSE_sum = 0;
          for(int n_set=0; n_set<n_sets; ++n_set){
            
            // split data into training and test sets
            int test_size = set_size;
            if(n_set == n_sets-1) test_size = order.size() - (n_sets-1) * set_size;
            int train_size = order.size() - test_size, i = 0, j = 0;
            NumericVector Y_train(train_size), Y_test_true(test_size), Y_test_predicted;
            NumericMatrix X_train(train_size, feature_size), X_test(test_size, feature_size);
            for(int index=0; index<order.size(); ++index){
              if( (index >= (n_set * set_size)) && (index < ((n_set + 1) * set_size))){
                Y_test_true[i] = Y_original[order[index]];
                X_test(i, _ ) = X_original(order[index], _ );
                ++i;
              }else{
                Y_train[j] = Y_original[order[index]];
                X_train(j, _ ) = X_original(order[index], _ );
                ++j;
              }
            }
            
            // fit to training data
            this->set_data(Y_train, X_train);
            
            // predict with test set and determine mse
            Y_test_predicted = this->predict_matrix(X_test); 
            MSE_sum += this->MSE(Y_test_predicted, Y_test_true);
          }
          
          // average
          curr_MSE = MSE_sum / n_sets;
          std::cout << inter << ", " << splits << ", " << t << ", " << s << ": MSE=" << curr_MSE << std::endl;
          
          // update optimal
          if(curr_MSE < MSE_min){
            MSE_min = curr_MSE;
            optimal_split = splits;
            optimal_t_try = t;
            optimal_split_try = s;
            optimal_inter = inter;
          }
        }
      }	
    }
  }
  
  // reset X&Y to original and fit with optimal pars
  this->n_splits = optimal_split;
  this->t_try = optimal_t_try;
  this->split_try = optimal_split_try;
  this->max_interaction = optimal_inter;
  this->set_data(Y_original, X_original);
  
  std::cout << "Optimal parameters: " << optimal_inter << ", " << optimal_split << ", " << optimal_t_try << ", " << optimal_split_try << ": MSE=" << MSE_min << std::endl;
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

double RandomPlantedForest::MSE(const NumericVector &Y_predicted, const NumericVector &Y_true){
  return sum(Rcpp::pow(Y_true - Y_predicted, 2)) / Y_true.size();
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

void RandomPlantedForest::get_parameters(){
  std::cout << "Parameters: n_trees=" <<  n_trees << ", n_splits=" << n_splits << ", max_interaction=" << max_interaction << ", t_try=" << t_try 
            << ", split_try=" << split_try << ", purified=" << purified << ", deterministic=" << deterministic << ", feature_size=" << feature_size 
            << ", sample_size=" << sample_size << std::endl;
}

void RandomPlantedForest::set_deterministic(bool deterministic){
  this->deterministic = deterministic;
  this->fit();
}

void RandomPlantedForest::set_parallel(bool parallelize){
  this->parallelize = parallelize;
}


// ----------------- Rcpp include  -----------------

RCPP_MODULE(mod_rpf) {
  
  class_<RandomPlantedForest>("RandomPlantedForest")
  .constructor<const NumericVector, const NumericMatrix, int, int, int, double>()
  .method("set_data", &RandomPlantedForest::set_data)
  .method("cross_validation", &RandomPlantedForest::cross_validation)
  .method("predict_matrix", &RandomPlantedForest::predict_matrix)
  .method("predict_vector", &RandomPlantedForest::predict_vector)
  .method("MSE", &RandomPlantedForest::MSE)
  .method("purify", &RandomPlantedForest::purify)
  .method("print", &RandomPlantedForest::print)
  .method("get_parameters", &RandomPlantedForest::get_parameters)
  .method("set_deterministic", &RandomPlantedForest::set_deterministic)
  .method("set_parallel", &RandomPlantedForest::set_parallel)
  ;
  
}