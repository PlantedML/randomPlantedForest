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
#include <thread>
#include <cassert>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;


//  ----------------- functions for converting R and Cpp types ----------------- 

/**
 * \brief Convert the std container set of type int into an IntegerVector 
 * from rcpp.
 * 
 * \param v the vector that is converted.
 */
Rcpp::IntegerVector from_std_set(std::set<int> v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

/**
 * \brief Convert the std container vector of type int into an IntegerVector 
 * from rcpp.
 * 
 * \param v the vector that is converted.
 */
Rcpp::IntegerVector from_std_vec(std::vector<int> v) {
    return Rcpp::IntegerVector(v.begin(), v.end());
}

/**
 * \brief Convert the std container vector of type double into a NumericVector 
 * from rcpp.
 * 
 * \param v the vector that is converted.
 */
Rcpp::NumericVector from_std_vec(std::vector<double> v) {
    return Rcpp::NumericVector(v.begin(), v.end());
}

/**
 * \brief Convert the nested std container vector containing a vector itself
 * of type double into a NumericMatrix from rcpp.
 * 
 * Predefines a NumericMatrix of respective size. Afterwards iterates over the 
 * outer vector, then transforms the inner vector using 'from_std_vec'-function 
 * and inserts row into NumericMatrix at the correct position.
 * 
 * \param v the vector of vectors that is converted.
 */
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

struct setComp{
  bool operator()(const std::set<int> &a, const std::set<int> &b){
    if(a == b) return false; // what if same?
    if(a.size() == b.size()){
      std::set<int>::iterator it2 = b.begin();
      for(std::set<int>::iterator it1=a.begin();it1!=a.end();it1++){
        if(*it1 < *it2){
          return true; 
        }else if(*it1 > *it2){
          return false;
        } 
        it2++;
      }
    }
    return (a.size() < b.size());
  }
};

typedef std::pair<double, double> Interval;

/**
 * \brief Structure to remember split correlated data like the intervals 
 * for each dimension, the individuals that are contained and the corresponding 
 * value per split.
 */
struct Leaf{
    std::vector<int> individuals;       /**< considered samples for each leaf */
    double value;                       /**< residual */
    std::vector<Interval> intervals;    /**< min/max for each feature of the interval */
};

/**
 * \brief Decision trees contain split data as leaves for respective splitting 
 * dimensions.
 */
class DecisionTree {
    
    friend class RandomPlantedForest;
    friend class ClassificationRPF;
    
    public:
        DecisionTree() {};
        DecisionTree(std::set<int> dims, std::vector<Leaf> first_leaves):
            split_dims(dims), leaves(first_leaves) {};
        DecisionTree(std::set<int> dims): split_dims(dims) {};
        std::set<int> get_split_dims() const;
        std::vector<Leaf> get_leaves() const;
        
    private:
        std::set<int> split_dims;       /**<  dimensions of the performed splits */
        std::vector<Leaf> leaves;       /**<  leaves of tree containing intervals and approximating value */
        // idea: save intervals as interval-tree with nodes and corresponding values
};

std::set<int> DecisionTree::get_split_dims() const{
    return split_dims;
}

std::vector<Leaf> DecisionTree::get_leaves() const{
    return leaves;
}

typedef std::map<std::set<int>, std::shared_ptr<DecisionTree>, setComp> TreeFamily;

const double INF = std::numeric_limits<double>::infinity();

namespace rpf{
  
  /**
   * \brief A split performed with a score at a leaf_index in tree_index.
   * 
   * Remembers data for the two news leaves.
   */
  struct Split {
      double min_sum;                 /**< minimal achievable sum of squared residuals */
      std::shared_ptr<DecisionTree> tree_index;   /**< pointer to tree */
      Leaf* leaf_index;               /**< pointer to leaf containing interval */
      int split_coordinate;           /**< coordinate for splitting */
      double split_point;             /**< splitpoint */
      std::vector<int> I_s;           /**< individuals smaller than splitpoint */
      std::vector<int> I_b;           /**< individuals bigger than splitpoint */
      double M_s;                     /**< mean or median of individuals smaller than splitpoin */
      double M_b;                     /**< mean or median of individuals bigger than splitpoint */
      std::vector<double> W_s;
      std::vector<double> W_b;
      std::vector<double> Y_s;          
      std::vector<double> Y_b; 
      Split(): min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), M_s(0.0), M_b(0.0) {};
  };

  template<typename T>
  struct Matrix{
    Matrix(){}
    Matrix(std::vector<int> dimensions, T initial_value = 0) 
      : dims(dimensions) {
      for(auto d: dimensions) n_entries *= d;
      entries = std::vector<T>(n_entries, initial_value);
      assert(n_entries == entries.size());
    }
    T& operator[](std::vector<int> indices){
      assert(indices.size() == dims.size());
      int index = 0;
      for(int i=0; i<indices.size(); ++i){
        int a = indices[i];
        for(int d=0; d<dims.size(); ++d){
          if(d>i){
            a *= dims[d];
          } 
        }
        index += a;
      }
      assert(index<n_entries);
      return entries[index];
    }
    std::vector<int> dims;
    int n_entries = 1;
    
  private:
    std::vector<T> entries;
  };

}

namespace NDGrid{
  
  typedef std::vector<int> Dimension;
  typedef std::vector< Dimension > Space;
  typedef std::vector< typename Dimension::iterator > Point;
  
  class NDGrid{
    
  private:
    Space space;
    bool first = true;
    Point current;
    
  public:
    
    NDGrid(){}
    
    NDGrid(std::vector<int> dims) : dimensions(dims){
      // fill space with dimensions
      for(const auto& dim: dims){
        Dimension d;
        for(int i=1; i<=dim; ++i){
          d.push_back(i);
        }
        space.push_back(d);
      }
    }
    
    Dimension dimensions;
    std::vector<int> getPoint(){
      std::vector<int> gridPoint(current.size());
      for(int i=0; i<current.size(); ++i){
        gridPoint[i] = *current[i];
      }
      return gridPoint;
    };
    
    // the loop over space and the function-pointer to call at each point
    bool nextPoint(){
      
      // get first point in N-dimensional space
      if(first){
        first = false;
        for( Space::iterator dims_it = space.begin() ; dims_it!=space.end() ; ++dims_it ){
          current.push_back( (*dims_it).begin() );
        }
        return false;
      }
      
      // go to next point in space
      Space::iterator dims_it = space.begin();
      Point::iterator cur_it = current.begin();
      for(  ; dims_it!=space.end() ; ++dims_it, ++cur_it ){
        
        // check if next in dimension is at the end
        if ( ++(*cur_it) == (*dims_it).end() ){
          
          // check if we have finished whole space
          if( dims_it == space.end() - 1 ){
            // reset setup for next time of iteration
            first = true;
            current.clear();
            // stop running now
            return true;
          }
          // reset this dimension to begin
          // and go to next dimension
          *cur_it = (*dims_it).begin();
        }else{
          // next point is valid
          break;
        }
      }
      return false;
    }
  };

}


// ----------------- helper functions -----------------

/**
 * \brief Check whether a tree with specified split_dims already exists in tree_family
 * 
 * \param split_dims defining the tree to be searched for.
 * \param tree_family the family to be tested whether containing the tree.
 */
std::shared_ptr<DecisionTree> treeExists(const std::set<int> split_dims, TreeFamily &tree_family){
    if(tree_family.find(split_dims) != tree_family.end()) return tree_family[split_dims];
    return nullptr;
}

/**
 * \brief Check whether a tree with resulting_dims for a split_coordinate is already in possible_splits
 * 
 * \param dim defining the dimension of the split.
 * \param possible_splits containing all possible splits.
 * \param resulting_dims as union set of split dimension and dimensions of tree which is splitted.
 */
bool possibleExists(const int dim, const std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, const std::set<int> &resulting_dims){
    for(auto& elem:possible_splits){
        if(elem.first == dim && elem.second->get_split_dims() == resulting_dims) return 1;
    }
    return 0;
}

/**
 * \brief Check whether a tree has a leaf with specific interval.
 * 
 * \param interval to be compared with.
 * \param tree to be searched for leaf with interval.
 */
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

/**
 * \brief Extract keys from a std::map as vector of arbitrary type.
 * 
 * \param map with arbitrary key and value type.
 */
template <typename KT, typename VT>
std::vector<KT> getKeys(std::map<KT, VT, setComp> m){
  std::vector<KT> keys;
  for(const auto& entry: m){
    keys.push_back(entry.first);
  }
  return keys;
}

void testSetComp(){
  TreeFamily m;
  m.insert({{1,2,3},nullptr});
  m.insert({{1,3},nullptr});
  m.insert({{1,3},nullptr});
  m.insert({{1,2,3,4},nullptr});
  m.insert({{2,3,4},nullptr});
  m.insert({{1,2},nullptr});
  getKeys(m);
}

/**
 * \brief Calculates median of the vector.
 * 
 * \param vec a vector of arbitrary type.
 */
template <typename VT>
VT calcMedian(std::vector<VT> vec){
  // sort vector
  std::sort(vec.begin(), vec.end());
  size_t s = vec.size();
  
  // differ between even and odd case
  if(s % 2 == 0)
    return (vec[(s - 1) / 2] + vec[s / 2]) / 2;
  return vec[s / 2];
}

/**
 * \brief Calculate mean of a vector.
 * 
 * \param vec a vector of arbitrary type.
 */
template <typename VT>
VT calcMean(std::vector<VT> vec){
  if(vec.empty()) return 0;
  return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

template <typename CT>
struct CreateTreeFamilies : public Worker {
  std::vector<Leaf> initial_leaves;
  CT* forest;
  CreateTreeFamilies(const std::vector<Leaf> initial_leaves, CT* forest) 
    : initial_leaves(initial_leaves), forest(forest) {}
  void operator()(std::size_t begin, std::size_t end){
    for(size_t n=begin; n<end; ++n){
      forest->CT::create_tree_family(initial_leaves, n);
    }
  };
};


// ----------------- main rpf class -----------------

/**
 * \brief Create a prediction model based on Random Forests for regression data sets.
 */
class RandomPlantedForest {
  
    friend struct CreateTreeFamilies<RandomPlantedForest>;
  
    public:
        RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                            const NumericVector parameters={1,50,30,10,0.4,0,0,0,0});
        RandomPlantedForest() {};
        void set_data(const NumericVector &samples_Y, const NumericMatrix &samples_X);
        NumericVector predict_matrix(const NumericMatrix &X, const NumericVector components = {0});
        NumericVector predict_vector(const NumericVector &X, const NumericVector components = {0});
        void purify();
        void new_purify();
        void print();
        void cross_validation(int n_sets=4, IntegerVector splits={5,50}, NumericVector t_tries={0.2,0.5,0.7,0.9}, IntegerVector split_tries={1,2,5,10});
        double MSE(const NumericVector &Y_predicted, const NumericVector &Y_true); 
        void get_parameters();
        void set_parameters(StringVector keys, NumericVector values);
        List get_model();
        
    protected:
        std::vector<std::vector<double>> X;         /**< Nested vector feature samples of size (sample_size x feature_size) */
        std::vector<double> Y;                      /**< Corresponding values for the feature samples */          
        int max_interaction;                        /**< Maximum level of interaction determining maximum number of split dimensions for a tree */
        int n_trees;                                /**< Number of trees generated per family */
        int n_splits;                               /**< Number of performed splits for each tree family */
        std::vector<int> n_leaves;                  /**< */
        double t_try;                               /**< */
        int split_try;                              /**< */
        int feature_size;                           /**< Number of feature dimension in X */
        int sample_size;                            /**< Number of samples of X */
        bool purify_forest;                         /**< Whether the forest should be purified */
        bool purified = false;                      /**< Track if forest is currently purified */
        bool deterministic = false;                 /**< Choose whether approach deterministic or random */
        bool parallelize = true;                    /**< Perform algorithm in parallel or serialized */
        bool cross_validate = false;                /**< Determines if cross validation is performed */
        std::vector<double> upper_bounds;           
        std::vector<double> lower_bounds;           
        std::vector<TreeFamily> tree_families;      /**<  random planted forest containing result */
        double predict_single(const std::vector<double> &X, std::set<int> component_index);
        void L2_loss(rpf::Split &split);
        virtual void fit();
        virtual void create_tree_family(std::vector<Leaf> initial_leaves, size_t n);
        virtual rpf::Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                              std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family);
};

void RandomPlantedForest::L2_loss(rpf::Split &split){
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);
  split.M_b = calcMean(split.Y_b); 
  std::for_each(split.Y_b.begin(), split.Y_b.end(), [&split](double val){ split.min_sum += pow(val - split.M_b, 2) - pow(val, 2); });
  std::for_each(split.Y_s.begin(), split.Y_s.end(), [&split](double val){ split.min_sum += pow(val - split.M_s, 2) - pow(val, 2); });
}

// constructor
RandomPlantedForest::RandomPlantedForest(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                                         const NumericVector parameters){

    // initialize class members 
    std::vector<double> pars = to_std_vec(parameters);
    if(pars.size() != 9){
        std::cout << "Wrong number of parameters - set to default." << std::endl;
        this->max_interaction = 1;
        this->n_trees = 50;
        this->n_splits = 30;
        this->split_try = 10;
        this->t_try = 0.4;
        this->purify_forest = 0;
        this->deterministic = 0;
        this->parallelize = 0;
        this->cross_validate = 0;
    }else{
        this->max_interaction = pars[0];
        this->n_trees = pars[1];
        this->n_splits = pars[2];
        this->split_try = pars[3];
        this->t_try = pars[4];
        this->purify_forest = pars[5];
        this->deterministic = pars[6];
        this->parallelize = pars[7];
        this->cross_validate = pars[8];
    }

    // set data and data related members
    this->set_data(samples_Y, samples_X);
}

// determine optimal split
rpf::Split RandomPlantedForest::calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                               std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family){

  rpf::Split curr_split, min_split;
  std::set<int> tree_dims;
  int k;
  bool splitable;
  size_t n = 0;
  double leaf_size, sample_point;
  
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
  
  // consider a fraction of possible splits
  while(n < n_candidates){
    
    // in the beginning not known if split viable
    splitable = false;
    
    // since size of possible splits changes, check if candidate in range
    if(possible_splits.empty()) break;
    if(split_candidates[n] >= possible_splits.size()) continue;
    
    auto candidate = possible_splits.begin();
    std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
    k = candidate->first - 1; // split dim of current candidate, converted to index starting at 0
    leaf_size = n_leaves[k];
    
    // Test if splitting in the current tree w.r.t. the coordinate "k" is an element of candidate tree
    tree_dims = candidate->second->split_dims;
    tree_dims.erase(k+1);
    tree_dims.erase(0);
    
    // consider only null tree or tree with same dims as candidate or with same dims excluding the splitting coordinate
    std::vector<std::shared_ptr<DecisionTree>> curr_trees{curr_family[std::set<int>{0}]};
    if(curr_family.find(candidate->second->split_dims) != curr_family.end()) curr_trees.push_back(curr_family[candidate->second->split_dims]);
    if(curr_family.find(tree_dims) != curr_family.end()) curr_trees.push_back(curr_family[tree_dims]);
    
    // go through all trees in current family
    for(auto& curr_tree: curr_trees){
      
      // skip if tree has no leaves
      if(curr_tree->leaves.size() == 0) continue; 
      
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
        splitable = true;
        
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
          
          // clear current split
          {
            curr_split.I_s.clear();
            curr_split.I_b.clear();
            curr_split.Y_s.clear();
            curr_split.Y_b.clear();
          }
          
          // get samples greater/smaller than samplepoint
          for(int individual: curr_individuals){
            if(X[individual][k] < sample_point){
              curr_split.Y_s.push_back(Y[individual]);
              curr_split.I_s.push_back(individual);
            }else{
              curr_split.Y_b.push_back(Y[individual]);
              curr_split.I_b.push_back(individual);
            }
          }
          
          // ensure individuals are sorted and unique
          std::sort(curr_split.I_s.begin(), curr_split.I_s.end());
          std::sort(curr_split.I_b.begin(), curr_split.I_b.end());
          curr_split.I_s.erase(std::unique(curr_split.I_s.begin(), curr_split.I_s.end()), curr_split.I_s.end());
          curr_split.I_b.erase(std::unique(curr_split.I_b.begin(), curr_split.I_b.end()), curr_split.I_b.end());
          
          // accumulate squared mean and get mean
          L2_loss(curr_split);

          // update split if squared sum is smaller
          if(curr_split.min_sum < min_split.min_sum){
            min_split.min_sum = curr_split.min_sum;
            min_split.tree_index = curr_tree;
            min_split.leaf_index =  &leaf;
            min_split.split_coordinate = k+1;
            min_split.split_point = sample_point;
            min_split.I_s = curr_split.I_s;
            min_split.I_b = curr_split.I_b;
            min_split.M_s = curr_split.M_s;
            min_split.M_b = curr_split.M_b;
            min_split.Y_s = curr_split.Y_s;
            min_split.Y_b = curr_split.Y_b;
          }
        }
      }
    }
    
    // if split viable, increase count, otherwise remove candidate
    if(splitable){
      ++n;
    }else{
      possible_splits.erase(candidate);
    }
  }
  
  return min_split;
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
    
    this->fit();
    
    if(cross_validate){
      this->cross_validation();
    }
}

void RandomPlantedForest::create_tree_family(std::vector<Leaf> initial_leaves, size_t n){
  
  TreeFamily curr_family;
  curr_family.insert(std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))); // save tree with one leaf in the beginning
  
  // store possible splits in map with splitting variable as key and pointer to resulting tree
  std::multimap<int, std::shared_ptr<DecisionTree>> possible_splits;
  for(int feature_dim = 1; feature_dim<=feature_size; ++feature_dim){
      // add pointer to resulting tree with split dimension as key
      curr_family.insert(std::make_pair(std::set<int>{feature_dim}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{feature_dim}))));
      possible_splits.insert(std::make_pair(feature_dim, curr_family[std::set<int>{feature_dim}]));
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
          curr_split.M_s << "/" << curr_split.M_b <<  std::endl;
      }

      // update possible splits
      for(int feature_dim = 1; feature_dim<=feature_size; ++feature_dim){ // consider all possible dimensions
        
        // ignore dim if same as split coordinate or in dimensions of old tree
        if(feature_dim == curr_split.split_coordinate || curr_split.tree_index->split_dims.count(feature_dim) > 0) continue;
        
        // create union of split coord, feature dim and dimensions of old tree
        std::set<int> curr_dims = curr_split.tree_index->split_dims;
        curr_dims.insert(curr_split.split_coordinate);
        curr_dims.insert(feature_dim);
        curr_dims.erase(0);
        
        // do not exceed maximum level of interaction
        if(curr_dims.size() > max_interaction) continue; 
        
        // skip if possible_split already exists
        if(possibleExists(feature_dim, possible_splits, curr_dims)) continue;
        
        // check if resulting tree already exists in family
        std::shared_ptr<DecisionTree> found_tree = treeExists(curr_dims, curr_family);
        
        // update possible_splits if not already existing
        if(found_tree){ // if yes add pointer
          possible_splits.insert(std::make_pair(feature_dim, found_tree));
        }else{ // if not create new tree
          curr_family.insert(std::make_pair(curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))));
          possible_splits.insert(std::make_pair(feature_dim, curr_family[curr_dims]));
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
      
      // update values of individuals of split interval with mean
      for(int individual: curr_split.leaf_index->individuals){ // todo: loop directly over I_s I_b
        if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
          samples_Y[individual] -= curr_split.M_s;
        }else{
          samples_Y[individual] -= curr_split.M_b;
        }
      }
      
      // construct new leaves
      Leaf leaf_s, leaf_b;
      {
        leaf_s.individuals = curr_split.I_s;
        leaf_b.individuals = curr_split.I_b;
        leaf_s.value = curr_split.M_s;
        leaf_b.value = curr_split.M_b;
        
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
      
      // construct split_dims of resulting tree when splitting in split_coordinate
      std::set<int> resulting_dims = curr_split.tree_index->split_dims;
      resulting_dims.insert(curr_split.split_coordinate);
      resulting_dims.erase(0);
      
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
        std::cout << "Current TreeFamily: (" << curr_family.size() << ") ";
        for(auto tree: curr_family){
          std::cout << "Dims = ";
          for(auto dim: tree.first) std::cout << dim << ", ";
          std::cout << "; " << "Number of Leaves = " << tree.second->leaves.size();
          std::cout << " / ";
        }
        std::cout << std::endl << std::endl;
      }
    }
  }
  
  // remove empty trees & clear individuals of each tree
  auto keys = getKeys(curr_family);
  for(auto& key: keys){
    if(curr_family[key]->leaves.size() == 0){
      curr_family.erase(key);
      continue;
    } 
    for(auto leaf: curr_family[key]->leaves){
      leaf.individuals.clear();
    }
  }
  
  tree_families[n] = curr_family;
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
      int n_threads = std::thread::hardware_concurrency()-1;
      for(int n = 0; n<n_trees; n+=n_threads){
        if(n>=(n_trees-n_threads)) n_threads = n_trees % n_threads;
        std::vector<std::thread> threads(n_threads);
        for(size_t t=0; t<n_threads; ++t){
          threads[t] = std::thread(&RandomPlantedForest::create_tree_family, this, std::ref(initial_leaves), n+t);
        }
        for(auto& t: threads){
          if(t.joinable()) t.join();
        }
      }
    }else{
      for(size_t n=0; n<n_trees; ++n){
        create_tree_family(initial_leaves, n);
      }
    }

    // optionally purify tree
    if(purify_forest){
      this->purify();
    }else{
      purified = false;
    }
}

void RandomPlantedForest::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries){
    
    bool cv_tmp = this->cross_validate;
    this->cross_validate = false;
  
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
    this->cross_validate = cv_tmp;
    
    std::cout << "Optimal parameters: " << optimal_inter << ", " << optimal_split << ", " << optimal_t_try << ", " << optimal_split_try << ": MSE=" << MSE_min << std::endl;
}

// predict single feature vector
double RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index){

    double total_res = 0;
    
    // consider all components
    if(component_index == std::set<int>{0}) {
        for(auto& tree_family: this->tree_families){
            for(auto& tree: tree_family){
                for(auto& leaf: tree.second->leaves){
                    bool valid = true;
                    for(auto& dim: tree.first){
                        if(!(leaf.intervals[std::max(0, dim-1)].first <= X[std::max(0, dim-1)]
                            && (leaf.intervals[std::max(0, dim-1)].second > X[std::max(0, dim-1)]
                            || leaf.intervals[std::max(0, dim-1)].second == upper_bounds[std::max(0, dim-1)]))){
                            valid = false;
                        }
                    }
                    if(valid) total_res += leaf.value;
                }
            }
            // todo: check if constant not needed with null tree
            // total_res += tree_family.constant;
        }
    }else{ // choose components for prediction
        for(auto& tree_family: this->tree_families){
            for(auto& tree: tree_family){
                
                // only consider trees with same dimensions as component_index
                if(tree.first != component_index) continue;
                
                std::vector<int> dims;
                for(auto dim: tree.first) {dims.push_back(dim);}
                
                for(auto& leaf: tree.second->leaves){
                    bool valid = true;
                    for(int i = 0; i<dims.size(); ++i){

                        int dim = dims[i];
                      
                        if(!(leaf.intervals[std::max(0, dim-1)].first <= X[i]
                            && (leaf.intervals[std::max(0, dim-1)].second > X[i]
                            || leaf.intervals[std::max(0, dim-1)].second == upper_bounds[std::max(0, dim-1)]))){
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
    for(auto tree: curr_family){
      if(tree.first.size() > curr_max) curr_max = tree.first.size();
    }
    
    while(curr_max >= 1){
      
      // go through split dimensions of all trees
      auto keys = getKeys(curr_family);
      std::vector<std::set<int>>::reverse_iterator key = keys.rbegin();
      while(key != keys.rend()){
        
        auto& curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;
        
        // check if number of dims same as current max_interaction
        if(curr_dims.size() == curr_max){
            
          // go through feature dims
          for(int feature_dim=1; feature_dim<=feature_size; ++feature_dim){
            
            // continue only if dim in current tree
            if(curr_tree->split_dims.count(feature_dim) != 0){
              
              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree
              
              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if(curr_max == 1){
                tree = curr_family[std::set<int>{0}];
              }else{
                if(!tree){
                  curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                  tree = curr_family[tree_dims];
                }
              }
              
              // go through leaves of current tree
              int n_leaves = curr_tree->leaves.size();
              for(int l=0; l<n_leaves; ++l){
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
        key++;
      }
      
      // update currently considered dimension size
      --curr_max;
    }
  }
}

void RandomPlantedForest::new_purify(){

  // go through all n_trees families 
  //for(const auto& curr_family: this->tree_families){
  auto& curr_family = this->tree_families[0];
  
    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);
    
    // go through all variables of the component
    for(int curr_dim=1; curr_dim<=feature_size; ++curr_dim){
      std::vector<double> bounds;
      
      // go through trees of family
      for(const auto& curr_tree: curr_family){

        // consider only relevant trees that have current dimension as variable
        if(!curr_tree.first.count(curr_dim)) continue;
        
        // go through leaves of tree
        for(const auto& curr_leaf: curr_tree.second->leaves){
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim-1].second);
        }
      }
      std::sort( bounds.begin(), bounds.end());
      bounds.erase( std::unique(bounds.begin(), bounds.end()), bounds.end());
      lim_list[curr_dim-1] = bounds;
    }
    
    std::vector< NDGrid::NDGrid> grids(curr_family.size() - 1); // ignore null tree?

    // initialize values and individuals for each tree in family
    std::vector<rpf::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<rpf::Matrix<double>> values(curr_family.size() - 1);
    
    // setup finer grid with individuals and values
    int tree_index = 0; // ignore null tree?
    for(const auto& curr_tree: curr_family)
    {
      std::cout << tree_index << ": " << std::endl;
      
      if(curr_tree.first == std::set<int>{0}) continue; // ignore null tree?
      
      // fill space with dimensions
      std::vector<int> dimensions;
      
      // get variables and individuals of tree
      for(const auto& var: curr_tree.first){ // go through dimensions of tree
        dimensions.push_back(lim_list[var - 1].size() - 1); // size - 1 ?
      }
      
      // setup grid for leaf indices
      auto grid =  NDGrid::NDGrid(dimensions);
      
      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = rpf::Matrix<int>(dimensions);
      values[tree_index] = rpf::Matrix<double>(dimensions);
      
      // loop over grid points
      while(!grid.nextPoint()){
        std::vector<int> gridPoint = grid.getPoint();
      
        // go through leaves of tree
        for(const auto& leaf: curr_tree.second->get_leaves()){
          bool in_leaf_val = false, in_leaf_ind = false;
          int dim_index = 0;
          for(const auto& dim: curr_tree.first){
            
            // consider values only if any in
            if( (leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) 
                 && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1]) ) in_leaf_val = true;
            
            // consider individuals only if all in
            bool all_true = true;
            for(const auto& x: X[dim]){
              if( !( (x >= lim_list[dim - 1][gridPoint[dim_index]]) 
                  && (x >= lim_list[dim - 1][gridPoint[dim_index] + 1]) ) ) all_true = false;
            }
            if(all_true) in_leaf_ind = true;
            
            ++dim_index;
          }
          
          // sum up individuals and values
          //if(in_leaf_val) values[tree_index][gridPoint] += double(leaf.value); 
          if(in_leaf_ind) individuals[tree_index][gridPoint] += 1; 
        }
      }
      
      
      while(!grid.nextPoint()){
        std::cout << individuals[tree_index][grid.getPoint()] << ", ";
      }
      std::cout << std::endl;
      
      ++tree_index;
    }
    
    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for(auto tree: curr_family){
      if(tree.first.size() > curr_max) curr_max = tree.first.size();
    }
    
    // create new trees
    std::vector<std::set<int>> new_trees;
    while(curr_max >= 1){
      
      // go through split dimensions of all trees
      auto keys = getKeys(curr_family);
      int tree_index = 0;
      std::vector<std::set<int>>::reverse_iterator key = keys.rbegin();
      while(key != keys.rend()){
        
        auto& curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;
        
        if(false){
          std::cout << "Old tree: ";
          for(const auto& dim: curr_dims) std::cout << dim << ", "; 
          std::cout << std::endl;
        }
        
        // check if number of dims same as current max_interaction
        if(curr_dims.size() == curr_max){
          
          // go through feature dims
          int dim_index = 0;
          for(int feature_dim=1; feature_dim<=feature_size; ++feature_dim){
            
            // continue only if dim in current tree
            if(curr_tree->split_dims.count(feature_dim) != 0){
              
              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree
              
              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if(curr_max == 1){
                tree = curr_family[std::set<int>{0}];
              }else{
                if(!tree){
                  if(false){
                    std::cout << "   New tree: ";
                  for(const auto& dim: tree_dims) std::cout << dim << ", "; 
                  std::cout << std::endl;
                  }
                  curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                  new_trees.push_back(tree_dims);
                  std::vector<int> matrix_dimensions = values[tree_index].dims;
                  matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);
                  values.push_back(rpf::Matrix<double>(matrix_dimensions));
                  // todo: add individuals to new trees
                  tree = curr_family[tree_dims];
                  std::cout << std::endl;
                }
              }
              dim_index++;
            }
          }
        }
        key++;
        tree_index++;
      }
      
      // update currently considered dimension size
      --curr_max;
    }
    
    // todo: purify
    
  //}
}

void RandomPlantedForest::print(){
    for(int n=0; n<n_trees; ++n){
        TreeFamily family = tree_families[n];
        // todo: check if constant needed
        // std::cout << n+1 << " TreeFamily: constant=" << family.constant << std::endl << std::endl;
        auto keys = getKeys(family);
        for(int m=0; m<keys.size(); ++m){
            DecisionTree tree = *(family[keys[m]]);
            std::cout << m+1 << " Tree: ";
            std::cout << "Dims=";
            for(const auto& dim: tree.split_dims) std::cout << dim << ",";
            std::cout << std::endl << "Leaves: (" << tree.leaves.size() << ")" << std::endl;
            for(const auto& leaf: tree.leaves){
                std::cout << "Intervals=";
                for(const auto& interval: leaf.intervals){
                    std::cout << interval.first << "," << interval.second << "/";
                }
                std::cout << " Value=" << leaf.value << std::endl;
            }  
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}

// print parameters of the model to the console
void RandomPlantedForest::get_parameters(){
  std::cout << "Parameters: n_trees=" <<  n_trees << ", n_splits=" << n_splits << ", max_interaction=" << max_interaction << ", t_try=" << t_try 
            << ", split_try=" << split_try << ", purified=" << purified << ", deterministic=" << deterministic << ", parallel=" << parallelize
            << ", feature_size=" << feature_size << ", sample_size=" << sample_size << std::endl;
}

/*  retrospectively change parameters of existing class object, 
    updates the model, so far only single valued parameters supported,
    for replacing training data use 'set_data', 
    note that changing cv does not trigger cross validation */
void RandomPlantedForest::set_parameters(StringVector keys, NumericVector values){
  if(keys.size() != values.size()) {
    std::cout << "Size of input vectors is not the same. " << std::endl;
    return;
  }
  
  for(int i=0; i<keys.size(); ++i){
    if(keys[i] == "deterministic"){
      this->deterministic = values[i];
    }else if(keys[i] == "parallel"){
      this->parallelize = values[i];
    }else if(keys[i] == "purify"){
      this->purify_forest = values[i];
    }else if(keys[i] == "n_trees"){
      this->n_trees = values[i];
    }else if(keys[i] == "n_splits"){
      this->n_splits= values[i];
    }else if(keys[i] == "t_try"){
      this->t_try = values[i];
    }else if(keys[i] == "split_try"){
      this->split_try = values[i];
    }else if(keys[i] == "max_interaction"){
      this->max_interaction = values[i];
    }else if(keys[i] == "cv"){
      this->cross_validate = values[i];
    }else{
      std::cout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}

List RandomPlantedForest::get_model(){
  List model;
  for(const auto family: tree_families){
    List variables, family_values, family_intervals;
    for(const auto tree: family){
      NumericVector tree_values;
      List tree_intervals;
      variables.push_back(from_std_set(tree.first));
      for(const auto leaf: tree.second->leaves){
        tree_values.push_back(leaf.value);
        NumericVector intervals;
        for(const auto interval: leaf.intervals){
          intervals.push_back(interval.first);
          intervals.push_back(interval.second);
        }
        NumericMatrix leaf_intervals(2, feature_size, intervals.begin());
        tree_intervals.push_back(leaf_intervals);
      }
      family_intervals.push_back(tree_intervals);
      family_values.push_back(tree_values);
    }
    model.push_back(List::create(Named("variables") = variables, _["values"] = family_values, _["intervals"] = family_intervals));
  }
  return(model);
}


// ----------------- rpf subclass for classification -----------------

/**
 * \brief Create a prediction model based on Random Forests for classification data sets.
 */
class ClassificationRPF : public RandomPlantedForest {
  
  friend struct CreateTreeFamilies<ClassificationRPF>;
  
  public:
    ClassificationRPF(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                      const String loss="L2", const NumericVector parameters={1,50,30,10,0.4,0,0,0,0,0,0.1});
    void set_parameters(StringVector keys, NumericVector values);
    
  private:
    double delta;
    double epsilon;
    enum LossType { L1, L2, median, logit, exponential};
    LossType loss; 
    void (ClassificationRPF::*calcLoss)(rpf::Split&);                                               // function pointer
    std::vector<double> weights; 
    void create_tree_family(std::vector<Leaf> initial_leaves, size_t n) override;
    rpf::Split calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                               std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family) override;
    void fit() override;
    void L1_loss(rpf::Split &split);
    void median_loss(rpf::Split &split);
    void logit_loss(rpf::Split &split);
    void exponential_loss(rpf::Split &split);
};

void ClassificationRPF::L1_loss(rpf::Split &split){
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);
  split.M_b = calcMean(split.Y_b); 
  std::for_each(split.Y_b.begin(), split.Y_b.end(), [&split](double val){ split.min_sum += abs(val - split.M_b) - abs(val); });
  std::for_each(split.Y_s.begin(), split.Y_s.end(), [&split](double val){ split.min_sum  += abs(val - split.M_s) - abs(val); });
}

void ClassificationRPF::median_loss(rpf::Split &split){
  split.min_sum = 0;
  split.M_s = calcMedian(split.Y_s);
  split.M_b = calcMedian(split.Y_b); 
  std::for_each(split.Y_b.begin(), split.Y_b.end(), [&split](double val){ split.min_sum += abs(val - split.M_b) - abs(val); });
  std::for_each(split.Y_s.begin(), split.Y_s.end(), [&split](double val){ split.min_sum += abs(val - split.M_s) - abs(val); });
}

void ClassificationRPF::logit_loss(rpf::Split &split){
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);  
  split.M_b = calcMean(split.Y_b); 
  double M_s = std::min(1 - delta, std::max(delta, split.M_s)); // ~ R3
  double M_b = std::min(1 - delta, std::max(delta, split.M_b)); // ~ R2
  double W_s_mean = calcMean(split.W_s), W_b_mean = calcMean(split.W_b);
  
  std::vector<double> W_s = split.W_s, W_b = split.W_b;
  for(auto &w: W_s) w = 1 / (1 + exp( -(w + log(M_s / (1 - M_s)) - W_s_mean))); // W_3 -> P_3
  for(auto &w: W_b) w = 1 / (1 + exp( -(w + log(M_b / (1 - M_b)) - W_b_mean))); // W_2 -> P_2
  
  for(int i=0; i<split.W_s.size(); ++i) split.min_sum += split.Y_s[i] * log(1/ (1 + exp(-split.W_s[i]))) + (1 - split.Y_s[i]) * log(1 - (1/ (1 + exp(-split.W_s[i])))) ; // ~ R_old
  for(int i=0; i<split.W_b.size(); ++i) split.min_sum += split.Y_b[i] * log(1/ (1 + exp(-split.W_b[i]))) + (1 - split.Y_b[i]) * log(1 - (1/ (1 + exp(-split.W_b[i])))) ; // ~ R_old
  
  for(int i=0; i<split.Y_s.size(); ++i) split.min_sum -= split.Y_s[i] * log(W_s[i]) + (1 - split.Y_s[i]) * log(1 - W_s[i]); // Y_3
  for(int i=0; i<split.Y_b.size(); ++i) split.min_sum -= split.Y_b[i] * log(W_b[i]) + (1 - split.Y_b[i]) * log(1 - W_b[i]); // Y_2
  
  if(std::isnan(split.min_sum)){
    split.min_sum = INF;
  }
}

void ClassificationRPF::exponential_loss(rpf::Split &split){
  split.min_sum = 0;
  
  double W_s_sum = std::accumulate(split.W_s.begin(), split.W_s.end(), 0.0);
  double W_b_sum = std::accumulate(split.W_b.begin(), split.W_b.end(), 0.0);
  
  double sum_s = 0, sum_b = 0;
  for(int i=0; i<split.Y_s.size(); ++i) sum_s += ((split.Y_s[i] + 1) / 2) * (split.W_s[i] / W_s_sum);
  for(int i=0; i<split.Y_b.size(); ++i) sum_b += ((split.Y_b[i] + 1) / 2) * (split.W_b[i] / W_b_sum);
  split.M_s = sum_s;
  split.M_b = sum_b;
  sum_s = std::min(1 - delta, std::max(delta, sum_s));
  sum_b = std::min(1 - delta, std::max(delta, sum_b));
  
  for(int i=0; i<split.Y_s.size(); ++i) split.min_sum += split.W_s[i] * exp(-0.5 * split.Y_s[i] * log(sum_s / (1 - sum_s)));
  for(int i=0; i<split.Y_b.size(); ++i) split.min_sum += split.W_b[i] * exp(-0.5 * split.Y_b[i] * log(sum_b / (1 - sum_b)));
  
  split.min_sum -= W_s_sum + W_b_sum;
  
  if(W_b_sum == 0 || std::isnan(split.min_sum)){
    split.min_sum = INF;
  }
}

// constructor with parameters split_try, t_try, purify_forest, deterministic, parallelize
ClassificationRPF::ClassificationRPF(const NumericVector &samples_Y, const NumericMatrix &samples_X,
                                     const String loss, const NumericVector parameters)
   : RandomPlantedForest{}{

  // initialize class members 
  std::vector<double> pars = to_std_vec(parameters);
  if(loss == "L1"){
    this->loss = LossType::L1;
    this->calcLoss = &ClassificationRPF::L1_loss;
  }else if(loss == "L2"){
    this->loss = LossType::L2; 
    this->calcLoss = &ClassificationRPF::L2_loss;
  }else if(loss == "median"){
    this->loss = LossType::median;
    this->calcLoss = &ClassificationRPF::median_loss;
  }else if(loss == "logit"){
    this->loss = LossType::logit;
    this->calcLoss = &ClassificationRPF::logit_loss;
  }else if(loss == "exponential"){
    this->loss = LossType::exponential;
    this->calcLoss = &ClassificationRPF::exponential_loss;
  }else{
    std::cout << "Unkown loss function." << std::endl;
  }
  if(pars.size() != 11){
    std::cout << "Wrong number of parameters - set to default." << std::endl;
    this->max_interaction = 1;
    this->n_trees = 50;
    this->n_splits = 30;
    this->split_try = 10;
    this->t_try = 0.4;
    this->purify_forest = 0;
    this->deterministic = 0;
    this->parallelize = 0;
    this->cross_validate = 0;
    this->delta = 0.1;
    this->epsilon = 0;
  }else{
    this->max_interaction = pars[0];
    this->n_trees = pars[1];
    this->n_splits = pars[2];
    this->split_try = pars[3];
    this->t_try = pars[4];
    this->purify_forest = pars[5];
    this->deterministic = pars[6];
    this->parallelize = pars[7];
    this->cross_validate = pars[8];
    this->delta = pars[9];
    this->epsilon = pars[10];
  }
  
  // set data and data related members
  this->set_data(samples_Y, samples_X);
}

// determine optimal split
rpf::Split ClassificationRPF::calcOptimalSplit(const std::vector<double> &Y, const std::vector<std::vector<double>> &X,
                                                 std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family){

  rpf::Split curr_split, min_split;
  std::set<int> tree_dims;
  int k;
  bool splitable;
  size_t n = 0;
  double leaf_size, sample_point;
  
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
  
  // consider a fraction of possible splits
  while(n < n_candidates){
    
    // in the beginning not known if split viable
    splitable = false;
    
    // since size of possible splits changes, check if candidate in range
    if(possible_splits.empty()) break;
    if(split_candidates[n] >= possible_splits.size()) continue;
    
    auto candidate = possible_splits.begin();
    std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
    k = candidate->first - 1; // split dim of current candidate, converted to index starting at 0
    leaf_size = n_leaves[k];
    
    // Test if splitting in the current tree w.r.t. the coordinate "k" is an element of candidate tree
    tree_dims = candidate->second->split_dims;
    tree_dims.erase(k+1);
    tree_dims.erase(0);
    
    // consider only null tree or tree with same dims as candidate or with same dims excluding the splitting coordinate
    std::vector<std::shared_ptr<DecisionTree>> curr_trees{curr_family[std::set<int>{0}]};
    if(curr_family.find(candidate->second->split_dims) != curr_family.end()) curr_trees.push_back(curr_family[candidate->second->split_dims]);
    if(curr_family.find(tree_dims) != curr_family.end()) curr_trees.push_back(curr_family[tree_dims]);
    
    // go through all trees in current family
    for(auto& curr_tree: curr_trees){
      
      // skip if tree has no leaves
      if(curr_tree->leaves.size() == 0) continue; 
      
      // go through all leaves of current tree
      for(auto& leaf: curr_tree->leaves){
        std::vector<int> curr_individuals = leaf.individuals; // consider individuals of current leaf
        
        // extract sample points according to individuals from X and Y
        std::map<double, double> unique_samples;
        for(int individual: curr_individuals){
          unique_samples[X[individual][k]] = Y[individual];
        }
        
        // check if number of sample points is within limit
        if(unique_samples.size() < 2 * leaf_size) continue;
        splitable = true;
        
        int start = 0, end = split_try;
        if(deterministic){
          start = 1;
          end = unique_samples.size() - 1;
        }
        
        // consider split_try-number of random samples
        for(int t=start; t<end; ++t){
          
          // get samplepoint
          auto sample_pos = unique_samples.begin();
          std::uniform_int_distribution<> distrib(leaf_size, unique_samples.size() - leaf_size + 1);
          if(deterministic){
            std::advance(sample_pos, t);
          }else{
            std::advance(sample_pos, distrib(gen)); // consider only sample points with offset
          }
          sample_point = sample_pos->first;
          
          // clear current split
          {
            curr_split.I_s.clear();
            curr_split.I_b.clear();
            curr_split.Y_s.clear();
            curr_split.Y_b.clear();
            curr_split.W_s.clear();
            curr_split.W_b.clear();
          }

          // get samples greater/smaller than samplepoint
          for(int individual: curr_individuals){
            if(X[individual][k] < sample_point){
              curr_split.Y_s.push_back(Y[individual]);
              curr_split.I_s.push_back(individual);
              curr_split.W_s.push_back(weights[individual]);
            }else{
              curr_split.Y_b.push_back(Y[individual]);
              curr_split.I_b.push_back(individual);
              curr_split.W_b.push_back(weights[individual]);
            }
          }
          
          // ensure individuals are sorted and unique
          std::sort(curr_split.I_s.begin(), curr_split.I_s.end());
          std::sort(curr_split.I_b.begin(), curr_split.I_b.end());
          curr_split.I_s.erase(std::unique(curr_split.I_s.begin(), curr_split.I_s.end()), curr_split.I_s.end());
          curr_split.I_b.erase(std::unique(curr_split.I_b.begin(), curr_split.I_b.end()), curr_split.I_b.end());
          
          // accumulate squared mean and get mean
          (this->*ClassificationRPF::calcLoss)(curr_split);

          // update split if squared sum is smaller
          if(curr_split.min_sum < min_split.min_sum){
            min_split.min_sum = curr_split.min_sum;
            min_split.tree_index = curr_tree;
            min_split.leaf_index =  &leaf;
            min_split.split_coordinate = k+1;
            min_split.split_point = sample_point;
            min_split.I_s = curr_split.I_s;
            min_split.I_b = curr_split.I_b;
            min_split.M_s = curr_split.M_s;
            min_split.M_b = curr_split.M_b;
            min_split.W_s = curr_split.W_s;
            min_split.W_b = curr_split.W_b;
            min_split.Y_s = curr_split.Y_s;
            min_split.Y_b = curr_split.Y_b;
          }
        }
      }
    }
    
    // if split viable, increase count, otherwise remove candidate
    if(splitable){
      ++n;
    }else{
      possible_splits.erase(candidate);
    }
  }
  
  return min_split;
}

void ClassificationRPF::create_tree_family(std::vector<Leaf> initial_leaves, size_t n){

  TreeFamily curr_family;
  curr_family.insert(std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))); // save tree with one leaf in the beginning
  
  // store possible splits in map with splitting variable as key and pointer to resulting tree
  std::multimap<int, std::shared_ptr<DecisionTree>> possible_splits;
  for(int feature_dim = 1; feature_dim<=feature_size; ++feature_dim){
    // add pointer to resulting tree with split dimension as key
    curr_family.insert(std::make_pair(std::set<int>{feature_dim}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{feature_dim}))));
    possible_splits.insert(std::make_pair(feature_dim, curr_family[std::set<int>{feature_dim}]));
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
  
  // initialize weights
  switch(this->loss){
    case LossType::logit:
      this->weights = std::vector<double>(sample_size, 0);
      break;
    case LossType::exponential:
      this->weights = std::vector<double>(sample_size, 1);
      break;
    default:
      this->weights = std::vector<double>(sample_size, 0);
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
      
      // update possible splits
      for(int feature_dim = 1; feature_dim<=feature_size; ++feature_dim){ // consider all possible dimensions
        
        // ignore dim if same as split coordinate or in dimensions of old tree
        if(feature_dim == curr_split.split_coordinate || curr_split.tree_index->split_dims.count(feature_dim) > 0) continue;
        
        // create union of split coord, feature dim and dimensions of old tree
        std::set<int> curr_dims = curr_split.tree_index->split_dims;
        curr_dims.insert(curr_split.split_coordinate);
        curr_dims.insert(feature_dim);
        curr_dims.erase(0);
        
        // do not exceed maximum level of interaction
        if(curr_dims.size() > max_interaction) continue; 
        
        // skip if possible_split already exists
        if(possibleExists(feature_dim, possible_splits, curr_dims)) continue;
        
        // check if resulting tree already exists in family
        std::shared_ptr<DecisionTree> found_tree = treeExists(curr_dims, curr_family);
        
        // update possible_splits if not already existing
        if(found_tree){ // if yes add pointer
          possible_splits.insert(std::make_pair(feature_dim, found_tree));
        }else{ // if not create new tree
          curr_family.insert(std::make_pair(curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))));
          possible_splits.insert(std::make_pair(feature_dim, curr_family[curr_dims]));
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
          
      // update values of individuals of split interval
      double update_s = 0, update_b = 0;  
      switch(this->loss){
        case LossType::L1: case LossType::L2: case LossType::median: {
          update_s = curr_split.M_s;
          update_b = curr_split.M_b;
          for(int individual: curr_split.leaf_index->individuals){
            if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              samples_Y[individual] -= update_s;
            }else{
              samples_Y[individual] -= update_b;
            }
          }
          break;
        }
        case LossType::logit: {
          double v_s = std::min(1 - epsilon, std::max(epsilon, curr_split.M_s));
          double v_b = std::min(1 - epsilon, std::max(epsilon, curr_split.M_b));
          update_s = log(v_s / (1 - v_s)) - calcMean(curr_split.W_s);
          update_b = log(v_b / (1 - v_b)) - calcMean(curr_split.W_b);
          for(int individual: curr_split.leaf_index->individuals){
            if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              weights[individual] += update_s;
            }else{
              weights[individual] += update_b;
            }
          }
          break;
        }
        case LossType::exponential: {
          double sum_s = std::min(1 - epsilon, std::max(epsilon, curr_split.M_s));
          double sum_b = std::min(1 - epsilon, std::max(epsilon, curr_split.M_b));
          if(curr_split.M_s != 0) update_s = log(sum_s / (1 - sum_s));
          if(curr_split.M_b != 0) update_b = log(sum_b / (1 - sum_b));
          for(int individual: curr_split.leaf_index->individuals){
            if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              if(std::isinf(update_s)){
                weights[individual] = 0;
              }else{
                weights[individual] *= exp(-0.5 * samples_Y[individual] * update_s);
              }
            }else{
              if(std::isinf(update_b)){
                weights[individual] = 0;
              }else{
                weights[individual] *= exp(-0.5 * samples_Y[individual] * update_b);
              }
            }
          }
          break;
        }
      }
      
      // construct new leaves
      Leaf leaf_s, leaf_b;
      {
        leaf_s.individuals = curr_split.I_s;
        leaf_b.individuals = curr_split.I_b;
        leaf_s.value = update_s;
        leaf_b.value = update_b;
        
        // initialize interval with split interval
        leaf_s.intervals = curr_split.leaf_index->intervals;
        leaf_b.intervals = curr_split.leaf_index->intervals;
        
        // interval of leaf with smaller individuals has new upper bound in splitting dimension
        leaf_s.intervals[curr_split.split_coordinate-1].second = curr_split.split_point;
        // interval of leaf with bigger individuals has new lower bound in splitting dimension
        leaf_b.intervals[curr_split.split_coordinate-1].first = curr_split.split_point;
      }

      // construct split_dims of resulting tree when splitting in split_coordinate
      std::set<int> resulting_dims = curr_split.tree_index->split_dims;
      resulting_dims.insert(curr_split.split_coordinate);
      resulting_dims.erase(0);
      
      // check if resulting tree already exists in family
      std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);
      
      // determine which tree is modified
      if(curr_split.tree_index->split_dims.count(curr_split.split_coordinate) ){ // if split variable is already in tree to be split
        // change values
        {
          leaf_s.value += curr_split.leaf_index->value + update_s;
          leaf_b.value += curr_split.leaf_index->value + update_b;
        }
        *curr_split.leaf_index = leaf_b; // replace old interval
        curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
      } else{ // otherwise 
        found_tree->leaves.push_back(leaf_s); //append new leaves
        found_tree->leaves.push_back(leaf_b);
      }
    } else{
      std::cout << "test" << std::endl;
    }
  }
  
  // remove empty trees & clear individuals of each tree
  auto keys = getKeys(curr_family);
  for(auto& key: keys){
    if(curr_family[key]->leaves.size() == 0){
      curr_family.erase(key);
      continue;
    } 
    for(auto leaf: curr_family[key]->leaves){
      leaf.individuals.clear();
    }
  }
  
  tree_families[n] = curr_family;
}

// fit forest to new data
void ClassificationRPF::fit(){
  
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
    CreateTreeFamilies<ClassificationRPF> create_tree_families(initial_leaves, this);
    parallelFor(0, n_trees, create_tree_families);
  }else{
    for(size_t n=0; n<n_trees; ++n){
      create_tree_family(initial_leaves, n);
    }
  }
  
  // optionally purify tree
  if(purify_forest){
    this->purify();
  }else{
    purified = false;
  }
}

/*  retrospectively change parameters of existing class object, 
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data', 
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(StringVector keys, NumericVector values){
  if(keys.size() != values.size()) {
    std::cout << "Size of input vectors is not the same. " << std::endl;
    return;
  }
  
  for(int i=0; i<keys.size(); ++i){
    if(keys[i] == "deterministic"){
      this->deterministic = values[i];
    }else if(keys[i] == "parallel"){
      this->parallelize = values[i];
    }else if(keys[i] == "purify"){
      this->purify_forest = values[i];
    }else if(keys[i] == "n_trees"){
      this->n_trees = values[i];
    }else if(keys[i] == "n_splits"){
      this->n_splits= values[i];
    }else if(keys[i] == "t_try"){
      this->t_try = values[i];
    }else if(keys[i] == "split_try"){
      this->split_try = values[i];
    }else if(keys[i] == "max_interaction"){
      this->max_interaction = values[i];
    }else if(keys[i] == "cv"){
      this->cross_validate = values[i];
    }else if(keys[i] == "loss"){
      if(keys[i] == "L1"){
        this->loss = LossType::L1;
        this->calcLoss = &ClassificationRPF::L1_loss;
      }else if(keys[i] == "L2"){
        this->loss = LossType::L2; 
        this->calcLoss = &ClassificationRPF::L2_loss;
      }else if(keys[i] == "median"){
        this->loss = LossType::median;
        this->calcLoss = &ClassificationRPF::median_loss;
      }else if(keys[i] == "logit"){
        this->loss = LossType::logit;
        this->calcLoss = &ClassificationRPF::logit_loss;
      }else if(keys[i] == "exponential"){
        this->loss = LossType::exponential;
        this->calcLoss = &ClassificationRPF::exponential_loss;
      }else{
        std::cout << "Unkown loss function." << std::endl;
      }
    }else if(keys[i] == "delta"){
      this->delta = values[i];
    }else if(keys[i] == "epsilon"){
      this->epsilon = values[i];
    }else{
      std::cout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}


// ----------------- Rcpp include  -----------------

RCPP_MODULE(mod_rpf) {

    class_<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericVector, const NumericMatrix, const NumericVector>()
      .method("set_data", &RandomPlantedForest::set_data)
      .method("cross_validation", &RandomPlantedForest::cross_validation)
      .method("predict_matrix", &RandomPlantedForest::predict_matrix)
      .method("predict_vector", &RandomPlantedForest::predict_vector)
      .method("MSE", &RandomPlantedForest::MSE)
      .method("purify", &RandomPlantedForest::purify)
      .method("new_purify", &RandomPlantedForest::new_purify)
      .method("print", &RandomPlantedForest::print)
      .method("get_parameters", &RandomPlantedForest::get_parameters)
      .method("set_parameters", &RandomPlantedForest::set_parameters)
      .method("get_model", &RandomPlantedForest::get_model)
    ;
  
    class_<ClassificationRPF>("ClassificationRPF")
      .derives<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericVector, const NumericMatrix, const String, const NumericVector>()
      .method("set_parameters", &ClassificationRPF::set_parameters)
    ;

    function("testSetComp", &testSetComp);
}
