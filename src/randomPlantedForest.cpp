#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <Rcpp.h>
#include <thread>
#include <assert.h>

using namespace Rcpp;

//  ----------------- functions for converting R and Cpp types ----------------- 

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(R::runif(0,1)*n); }


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
  for(unsigned int row=0; row<v.size(); ++row){
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


//  ----------------- overload of vector operators ----------------- 

// enable to subtract scalar from vector entries
template<typename T, typename S>
std::vector<T> operator-(const std::vector<T>& vec_a, S val){
  std::vector<T> res = vec_a;
  for(auto& entry: res) entry -= val;
  return res;
}

// enable to subtract scalar from vector entries inplace
template<typename T, typename S>
void operator-=( std::vector<T>& vec, S val){
  for(auto& entry: vec) entry -= val;
}

// enable to add scalar to vector entries
template<typename T, typename S>
std::vector<T> operator+(const std::vector<T>& vec_a, S val){
  std::vector<T> res = vec_a;
  for(auto& entry: res) entry += val;
  return res;
}

// enable to add scalar to vector entries inplace
template<typename T, typename S>
void operator+=(std::vector<T>& vec, S val){
  for(auto& entry: vec) entry += val;
}

// enable to multiply vector entries by scalar
template<typename T, typename S>
std::vector<T> operator*(const std::vector<T>& vec_a, S val){
  std::vector<T> res = vec_a;
  for(auto& entry: res) entry *= val;
  return res;
}

// enable to multiply vector entries by scalar
template<typename T, typename S>
std::vector<T> operator*(S val, const std::vector<T>& vec_a){
  std::vector<T> res = vec_a;
  for(auto& entry: res) entry *= val;
  return res;
}

// enable to multiply vector entries by scalar inplace
template<typename T, typename S>
void operator*=(std::vector<T>& vec, S val){
  for(auto& entry: vec) entry *= val;
}

// enable to divide vector entries by scalar
template<typename T, typename S>
std::vector<T> operator/(const std::vector<T>& vec_a, S val){
  std::vector<T> res = vec_a;
  for(auto& entry: res) entry /= val;
  return res;
}

// enable to divide vector entries by scalar inplace
template<typename T, typename S>
void operator/=( std::vector<T>& vec, S val){
  for(auto& entry: vec) entry /= val;
}

// element-wise exp() of vector
template<typename T>
std::vector<T> exp(const std::vector<T>& vec){
  std::vector<T> res = vec;
  for(auto& entry: res) entry = exp(entry);
  return res;
}

// enable to add two vectors
template<typename T>
std::vector<T> operator+(const std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");

  std::vector<T> res = vec_a;
  for(unsigned int i=0;i<vec_b.size(); ++i) res[i] += vec_b[i];
  return res;
}

// enable to add two vectors inplace
template<typename T>
void operator+=(std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");
  
  for(unsigned int i=0;i<vec_b.size(); ++i) vec_a[i] += vec_b[i];
}

// enable to subtract two vectors
template<typename T>
std::vector<T> operator-(const std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");
  
  std::vector<T> res = vec_a;
  for(unsigned int i=0;i<vec_b.size(); ++i) res[i] -= vec_b[i];
  return res;
}

// enable to subtract two vectors inplace
template<typename T>
void operator-=(std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");
  
  for(unsigned int i=0;i<vec_b.size(); ++i) vec_a[i] -= vec_b[i];
}

// enable to multiply two vectors
template<typename T>
std::vector<T> operator*(const std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");

  std::vector<T> res = vec_a;
  for(unsigned int i=0;i<vec_b.size(); ++i) res[i] *= vec_b[i];
  return res;
}

// enable to multiply two vectors inplace
template<typename T>
void operator*=(std::vector<T>& vec_a, const std::vector<T>& vec_b){
  if(vec_a.size() != vec_b.size()) throw std::invalid_argument("The two vectors are not of same size.");
  
  for(unsigned int i=0;i<vec_b.size(); ++i) vec_a[i] *= vec_b[i];
}


// ----------------- custom data types ----------------- 

struct setComp{
  bool operator()(const std::set<int>& a, const std::set<int>& b) const{
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

double eps = 1e-16;

/**
 * \brief Structure to remember split correlated data like the intervals 
 * for each dimension, the individuals that are contained and the corresponding 
 * value per split.
 */
struct Leaf{
    std::vector<int> individuals;       /**< considered samples for each leaf */
    std::vector<double> value;          /**< residual */
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
        DecisionTree(std::set<int> dims, std::vector<Leaf>& first_leaves):
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
      double M_sp;
      double M_bp;
      std::vector<int> I_s;           /**< individuals smaller than splitpoint */
      std::vector<int> I_b;           /**< individuals bigger than splitpoint */
      std::vector<double> M_s;        /**< mean or median of individuals smaller than splitpoin */
      std::vector<double> M_b;        /**< mean or median of individuals bigger than splitpoint */
      std::vector<std::vector<double>> W_s;
      std::vector<std::vector<double>> W_b;
      std::vector<std::vector<double>> Y_s;          
      std::vector<std::vector<double>> Y_b; 
      Split(): min_sum(INF), tree_index(nullptr), leaf_index(nullptr), split_coordinate(1), split_point(0), M_sp(0), M_bp(0) {};
  };

  template<typename T>
  struct Matrix{
    Matrix(){}
    Matrix(std::vector<int> dimensions, T initial_value = 0) 
      : dims(dimensions) {
      for(auto d: dimensions) n_entries *= d;
      entries = std::vector<T>(n_entries, initial_value);
      if(n_entries != entries.size()) throw std::invalid_argument("Invalid matrix size.");
    }
    T& operator[](std::vector<int>& indices){
      if(indices.size() != dims.size()) throw std::invalid_argument("Invalid index.");
      int index = 0;
      for(int i=indices.size()-1; i>0; --i){
        int a = indices[i];
        for(unsigned int d=0; d<dims.size(); ++d){
          if(d<i){
            a *= dims[d];
          }
        }
        index += a;
      }
      index += indices[0];
      if(index < n_entries) throw std::invalid_argument("Invalid index");
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
        for(int i=0; i<dim; ++i){
          d.push_back(i);
        }
        space.push_back(d);
      }
    }
    
    Dimension dimensions;
    std::vector<int> getPoint(){
      std::vector<int> gridPoint(current.size());
      for(unsigned int i=0; i<current.size(); ++i){
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
std::shared_ptr<DecisionTree> treeExists(const std::set<int>& split_dims, TreeFamily& tree_family){
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
bool possibleExists(const int dim, const std::multimap<int, std::shared_ptr<DecisionTree>>& possible_splits, const std::set<int>& resulting_dims){
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
        for(unsigned int i=0; i<intervals.size(); ++i){
            if(leaf.intervals[i] != intervals[i]){
                same_intervals = false;
                break;
            }
        }
        
        if(same_intervals) exists = true;
    }
    return exists;
}

/**
 * \brief Extract keys from a std::map as vector of arbitrary type.
 * 
 * \param map with arbitrary key and value type.
 */
template <typename KT, typename VT>
std::vector<KT> getKeys(std::map<KT, VT, setComp>& m){
  std::vector<KT> keys;
  for(const auto& entry: m){
    keys.push_back(entry.first);
  }
  return keys;
}

template <typename VT>
std::vector<std::vector<VT>> transpose(std::vector<std::vector<VT>>& mat){
  
  if(mat.size()<=0) throw std::invalid_argument("Matrix is empty.");
  int value_size = mat[0].size();
  
  std::vector<std::vector<VT>> columns(value_size);
  for(auto vec: mat){
    for(int n_col=0; n_col<value_size; ++n_col){
      columns[n_col].push_back(vec[n_col]);
    }
  }
  return columns;
}

/**
 * \brief Calculates median of the vector.
 * 
 * \param vec a vector of arbitrary type.
 */
template <typename VT>
VT calcMedian(std::vector<VT>& vec){
  
  // sort vector
  std::sort(vec.begin(), vec.end());
  size_t s = vec.size();
  
  // differ between even and odd case
  if(s % 2 == 0)
    return (vec[s / 2 -1] + vec[s / 2]) / 2;
  return vec[s / 2];
}

/**
 * \brief Calculates row- or columnwise median of a matrix.
 * 
 * \param mat a matrix of arbitrary type.
 */
template <typename VT>
std::vector<VT> calcMedian(std::vector<std::vector<VT>>& mat, bool colwise = true){
  
  std::vector<VT> res;
  
  if(colwise){
    std::vector<std::vector<VT>> columns = transpose(mat);
    
    for(auto vec: columns){
      res.push_back(calcMedian(vec));
    } 
  }else{
    for(auto vec: mat){
      res.push_back(calcMedian(vec));
    } 
  }
  
  return res;
}

/**
 * \brief Calculate mean of a vector.
 * 
 * \param vec a vector of arbitrary type.
 */
template <typename VT>
VT calcMean(std::vector<VT>& vec){
  if(vec.empty()) return 0;
  return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

/**
 * \brief Calculate  row- or columnwise mean of a matrix.
 * 
 * \param mat a matrix of arbitrary type.
 */
template <typename VT>
std::vector<VT> calcMean(std::vector<std::vector<VT>>& mat, bool colwise = true){
  
  if(mat.size() == 0) throw std::invalid_argument("calcMean: Matrix empty - no data provided.");

  int colSize = mat[0].size(), rowSize = mat.size();
  std::vector<VT> res( std::max( colSize, rowSize), 0);
  
  if(colwise){
    res = std::vector<VT>( colSize, 0);
    for(int col=0; col < colSize; ++col){
      for(int row=0; row < rowSize; ++row){
        res[col] += mat[row][col];
      }
      res[col] /= rowSize;
    }
  }else{
    res = std::vector<VT>( rowSize, 0);
    for(int row=0; row < rowSize; ++row){
      for(int col=0; col < colSize; ++col){
        res[row] += mat[row][col];
      }
      res[row] /= colSize;
    } 
  }
  
  return res;
}


// ----------------- main rpf class -----------------

/**
 * \brief Create a prediction model based on Random Forests for regression data sets.
 */
class RandomPlantedForest {
  
    public:
        RandomPlantedForest(const NumericMatrix& samples_Y, const NumericMatrix& samples_X,
                            const NumericVector parameters={1,50,30,10,0.4,0,0,0,0});
        RandomPlantedForest() {};
        void set_data(const NumericMatrix& samples_Y, const NumericMatrix& samples_X);
        NumericMatrix predict_matrix(const NumericMatrix& X, const NumericVector components = {0});
        NumericMatrix predict_vector(const NumericVector& X, const NumericVector components = {0});
        void purify();
        void new_purify();
        void print();
        void cross_validation(int n_sets=4, IntegerVector splits={5,50}, NumericVector t_tries={0.2,0.5,0.7,0.9}, IntegerVector split_tries={1,2,5,10});
        double MSE(const NumericMatrix& Y_predicted, const NumericMatrix& Y_true); 
        void get_parameters();
        void set_parameters(StringVector keys, NumericVector values);
        List get_model();
        
    protected:
        double MSE_vec(const NumericVector& Y_predicted, const NumericVector& Y_true);
        std::vector<std::vector<double>> X;         /**< Nested vector feature samples of size (sample_size x feature_size) */
        std::vector<std::vector<double>> Y;         /**< Corresponding values for the feature samples */          
        int max_interaction;                        /**< Maximum level of interaction determining maximum number of split dimensions for a tree */
        int n_trees;                                /**< Number of trees generated per family */
        int n_splits;                               /**< Number of performed splits for each tree family */
        std::vector<int> n_leaves;                  /**< */
        double t_try = 0.4;                         /**< */
        int split_try = 10;                         /**< */
        int value_size = 1;
        int feature_size = 0;                       /**< Number of feature dimension in X */
        int sample_size = 0;                        /**< Number of samples of X */
        bool purify_forest = 0;                     /**< Whether the forest should be purified */
        bool purified = false;                      /**< Track if forest is currently purified */
        bool deterministic = false;                 /**< Choose whether approach deterministic or random */
        bool parallelize = false;                   /**< Perform algorithm in parallel or serialized */
        bool cross_validate = false;                /**< Determines if cross validation is performed */
        std::vector<double> upper_bounds;           
        std::vector<double> lower_bounds;           
        std::vector<TreeFamily> tree_families;      /**<  random planted forest containing result */
        std::vector<double> predict_single(const std::vector<double> &X, std::set<int> component_index);
        void L2_loss(rpf::Split &split);
        virtual void fit();
        virtual void create_tree_family(std::vector<Leaf> initial_leaves, size_t n);
        virtual rpf::Split calcOptimalSplit(const std::vector<std::vector<double>>& Y, const std::vector<std::vector<double>>& X,
                              std::multimap<int, std::shared_ptr<DecisionTree>>& possible_splits, TreeFamily& curr_family);
};

void RandomPlantedForest::L2_loss(rpf::Split& split){
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);
  split.M_b = calcMean(split.Y_b);
  
  for(unsigned int p=0; p<split.Y_s[0].size(); ++p){
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      split.min_sum += pow(split.Y_s[individual][p] - split.M_s[p], 2) - pow(split.Y_s[individual][p], 2);
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      split.min_sum += pow(split.Y_b[individual][p] - split.M_b[p], 2) - pow(split.Y_b[individual][p], 2);
    }
  }
}

// constructor
RandomPlantedForest::RandomPlantedForest(const NumericMatrix& samples_Y, const NumericMatrix& samples_X,
                                         const NumericVector parameters){

    // Ensure correct Rcpp RNG state
    Rcpp::RNGScope scope; 
  
    // initialize class members 
    std::vector<double> pars = to_std_vec(parameters);
    if(pars.size() != 9){
        Rcout << "Wrong number of parameters - set to default." << std::endl;
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
rpf::Split RandomPlantedForest::calcOptimalSplit(const std::vector<std::vector<double>>& Y, const std::vector<std::vector<double>>& X,
                                               std::multimap<int, std::shared_ptr<DecisionTree>>& possible_splits, TreeFamily& curr_family){

  rpf::Split curr_split, min_split;
  std::set<int> tree_dims;
  int k;
  bool splitable;
  size_t n = 0;
  double leaf_size, sample_point;
  
  // sample possible splits
  unsigned int n_candidates = ceil(t_try * possible_splits.size()); // number of candidates that will be considered
  std::vector<int> split_candidates(possible_splits.size());
  std::iota(split_candidates.begin(), split_candidates.end(), 0); // consecutive indices of possible candidates
  
  if(!deterministic){
    std::random_shuffle(split_candidates.begin(), split_candidates.end(), randWrapper); // shuffle for random order
  }
  
  // consider a fraction of possible splits
  int iter = 100;
  while(n < n_candidates){
    
    iter--;
    
    // in the beginning not known if split viable
    splitable = false;
    
    // since size of possible splits changes, check if candidate in range
    if(possible_splits.empty() || iter == 0) break;
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
        std::vector<double> unique_samples(curr_individuals.size());
        for(unsigned int i=0; i<curr_individuals.size(); ++i){
          unique_samples[i] = X[curr_individuals[i]][k];
        }
        std::sort(unique_samples.begin(), unique_samples.end());
        unique_samples.erase(std::unique(unique_samples.begin(), unique_samples.end()), unique_samples.end());
        
        // check if number of sample points is within limit
        if(unique_samples.size() < 2*leaf_size) continue;
        splitable = true;
        
        int start = 0;
        int end = std::min(split_try, int(unique_samples.size() - 1));
        if(deterministic){
          start = 1;
          end = unique_samples.size() - 1;
        }
        
        // consider split_try-number of random samples
        for(int t = start; t<end; ++t){
          
          // get samplepoint
          auto sample_pos = unique_samples.begin();
          std::uniform_int_distribution<> distrib(leaf_size, unique_samples.size() - leaf_size);
          if(deterministic){
            std::advance(sample_pos, t);
          }else{
            int offset = R::runif(leaf_size, unique_samples.size() - leaf_size );
            std::advance(sample_pos, offset); // consider only sample points with offset
          }
          sample_point = *sample_pos;
          
          // clear current split
          {
            curr_split.I_s.clear();
            curr_split.I_b.clear();
            curr_split.Y_s.clear();
            curr_split.Y_b.clear();
            curr_split.I_s.reserve(curr_individuals.size());
            curr_split.I_b.reserve(curr_individuals.size());
            curr_split.Y_s.reserve(curr_individuals.size());
            curr_split.Y_b.reserve(curr_individuals.size());
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
            min_split = curr_split;
            min_split.tree_index = curr_tree;
            min_split.leaf_index =  &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point = sample_point;
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

void RandomPlantedForest::set_data(const NumericMatrix& samples_Y, const NumericMatrix& samples_X){

    this->Y = to_std_vec(samples_Y);
    this->X = to_std_vec(samples_X);
    
    // Check for correct input
    if(Y.size() == 0) throw std::invalid_argument("Y empty - no data provided.");
    if(X.size() == 0) throw std::invalid_argument("X empty - no data provided.");
    this->feature_size = X[0].size();
    this->value_size = Y[0].size(); // multiclass
    for(const auto &vec:X){
      if(vec.size() != feature_size) throw std::invalid_argument("Feature dimensions of X not uniform.");
    }
    if(Y.size() != X.size()) throw std::invalid_argument("X and Y are not of the same length!");
    
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
        this->upper_bounds[i] = maxVal + 2 * eps; // to consider samples at max value
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
    Rcout << "Initial Possible Splits: ";
    for(auto split: possible_splits){
      Rcout << split.first << "-";
      for(auto dim: split.second->split_dims) Rcout << dim << ",";
      Rcout << "; ";
    }
    Rcout << std::endl;
  }
  
  // sample data points with replacement
  int sample_index;
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;
  
  // deterministic
  if(deterministic){
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  }else{
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);
    
    for(size_t i=0; i<sample_size; ++i){
      sample_index = R::runif(0, sample_size - 1);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }
  
  // modify existing or add new trees through splitting
  rpf::Split curr_split;
  for(size_t split_count=0; split_count<n_splits; ++split_count){
    
    // find optimal split
    curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family);
    
    // continue only if we get a significant result
    if(!std::isinf(curr_split.min_sum)){
      
      if(false){
        Rcout << "Current Optimal Split: " << curr_split.min_sum << "; " << curr_split.split_coordinate << "- ";
        for(auto dim: curr_split.tree_index->split_dims) Rcout << dim << ", ";
        //Rcout << "; " << curr_split.I_s.size() << "/" << curr_split.I_b.size() << "=" << curr_split.I_s.size()+curr_split.I_b.size() << "; " <<
        //  curr_split.M_s << "/" << curr_split.M_b <<  std::endl;
      }

      // update possible splits
      for(int feature_dim=1; feature_dim<=feature_size; ++feature_dim){ // consider all possible dimensions
        
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
          Rcout << "Updated Possible Splits: " << std::endl;
          for(auto split: possible_splits){
            Rcout << split.first << "-";
            for(auto dim: split.second->split_dims) Rcout << dim << ",";
            Rcout << "; ";
          }
          Rcout << std::endl;
        }
      }
      
      // update values of individuals of split interval with mean
      for(int individual: curr_split.leaf_index->individuals){ // todo: loop directly over I_s I_b
        if(samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point){
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
        Rcout << "First leaf: intervals=" << std::endl;
        for(auto interval: leaf_s.intervals) Rcout << interval.first << "," << interval.second << ";";
        Rcout << "individuals=";
        for(auto i: leaf_s.individuals) Rcout << i << ",";
        Rcout << std::endl << "Second leaf: intervals=" << std::endl;
        for(auto interval: leaf_b.intervals) Rcout << interval.first << "," << interval.second << ";";
        Rcout << "individuals=";
        for(auto i: leaf_b.individuals) Rcout << i << ",";
        Rcout << std::endl;
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
        Rcout << "Current TreeFamily: (" << curr_family.size() << ") ";
        for(auto tree: curr_family){
          Rcout << "Dims = ";
          for(auto dim: tree.first) Rcout << dim << ", ";
          Rcout << "; " << "Number of Leaves = " << tree.second->leaves.size();
          Rcout << " / ";
        }
        Rcout << std::endl << std::endl;
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
        initial_leaf.value = std::vector<double>(value_size, 0);
        initial_leaf.individuals = initial_individuals;
        initial_leaf.intervals = initial_intervals;
    }
    std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf
    
    // initialize tree families
    this->tree_families = std::vector<TreeFamily>(n_trees); 

    // iterate over families of trees and modify
    if(parallelize){
      int n_threads = std::thread::hardware_concurrency() - 1;
      for(int n = 0; n<n_trees; n+=n_threads){
        if(n >= (n_trees - n_threads)) n_threads = n_trees % n_threads;
        std::vector<std::thread> threads(n_threads);
        for(size_t t=0; t<n_threads; ++t){
          Rcout << (n + t) << "/" << n_trees << std::endl;
          threads[t] = std::thread(&RandomPlantedForest::create_tree_family, this, std::ref(initial_leaves), n + t);
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
    
    /*
    bool cv_tmp = this->cross_validate;
    this->cross_validate = false;
  
    if(deterministic) {
      Rcout << "Note: Set model to non-deterministic. " << std::endl;
      deterministic = false; 
    }
    
    std::set<int> splits_vec = to_std_set(splits);
    std::vector<int> split_tries_vec = to_std_vec(split_tries);
    std::vector<double> t_tries_vec = to_std_vec(t_tries);
    
    if(splits_vec.size()!=2) {Rcout << "Min and max needed for number of splits." << std::endl; return;}
    
    // remember optimal parameter set and MSE
    double  MSE_sum = 0, curr_MSE = 0, MSE_min = INF, optimal_split = INF, optimal_t_try = INF, optimal_split_try = INF;
    int optimal_inter = 1;
    
    std::vector<int> order(sample_size);
    std::iota(order.begin(), order.end(), 0);
    std::random_shuffle(order.begin(), order.end(), randWrapper);
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
          			    Rcout << inter << ", " << splits << ", " << t << ", " << s << ": MSE=" << curr_MSE << std::endl;
          			    
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
    
    Rcout << "Optimal parameters: " << optimal_inter << ", " << optimal_split << ", " << optimal_t_try << ", " << optimal_split_try << ": MSE=" << MSE_min << std::endl;
    */
}

// predict single feature vector
std::vector<double> RandomPlantedForest::predict_single(const std::vector<double>& X, std::set<int> component_index){

    std::vector<double> total_res = std::vector<double>(value_size, 0);
    
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
                    for(unsigned int i = 0; i<dims.size(); ++i){

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

    return total_res / n_trees;
}

// predict multiple feature vectors
Rcpp::NumericMatrix RandomPlantedForest::predict_matrix(const NumericMatrix& X, const NumericVector components){
    std::vector<std::vector<double>> feature_vec = to_std_vec(X);
    std::set<int> component_index = to_std_set(components);
    std::vector<std::vector<double>> predictions;
    
    // todo: sanity check for X
    if( feature_vec.empty() ) throw std::invalid_argument("Feature vector is empty.");
    if( component_index == std::set<int>{0} && feature_vec[0].size() != this->feature_size ) throw std::invalid_argument("Feature vector has wrong dimension.");
    if( component_index != std::set<int>{0} && component_index.size() != feature_vec[0].size() ) throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");
    
    for(auto& vec: feature_vec){
        predictions.push_back(predict_single(vec, component_index));
    }

    return from_std_vec(predictions);
}

Rcpp::NumericMatrix RandomPlantedForest::predict_vector(const NumericVector& X, const NumericVector components){
    std::vector<double> feature_vec = to_std_vec(X);
    std::set<int> component_index = to_std_set(components);
    std::vector<std::vector<double>> predictions;
    Rcpp::NumericMatrix res;
    
    // todo: sanity check for X
    if(feature_vec.empty()) { Rcout << "Feature vector is empty." << std::endl; return res; }
    
    if(component_index == std::set<int>{0} && feature_vec.size() != this->feature_size){ 
        Rcout << "Feature vector has wrong dimension." << std::endl; return res;
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

double RandomPlantedForest::MSE_vec(const NumericVector& Y_predicted, const NumericVector& Y_true){
    return sum(Rcpp::pow(Y_true - Y_predicted, 2)) / Y_true.size();
}

double RandomPlantedForest::MSE(const NumericMatrix& Y_predicted, const NumericMatrix& Y_true){
    // todo: multiclass
    double sum = 0;
    int Y_size = Y_predicted.size();
    
    for(int i=0; i<Y_size; ++i){
      sum += MSE_vec(Y_predicted( i , _ ), Y_true( i , _ ));
    }

    return sum / Y_size;
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
                for(int i=0; i<value_size; ++i) new_leaf.value[i] = -curr_leaf.value[i] * multiplier; // update value of new leaf
                
                // append new leaf
                if(!leafExists(new_leaf.intervals, curr_tree)) curr_tree->leaves.push_back(new_leaf);
                for(int i=0; i<value_size; ++i) new_leaf.value[i] = curr_leaf.value[i] * multiplier; // update value of new leaf
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
  
  /*
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
    double epsilon = 0.01;
    bounds.push_back(bounds[bounds.size() - 1] + epsilon); // toDo: 
    lim_list[curr_dim - 1] = bounds;
  }
  
  // initialize values and individuals for each tree in family
  std::vector< NDGrid::NDGrid> grids(curr_family.size() - 1); // ignore null tree?
  std::vector<rpf::Matrix<int>> individuals(curr_family.size() - 1);
  std::vector<rpf::Matrix<double>> values(curr_family.size() - 1);
  std::vector<std::set<int>> variables(curr_family.size() - 1); // todo: remove after checking if order of insertion correct
  
  if(false){
    // show original leafes
    for(const auto& curr_tree: curr_family){
      for(auto dim: curr_tree.first) Rcout << dim << ", "; // << std::endl;
      Rcout << ": " << std::endl;
      for(const auto& leaf: curr_tree.second->get_leaves()){
        Rcout << leaf.individuals.size() << ", ";
      }
      Rcout << std::endl;
    }
  }
  
  //  ------------- setup finer grid  -------------
  
  int tree_index = 0;
  for(const auto& curr_tree: curr_family){
    
    if(false){
      Rcout << tree_index << ": ";
      for(auto dim: curr_tree.first) Rcout << dim << ", "; // << std::endl;
      Rcout  << std::endl;
    }
    
    if(false){
      Rcout << "Lim-list: " << std::endl;;
      for(auto dim: curr_tree.first){
        for(auto l: lim_list[dim - 1]){
          Rcout << l << ", ";
        }
        Rcout  << std::endl;
      }
      Rcout  << std::endl;
    }
    
    if(curr_tree.first == std::set<int>{0}) continue; // ignore null tree?
    
    // fill space with dimensions
    std::vector<int> dimensions;
    for(const auto& dim: curr_tree.first){
      dimensions.push_back(lim_list[dim - 1].size() - 1); // size - 1 ?
    }
    
    // setup grid for leaf indices
    auto grid =  NDGrid::NDGrid(dimensions);
    
    // initialize data for current tree
    grids[tree_index] = grid;
    individuals[tree_index] = rpf::Matrix<int>(dimensions);
    values[tree_index] = rpf::Matrix<double>(dimensions);
    variables[tree_index] = curr_tree.first;
    
    //if(curr_tree.first.size() != 4) continue; // todo: remove
    
    if(false){
      Rcout << "Matrix dimensions: ";
      for(const auto& d: dimensions) Rcout <<  d << ", ";
      Rcout << std::endl;
    }
    
    // fill grid points with individuals and values
    while(!grid.nextPoint()){
      
      std::vector<int> gridPoint = grid.getPoint();
      
      if(false){
        Rcout << "Gridpoint: ";
        for(const auto& p: gridPoint) Rcout << p << ", ";
        //Rcout << std::endl;
      }
      
      bool in_leaf = true;
      
      // go through sample points to sum up individuals
      int i = 0;
      for(const auto& feature_vec: X){
        int dim_index = 0;
        //Rcout << "Datapoint " << i << ": ";
        in_leaf = true;
        for(const auto& dim: curr_tree.first){
          double val = feature_vec[dim - 1];
          //Rcout << val << " in (" << lim_list[dim - 1][gridPoint[dim_index]] << "," << lim_list[dim - 1][gridPoint[dim_index] + 1] << "), ";
          if( !( (val >= lim_list[dim - 1][gridPoint[dim_index]]) 
                   && (val < lim_list[dim - 1][gridPoint[dim_index] + 1]) ) ) in_leaf = false;
          ++dim_index;
        }
        //Rcout << in_leaf << " / ";
        
        // consider individuals only if all in
        if(in_leaf) individuals[tree_index][gridPoint] += 1;
        
        ++i;
      }
      //Rcout << std::endl << std::endl;
      
      // go through leaves of tree to sum up values
      for(const auto& leaf: curr_tree.second->get_leaves()){
        
        in_leaf = true;
        int dim_index = 0;
        for(const auto& dim: curr_tree.first){
          
          // consider values only if all in
          if( !( (leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) 
                   && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1]) ) ) in_leaf = false;
          
          ++dim_index;
        }
        
        // sum up values
        if(in_leaf) values[tree_index][gridPoint] += leaf.value; // todo: multiclass
      }
    }
    
    if(false){
      int sum = 0;
      while(!grid.nextPoint()){
        sum += individuals[tree_index][grid.getPoint()];
      }
      Rcout << tree_index << ": " << sum << ", " << std::endl;
    }
    
    ++tree_index;
  }
  
  // ------------- create new trees -------------
  
  // insert null tree
  grids.insert(grids.begin(), NDGrid::NDGrid());
  values.insert(values.begin(), rpf::Matrix<double>(std::vector<int>{1}));
  individuals.insert(individuals.begin(), rpf::Matrix<int>(std::vector<int>{1}));
  variables.insert(variables.begin(), std::set<int>{0});
  
  // recap maximum number of dimensions of current family
  unsigned int curr_max = 0;
  for(const auto& tree: curr_family){
    if(tree.first.size() > curr_max) curr_max = tree.first.size();
  }
  
  auto keys = getKeys(curr_family);
  while(curr_max > 1){
    
    // go through split dimensions of all trees
    for(std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key){
      
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
            if(!tree){
              
              // get index of old and new tree
              auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
              curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
              auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
              
              // remove matrix dimension of respective variable
              std::vector<int> matrix_dimensions = values[old_tree_index].dims;
              matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);
              
              if(false){
                std::cout << "   New tree: ";
                for(const auto& dim: tree_dims) std::cout << dim << ", "; 
                std::cout << std::endl <<  "   Removing dimension " << dim_index << " from {";
                for(auto dim: values[old_tree_index].dims) std::cout << dim << ", "; 
                std::cout <<  "} results in: {";
                for(auto dim: matrix_dimensions) std::cout << dim << ", "; 
                std::cout << "}" << std::endl;
              }
              
              // initialize data for new tree
              auto grid =  NDGrid::NDGrid(matrix_dimensions);
              grids.insert(grids.begin() + tree_index, grid);
              values.insert(values.begin() + tree_index, rpf::Matrix<double>(matrix_dimensions));
              individuals.insert(individuals.begin() + tree_index, rpf::Matrix<int>(matrix_dimensions));
              variables.insert(variables.begin() + tree_index, tree_dims);
              
              // fill individuals of new trees
              while(!grid.nextPoint()){
                
                std::vector<int> gridPoint = grid.getPoint();
                bool in_leaf = true;
                
                // go through sample points to sum up individuals
                for(const auto& feature_vec: X){
                  int dim_index = 0;
                  in_leaf = true;
                  for(const auto& dim: tree_dims){
                    double val = feature_vec[dim - 1];
                    if( !( (val >= lim_list[dim - 1][gridPoint[dim_index]]) 
                             && (val < lim_list[dim - 1][gridPoint[dim_index] + 1]) ) ) in_leaf = false;
                    ++dim_index;
                  }
                  
                  // consider individuals only if all in
                  if(in_leaf) individuals[tree_index][gridPoint] += 1;
                }
              }
              
              if(false){
                while(!grid.nextPoint()){
                  std::cout << individuals[tree_index][grid.getPoint()] << ", ";
                }
                std::cout << std::endl;
              }
            }
            
            dim_index++;
          }
        }
      }
    }
    
    // update currently considered dimension size
    --curr_max;
  }
  
  
  // ------------- purify -------------
  
  // measure tolerance and number of iterations
  std::vector<double> tol(curr_family.size(), 1);
  int iter;
  
  // iterate backwards through tree family
  int curr_tree_index = curr_family.size() - 1;
  for(TreeFamily::reverse_iterator curr_tree = curr_family.rbegin(); curr_tree != curr_family.rend(); ++curr_tree ) { 
    iter = 0;
    std::set<int> curr_dims = curr_tree->second->get_split_dims();
     
    // do not purify null
    if(curr_dims == std::set<int>{0}) continue;
    
    // repeat until tolerance small enough and (?) maximum number of iterations reached
    while( (tol[curr_tree_index] > 0.00000000001) && (iter < 100) ){
     
      // go through feature dims
      int curr_dim_index = 0;
      for(const auto& feature_dim: curr_dims){
        
        // get tree that has same variables as curr_tree minus j-variable
        std::set<int> tree_dims = curr_dims;
        tree_dims.erase(tree_dims.find(feature_dim));
        int tree_index = 0; // if tree not exist, set to null tree 
        if(curr_family.find(tree_dims) != curr_family.end()) tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims)) - 1;
        
        // update values
        if(grids[curr_tree_index].dimensions.size() == 1){ // one dimensional case
          
          int sum_ind = 0, avg = 0;
         
          // get sum of individuals
          for(int i=0; i<individuals[curr_tree_index].n_entries; ++i) sum_ind += individuals[curr_tree_index][{i}];
          if(sum_ind == 0) continue;
         
          // calc avg
          for(int i=0; i<individuals[curr_tree_index].n_entries; ++i) avg += (individuals[curr_tree_index][{i}] * values[curr_tree_index][{i}]) / sum_ind;
         
          // update values of one dimensional and null tree
          for(int i=0; i<values[curr_tree_index].n_entries; ++i) values[curr_tree_index][{i}] -= avg;
          values[tree_index][{0}] += avg;
         
        }else{ // higher dimensional case
         
          // setup new grid without dimension j
          std::vector<int> new_dimensions = grids[curr_tree_index].dimensions;
          int j_dim = new_dimensions[curr_dim_index];
          new_dimensions.erase(new_dimensions.begin() + curr_dim_index);
          NDGrid::NDGrid grid =  NDGrid::NDGrid(new_dimensions);
          
          // go through values without dimension j
          while(!grid.nextPoint()){
           
            auto gridPoint = grid.getPoint();
            
            int sum_ind = 0, avg = 0;
            
            // go through slice to sum up individuals
            for(int j = 0; j<j_dim; ++j){
              gridPoint = grid.getPoint();
              gridPoint.push_back(j);
               
              // get sum of individuals
              sum_ind += individuals[curr_tree_index][gridPoint];
            }
            
            // go through slice to calc avg
            for(int j = 0; j<j_dim; ++j){
              gridPoint = grid.getPoint();
              gridPoint.push_back(j);
               
              // calc avg
              avg += (individuals[curr_tree_index][gridPoint] * values[curr_tree_index][gridPoint]) / sum_ind;
            }
            
            // go through slice to update values
            for(int j = 0; j<j_dim; ++j){
              gridPoint = grid.getPoint();
              gridPoint.push_back(j);
               
              // update values of current slice
              values[curr_tree_index][gridPoint] -= avg;
            }
              
            // update lower dimensional tree
            values[tree_index][grid.getPoint()] += avg;
          }
        }
        
        ++curr_dim_index;
      }
      
      // update tolerance
      if(variables[curr_tree_index].size() == 1){
        tol[curr_tree_index] = 1; // todo
      }else{
        tol[curr_tree_index] = 1;
      }
      
      ++iter;
    }
    
    --curr_tree_index;
  } 
  
    
  // ------------- todo: convert to rpf class -------------
  
  /*
  for(int tree_index=0; i<variables.size(); ++tree_index){

    NDGrid::NDGrid grid = grids[tree_index];
    DecisionTree& curr_tree = curr_family.find(variables[tree_index]);
    

    while(!grid.nextPoint()){
      
      std::vector<int> gridPoint = grid.getPoint();
      
      Leaf new_leaf;
      new_leaf.value = values[tree_index][gridPoint]; // idea: save values[tree_index] separately to save access time
      curr_tree
    }
  }
  */
  
  //}
}

void RandomPlantedForest::print(){
    for(int n=0; n<n_trees; ++n){
        TreeFamily family = tree_families[n];
        // Rcout << n+1 << " TreeFamily: constant=" << family.constant << std::endl << std::endl;
        auto keys = getKeys(family);
        for(int m=0; m<keys.size(); ++m){
            DecisionTree tree = *(family[keys[m]]);
            Rcout << m+1 << " Tree: ";
            Rcout << "Dims=";
            for(const auto& dim: tree.split_dims) Rcout << dim << ",";
            Rcout << std::endl << "Leaves: (" << tree.leaves.size() << ")" << std::endl;
            for(const auto& leaf: tree.leaves){
                Rcout << "Intervals=";
                for(const auto& interval: leaf.intervals){
                    Rcout << interval.first << "," << interval.second << "/";
                }
                Rcout << " Value=";
                for(const auto& val: leaf.value) Rcout << val << ", ";
                Rcout << std::endl;
            }  
            Rcout << std::endl;
        }
        Rcout << std::endl << std::endl;
    }
}

// print parameters of the model to the console
void RandomPlantedForest::get_parameters(){
  Rcout << "Parameters: n_trees=" <<  n_trees << ", n_splits=" << n_splits << ", max_interaction=" << max_interaction << ", t_try=" << t_try 
            << ", split_try=" << split_try << ", purified=" << purified << ", deterministic=" << deterministic << ", parallel=" << parallelize
            << ", feature_size=" << feature_size << ", sample_size=" << sample_size << std::endl;
}

/*  retrospectively change parameters of existing class object, 
    updates the model, so far only single valued parameters supported,
    for replacing training data use 'set_data', 
    note that changing cv does not trigger cross validation */
void RandomPlantedForest::set_parameters(StringVector keys, NumericVector values){
  if(keys.size() != values.size()) {
    Rcout << "Size of input vectors is not the same. " << std::endl;
    return;
  }
  
  for(unsigned int i=0; i<keys.size(); ++i){
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
      Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}

List RandomPlantedForest::get_model(){
  List model;
  for(const auto family: tree_families){
    List variables, family_values, family_intervals;
    for(const auto tree: family){
      List tree_values;
      List tree_intervals;
      variables.push_back(from_std_set(tree.first));
      for(const auto leaf: tree.second->leaves){
        NumericMatrix leaf_values;
        for(const auto& val: leaf.value){
          leaf_values.push_back(val);
        }
        tree_values.push_back(leaf_values);
        
        NumericVector intervals;
        for(const auto& interval: leaf.intervals){
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

  public:
    ClassificationRPF(const NumericMatrix& samples_Y, const NumericMatrix& samples_X,
                      const String loss="L2", const NumericVector parameters={1,50,30,10,0.4,0,0,0,0,0,0.1});
    void set_parameters(StringVector keys, NumericVector values);
    
  private:
    double delta;
    double epsilon;
    enum LossType { L1, L2, median, logit, logit_2, logit_3, exponential, exponential_2, exponential_3};
    LossType loss; 
    void (ClassificationRPF::*calcLoss)(rpf::Split&);  
    void create_tree_family(std::vector<Leaf> initial_leaves, size_t n) override;
    void fit() override;
    rpf::Split calcOptimalSplit(const std::vector<std::vector<double>>& Y, const std::vector<std::vector<double>>& X,
                               std::multimap<int, std::shared_ptr<DecisionTree>>& possible_splits, TreeFamily& curr_family, 
                               std::vector<std::vector<double>>& weights);
    void L1_loss(rpf::Split &split);
    void median_loss(rpf::Split &split);
    void logit_loss(rpf::Split &split);
    void logit_loss_2(rpf::Split &split);
    void logit_loss_3(rpf::Split &split);
    void exponential_loss(rpf::Split &split);
    void exponential_loss_2(rpf::Split &split);
    void exponential_loss_3(rpf::Split &split);
};

void ClassificationRPF::L1_loss(rpf::Split& split){
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);
  split.M_b = calcMean(split.Y_b); 
  
  for(unsigned int p=0; p<split.Y_s[0].size(); ++p){
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      split.min_sum += abs(split.Y_s[individual][p] - split.M_s[p]) - abs(split.Y_s[individual][p]); 
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      split.min_sum += abs(split.Y_b[individual][p] - split.M_b[p]) - abs(split.Y_b[individual][p]);
    }
  }
}

void ClassificationRPF::median_loss(rpf::Split& split){
  split.min_sum = 0;
  split.M_s = calcMedian(split.Y_s);
  split.M_b = calcMedian(split.Y_b); 
  
  for(unsigned int p=0; p<split.Y_s[0].size(); ++p){
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      split.min_sum += abs(split.Y_s[individual][p] - split.M_s[p]) - abs(split.Y_s[individual][p]);
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      split.min_sum += abs(split.Y_b[individual][p] - split.M_b[p]) - abs(split.Y_b[individual][p]); 
    }
  }
}
  
void ClassificationRPF::logit_loss(rpf::Split& split){
    
    split.min_sum = 0;
    split.M_s = calcMean(split.Y_s);  
    split.M_b = calcMean(split.Y_b);
    split.M_sp = 1 - std::accumulate(split.M_s.begin(),split.M_s.end(), 0.0);
    split.M_bp = 1 - std::accumulate(split.M_b.begin(),split.M_b.end(), 0.0);

    std::vector<double> M_s = split.M_s;
    std::vector<double> M_b = split.M_b;

    std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
    std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
    
    double M_sp = std::max(delta, split.M_sp);
    double M_bp = std::max(delta, split.M_bp);

    std::vector<double> W_s_mean = calcMean(split.W_s);
    std::vector<double> W_b_mean = calcMean(split.W_b);

    std::vector<std::vector<double>> W_s = split.W_s, W_b = split.W_b, W_s_new = split.W_s, W_b_new = split.W_b;  
    
    for(unsigned int p=0; p<W_s[0].size(); ++p){ 
      for(unsigned int individual=0; individual<W_s.size(); ++individual){
        W_s[individual][p] = exp(W_s[individual][p]); 
        W_s_new[individual][p] = exp(W_s_new[individual][p] + log(M_s[p] / M_sp) - W_s_mean[p]); 
      }
      for(unsigned int individual=0; individual<W_b.size(); ++individual){
        
        W_b[individual][p] = exp(W_b[individual][p]); 
        W_b_new[individual][p] = exp(W_b_new[individual][p] + log(M_b[p] / M_bp) - W_b_mean[p]);
      } 
    }   

    for(unsigned int p=0; p<W_b[0].size(); ++p){
      for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
        split.min_sum += split.Y_s[individual][p] * log(W_s[individual][p] / (1 + std::accumulate(W_s[individual].begin(),W_s[individual].end(),0.0)) ); // ~ R_old
        split.min_sum -= split.Y_s[individual][p] * log(W_s_new[individual][p] / (1 + std::accumulate(W_s_new[individual].begin(),W_s_new[individual].end(),0.0)) ); // ~ R_new
      }
      for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
        split.min_sum += split.Y_b[individual][p] * log(W_b[individual][p] / (1 + std::accumulate(W_b[individual].begin(),W_b[individual].end(),0.0)) ); // ~ R_old
        split.min_sum -= split.Y_b[individual][p] * log(W_b_new[individual][p] / (1 + std::accumulate(W_b_new[individual].begin(),W_b_new[individual].end(),0.0)) ); // ~ R_new
      }
    }     

    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      split.min_sum += (1 - std::accumulate(split.Y_s[individual].begin(), split.Y_s[individual].end(), 0.0) * log(1 / (1 + std::accumulate(W_s[individual].begin(), W_s[individual].end(), 0.0)) )); // ~ R_old
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      split.min_sum += (1 - std::accumulate(split.Y_b[individual].begin(), split.Y_b[individual].end(), 0.0) * log(1 / (1 + std::accumulate(W_b[individual].begin(), W_b[individual].end(), 0.0)) )); // ~ R_old
    }
    
    if(std::isnan(split.min_sum)){
      split.min_sum = INF;
    }
  }

void ClassificationRPF::logit_loss_2(rpf::Split& split){
  
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);  
  split.M_b = calcMean(split.Y_b);
  
  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;
  
  std::vector<double> M_s2 = split.M_s;
  std::vector<double> M_b2 = split.M_b;
  
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
  
  std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(delta, 1-M); });
  std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(delta, 1-M); });
  
  std::vector<double> W_s_mean = calcMean(split.W_s);
  std::vector<double> W_b_mean = calcMean(split.W_b);
  
  std::vector<std::vector<double>> W_s = split.W_s, W_b = split.W_b, W_s_new = split.W_s, W_b_new = split.W_b;  
  
  for(int p=0; p<value_size; ++p){ 
    for(unsigned int individual=0; individual<W_s.size(); ++individual){
      W_s[individual][p] = exp(W_s[individual][p]); 
      W_s_new[individual][p] = exp(W_s_new[individual][p] + log(M_s[p] / M_s2[p]) - W_s_mean[p]); 
    }
    for(unsigned int individual=0; individual<W_b.size(); ++individual){
      
      W_b[individual][p] = exp(W_b[individual][p]); 
      W_b_new[individual][p] = exp(W_b_new[individual][p] + log(M_b[p] / M_b2[p]) - W_b_mean[p]); 
    } 
  }
  
  for(int p=0; p<value_size; ++p){
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      split.min_sum += split.Y_s[individual][p] * log(W_s[individual][p] / (1 + W_s[individual][p] )); // ~ R_old
      split.min_sum -= split.Y_s[individual][p] * log(W_s_new[individual][p] / (1 + W_s_new[individual][p] )); // ~ R_new
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      split.min_sum += split.Y_b[individual][p] * log(W_b[individual][p] / (1 + W_b[individual][p] )); // ~ R_old
      split.min_sum -= split.Y_b[individual][p] * log(W_b_new[individual][p] / (1 + W_b_new[individual][p] )); // ~ R_new
    }
  }
  
  if(std::isnan(split.min_sum)){
    split.min_sum = INF;
  } 
}

void ClassificationRPF::logit_loss_3(rpf::Split& split){
  
  split.min_sum = 0;
  split.M_s = calcMean(split.Y_s);  
  split.M_b = calcMean(split.Y_b);
  split.M_sp = 1 - std::accumulate(split.M_s.begin(),split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(),split.M_b.end(), 0.0); 
 
  std::vector<double> M_s = split.M_s;
  std::vector<double> M_b = split.M_b;
  
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(delta, M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(delta, M); });
  
  std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = log(M); });
  std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = log(M); });

  double M_sp = std::max(delta, split.M_sp);
  double M_bp = std::max(delta, split.M_bp);    
  
  M_sp = log(M_sp);
  M_bp = log(M_bp);

  double sum_s = (std::accumulate(M_s.begin(), M_s.end(),0.0)+M_sp)/(M_s.size()+1);
  double sum_b = (std::accumulate(M_b.begin(), M_b.end(),0.0)+M_bp)/(M_b.size()+1);

  std::vector<double> W_s_mean = calcMean(split.W_s);
  std::vector<double> W_b_mean = calcMean(split.W_b);
  
  std::vector<std::vector<double>> W_s = split.W_s, W_b = split.W_b, W_s_new = split.W_s, W_b_new = split.W_b;
  
  std::vector<std::vector<double>> Y_s = split.Y_s;
  std::vector<std::vector<double>> Y_b = split.Y_b;
  
  for(unsigned int p=0; p<W_s[0].size(); ++p){ 
    for(unsigned int individual=0; individual<W_s.size(); ++individual){
      
      W_s_new[individual][p] = W_s_new[individual][p] + M_s[p] - sum_s - W_s_mean[p]; 
    }
    for(unsigned int individual=0; individual<W_b.size(); ++individual){
      
      W_b_new[individual][p] = W_b_new[individual][p] + M_b[p] - sum_b - W_b_mean[p]; 
    } 
  }
  
  std::vector<double> W_sp;
  std::vector<double> W_bp;
  std::vector<double> W_sp_new;
  std::vector<double> W_bp_new;
  
  std::vector<double> Y_sp;
  std::vector<double> Y_bp;
  
  for(unsigned int individual=0; individual<W_s.size(); ++individual){
    
    W_sp.push_back(- accumulate(W_s[individual].begin(), W_s[individual].end(), 0.0)); 
    W_sp_new.push_back(- accumulate(W_s_new[individual].begin(), W_s_new[individual].end(), 0.0));
    
    Y_sp.push_back(1-accumulate(Y_s[individual].begin(), Y_s[individual].end(),0.0));
  }
  
  for(unsigned int individual=0; individual<W_b.size(); ++individual){
    
    W_bp.push_back(- accumulate(W_b[individual].begin(), W_b[individual].end(), 0.0)); 
    W_bp_new.push_back(- accumulate(W_b_new[individual].begin(), W_b_new[individual].end(), 0.0));
    
    Y_bp.push_back(1-accumulate(Y_b[individual].begin(), Y_b[individual].end(),0.0));
  }
  
  W_s = transpose(W_s);
  W_s.push_back(W_sp);
  W_s = transpose(W_s);
  
  W_b = transpose(W_b);
  W_b.push_back(W_bp);
  W_b = transpose(W_b);
  
  W_s_new = transpose(W_s_new);
  W_s_new.push_back(W_sp_new);
  W_s_new = transpose(W_s_new);
  
  W_b_new = transpose(W_b_new);
  W_b_new.push_back(W_bp_new);
  W_b_new = transpose(W_b_new);
  
  Y_s=transpose(Y_s);
  Y_s.push_back(Y_sp);
  Y_s = transpose(Y_s);
  
  Y_b = transpose(Y_b);
  Y_b.push_back(Y_bp);
  Y_b = transpose(Y_b);
  
  for(unsigned int p=0; p<W_s[0].size(); ++p){ 
    for(unsigned int individual=0; individual<W_s.size(); ++individual){
      
      W_s[individual][p] = exp(W_s[individual][p]);         
      W_s_new[individual][p] = exp(W_s_new[individual][p]); 
    }
    for(unsigned int individual=0; individual<W_b.size(); ++individual){
      
      W_b[individual][p] = exp(W_b[individual][p]); 
      W_b_new[individual][p] = exp(W_b_new[individual][p]); 
    } 
  }
  
  for(unsigned int p=0; p<W_b[0].size(); ++p){
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      split.min_sum += split.Y_s[individual][p] * log(W_s[individual][p] / (std::accumulate(W_s[individual].begin(),W_s[individual].end(),0.0)) ); // ~ R_old
      split.min_sum -= split.Y_s[individual][p] * log(W_s_new[individual][p] / (std::accumulate(W_s_new[individual].begin(),W_s_new[individual].end(),0.0)) ); // ~ R_new
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      split.min_sum += split.Y_b[individual][p] * log(W_b[individual][p] / (std::accumulate(W_b[individual].begin(),W_b[individual].end(),0.0)) ); // ~ R_old
      split.min_sum -= split.Y_b[individual][p] * log(W_b_new[individual][p] / (std::accumulate(W_b_new[individual].begin(),W_b_new[individual].end(),0.0)) ); // ~ R_new
    }
  } 
  
  if(std::isnan(split.min_sum)){
    split.min_sum = INF;
  } 
}

void ClassificationRPF::exponential_loss(rpf::Split& split){
  
  split.min_sum = 0;
  split.M_s =  std::vector<double>(value_size, 0);
  split.M_b =  std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);
  
  for(int p=0; p<value_size; ++p){
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      W_s_sum[p] += split.W_s[individual][p];
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      W_b_sum[p] += split.W_b[individual][p];
    }
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      sum_s[p] += ((split.Y_s[individual][p] + 1) / 2) * (split.W_s[individual][p] / W_s_sum[p]);
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      sum_b[p] += ((split.Y_b[individual][p] + 1) / 2) * (split.W_b[individual][p] / W_b_sum[p]);
    }
    
    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];
    
    sum_s[p] = std::max(delta, sum_s[p]);
    sum_b[p] = std::max(delta, sum_b[p]);
  }
  
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  
  double sum_sp = std::max(delta, split.M_sp);
  double sum_bp = std::max(delta, split.M_bp);
  
  for(unsigned int p=0; p<split.Y_s[0].size(); ++p){
    
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      split.min_sum += split.W_s[individual][p] * exp(-0.5 * split.Y_s[individual][p] * log(sum_s[p] / sum_sp));
    } 
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      split.min_sum += split.W_b[individual][p] * exp(-0.5 * split.Y_b[individual][p] * log(sum_b[p] / sum_bp));
    }
    
    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }
  
  // check if valid result
  for(const auto& s: W_s_sum) if(s == 0) split.min_sum = INF;
  for(const auto& s: W_b_sum) if(s == 0) split.min_sum = INF;
  if(std::isnan(split.min_sum)) split.min_sum = INF; 
}

void ClassificationRPF::exponential_loss_2(rpf::Split& split){
  
  split.min_sum = 0;
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);
  std::vector<double> sum_s2(value_size, 0);
  std::vector<double> sum_b2(value_size, 0);
  
  for(int p=0; p<value_size; ++p){
    
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      W_s_sum[p] += split.W_s[individual][p];
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      W_b_sum[p] += split.W_b[individual][p];
    }
    
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      sum_s[p] += ((split.Y_s[individual][p] + 1) / 2) * (split.W_s[individual][p] / W_s_sum[p]);
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      sum_b[p] += ((split.Y_b[individual][p] + 1) / 2) * (split.W_b[individual][p] / W_b_sum[p]);
    }
    
    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];
    
    sum_s2[p] = std::max(delta, 1 - sum_s[p]);
    sum_b2[p] = std::max(delta, 1 - sum_s[p]);
    
    sum_s[p] = std::max(delta, sum_s[p]);
    sum_b[p] = std::max(delta, sum_b[p]);
  }
  
  for(int p=0; p<value_size; ++p){
    
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      split.min_sum += split.W_s[individual][p] * exp(-0.5 * split.Y_s[individual][p] * log(sum_s[p] / sum_s2[p]));
    } 
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      split.min_sum += split.W_b[individual][p] * exp(-0.5 * split.Y_b[individual][p] * log(sum_b[p] / sum_b2[p]));
    }
    
    split.min_sum -= W_s_sum[p] + W_b_sum[p]; 
  }

  // check if valid result
  for(const auto& s: W_s_sum) if(s == 0) split.min_sum = INF;
  for(const auto& s: W_b_sum) if(s == 0) split.min_sum = INF;
  if(std::isnan(split.min_sum)) split.min_sum = INF;
}

void ClassificationRPF::exponential_loss_3(rpf::Split& split){
  
  split.min_sum = 0;
  split.M_s =  std::vector<double>(value_size, 0);
  split.M_b =  std::vector<double>(value_size, 0);
  std::vector<double> W_s_sum(value_size, 0);
  std::vector<double> W_b_sum(value_size, 0);
  std::vector<double> sum_s(value_size, 0);
  std::vector<double> sum_b(value_size, 0);
  
  for(int p=0; p<value_size; ++p){
    for(unsigned int individual=0; individual<split.W_s.size(); ++individual){
      W_s_sum[p] += split.W_s[individual][p];
    }
    for(unsigned int individual=0; individual<split.W_b.size(); ++individual){
      W_b_sum[p] += split.W_b[individual][p];
    }
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      sum_s[p] += ((split.Y_s[individual][p] + 1) / 2) * (split.W_s[individual][p] / W_s_sum[p]);
    }
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      sum_b[p] += ((split.Y_b[individual][p] + 1) / 2) * (split.W_b[individual][p] / W_b_sum[p]);
    }
    
    split.M_s[p] = sum_s[p];
    split.M_b[p] = sum_b[p];
    sum_s[p] = std::max(delta, sum_s[p]);
    sum_b[p] = std::max(delta, sum_b[p]);
    sum_s[p] = log(sum_s[p]);
    sum_b[p] = log(sum_b[p]);
  }
  
  split.M_sp = 1 - std::accumulate(split.M_s.begin(), split.M_s.end(), 0.0);
  split.M_bp = 1 - std::accumulate(split.M_b.begin(), split.M_b.end(), 0.0);
  
  double sum_sp = std::max(delta, split.M_sp);
  double sum_bp = std::max(delta, split.M_bp);
  
  sum_sp = log(sum_sp);
  sum_bp = log(sum_bp);

  sum_sp += std::accumulate(sum_s.begin(), sum_s.end(), 0.0);
  sum_bp += std::accumulate(sum_b.begin(), sum_b.end(), 0.0);
  
  sum_sp = sum_sp/(sum_s.size()+1);
  sum_bp = sum_bp/(sum_b.size()+1);
  
  for(unsigned int p=0; p<split.Y_s[0].size(); ++p){
    
    for(unsigned int individual=0; individual<split.Y_s.size(); ++individual){
      split.min_sum += split.W_s[individual][p] * exp(-0.5 * split.Y_s[individual][p] * (sum_s[p] - sum_sp));
    } 
    for(unsigned int individual=0; individual<split.Y_b.size(); ++individual){
      split.min_sum += split.W_b[individual][p] * exp(-0.5 * split.Y_b[individual][p] * (sum_b[p] - sum_bp));
    }
    
    split.min_sum -= W_s_sum[p] + W_b_sum[p];
  }
  
  // check if valid result
  for(const auto& s: W_s_sum) if(s == 0) split.min_sum = INF;
  for(const auto& s: W_b_sum) if(s == 0) split.min_sum = INF;
  if(std::isnan(split.min_sum)) split.min_sum = INF; 
}

// constructor with parameters split_try, t_try, purify_forest, deterministic, parallelize
ClassificationRPF::ClassificationRPF(const NumericMatrix& samples_Y, const NumericMatrix& samples_X,
                                     const String loss, const NumericVector parameters)
   : RandomPlantedForest{}{

  // Ensure correct Rcpp RNG state
  Rcpp::RNGScope scope; 

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
   }else if(loss == "logit_2"){
     this->loss = LossType::logit_2;
     this->calcLoss = &ClassificationRPF::logit_loss_2;
   }else if(loss == "logit_3"){
     this->loss = LossType::logit_3;
     this->calcLoss = &ClassificationRPF::logit_loss_3;
   }else if(loss == "exponential"){
     this->loss = LossType::exponential;
     this->calcLoss = &ClassificationRPF::exponential_loss;
   }else if(loss == "exponential_2"){
     this->loss = LossType::exponential_2;
     this->calcLoss = &ClassificationRPF::exponential_loss_2;
   }else if(loss == "exponential_3"){
     this->loss = LossType::exponential_3;
     this->calcLoss = &ClassificationRPF::exponential_loss_3;
   }else{
     Rcout << "Unkown loss function, set to default (L2)." << std::endl;
     this->loss = LossType::L2;
     this->calcLoss = &ClassificationRPF::L2_loss;
   }
  if(pars.size() != 11){
    Rcout << "Wrong number of parameters - set to default." << std::endl;
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
rpf::Split ClassificationRPF::calcOptimalSplit(const std::vector<std::vector<double>>& Y, const std::vector<std::vector<double>>& X,
                                                 std::multimap<int, std::shared_ptr<DecisionTree>>& possible_splits, TreeFamily& curr_family, std::vector<std::vector<double>>& weights){

  rpf::Split curr_split, min_split;
  std::set<int> tree_dims;
  int k;
  bool splitable;
  size_t n = 0;
  double leaf_size, sample_point;
  
  // sample possible splits
  int n_candidates = ceil(t_try*possible_splits.size()); // number of candidates that will be considered
  std::vector<int> split_candidates(possible_splits.size());
  std::iota(split_candidates.begin(), split_candidates.end(), 0); // consecutive indices of possible candidates
  
  if(!deterministic){
    std::random_shuffle(split_candidates.begin(), split_candidates.end(), randWrapper); // shuffle for random order
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
        std::vector<double> unique_samples(curr_individuals.size());
        for(unsigned int i=0; i<curr_individuals.size(); ++i){
          unique_samples[i] = X[curr_individuals[i]][k];
        }
        std::sort(unique_samples.begin(), unique_samples.end());
        unique_samples.erase(std::unique(unique_samples.begin(), unique_samples.end()), unique_samples.end());
        
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
          std::uniform_int_distribution<> distrib(leaf_size, unique_samples.size() - leaf_size);
          if(deterministic){
            std::advance(sample_pos, t);
          }else{
            int offset = R::runif(leaf_size, unique_samples.size() - leaf_size );
            std::advance(sample_pos, offset); // consider only sample points with offset
          }
          sample_point = *sample_pos;
          
          // clear current split
          {
            curr_split.I_s.clear();
            curr_split.I_b.clear();
            curr_split.Y_s.clear();
            curr_split.Y_b.clear();
            curr_split.W_s.clear();
            curr_split.W_b.clear();
            curr_split.M_s =  std::vector<double>(value_size, 0);
            curr_split.M_b =  std::vector<double>(value_size, 0);
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
            min_split = curr_split;
            min_split.tree_index = curr_tree;
            min_split.leaf_index =  &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point = sample_point;
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
    Rcout << "Initial Possible Splits: ";
    for(auto split: possible_splits){
      Rcout << split.first << "-";
      for(auto dim: split.second->split_dims) Rcout << dim << ",";
      Rcout << "; ";
    }
    Rcout << std::endl;
  }
  
  // sample data points with replacement
  int sample_index;
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;
  
  // deterministic
  if(deterministic){
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  }else{
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);
    
    for(size_t i=0; i<sample_size; ++i){
      sample_index = R::runif(0, sample_size - 1);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }
  
  // initialize weights
  // function pointer
  std::vector<std::vector<double>> weights; 
  switch(this->loss){
  case LossType::logit: case LossType::logit_2: case LossType::logit_3:
    weights = std::vector<std::vector<double>>(sample_size);
    for(auto& W: weights) W = std::vector<double>(value_size, 0);
    break;
  case LossType::exponential: case LossType::exponential_2: case LossType::exponential_3:
    weights = std::vector<std::vector<double>>(sample_size);
    for(auto& W: weights) W = std::vector<double>(value_size, 1);
    break;
  default:
    weights = std::vector<std::vector<double>>(sample_size);
    for(auto& W: weights) W = std::vector<double>(value_size, 0);
  }

  // modify existing or add new trees through splitting
  rpf::Split curr_split;
  for(size_t split_count=0; split_count<n_splits; ++split_count){
    
    // find optimal split
    curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family, weights);
    
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
          Rcout << "Updated Possible Splits: " << std::endl;
          for(auto split: possible_splits){
            Rcout << split.first << "-";
            for(auto dim: split.second->split_dims) Rcout << dim << ",";
            Rcout << "; ";
          }
          Rcout << std::endl;
        }
      }
          
      // update values of individuals of split interval
      std::vector<double> update_s = curr_split.M_s, update_b = curr_split.M_b;  
      switch(this->loss){
      case LossType::L1: case LossType::L2: case LossType::median: {
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

        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;
        
        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });
        
        double M_sp = std::max(epsilon, curr_split.M_sp);
        double M_bp = std::max(epsilon, curr_split.M_bp);
        
        std::vector<double> W_s_mean = calcMean(curr_split.W_s);
        std::vector<double> W_b_mean = calcMean(curr_split.W_b);
        
        for(unsigned int p=0; p<M_s.size(); ++p){
          update_s[p] = log(M_s[p] / M_sp) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_bp) - W_b_mean[p];
        }
        
        for(int individual: curr_split.leaf_index->individuals){
          if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
            weights[individual] += update_s;
          }else{
            weights[individual] += update_b;
          }
        }
        
        break;
      }
      case LossType::logit_2: {
        
        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;
        
        std::vector<double> M_s2 = curr_split.M_s;
        std::vector<double> M_b2 = curr_split.M_b;
        
        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });
        
        std::for_each(M_s2.begin(), M_s2.end(), [this](double &M) { M = std::max(epsilon, 1-M); });
        std::for_each(M_b2.begin(), M_b2.end(), [this](double &M) { M = std::max(epsilon, 1-M); });
        
        std::vector<double> W_s_mean = calcMean(curr_split.W_s);
        std::vector<double> W_b_mean = calcMean(curr_split.W_b);
        
        for(int p=0; p<value_size; ++p){
          update_s[p] = log(M_s[p] / M_s2[p]) - W_s_mean[p];
          update_b[p] = log(M_b[p] / M_b2[p]) - W_b_mean[p];
        }
        
        for(int individual: curr_split.leaf_index->individuals){
          if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
            weights[individual] += update_s;
          }else{
            weights[individual] += update_b;
          }
        }
        
        break;
      }
      case LossType::logit_3: {
        
        std::vector<double> M_s = curr_split.M_s;
        std::vector<double> M_b = curr_split.M_b;
        
        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = std::max(epsilon, M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = std::max(epsilon, M); });
        
        std::for_each(M_s.begin(), M_s.end(), [this](double &M) { M = log(M); });
        std::for_each(M_b.begin(), M_b.end(), [this](double &M) { M = log(M); });

        double M_sp = std::max(epsilon, curr_split.M_sp);
        double M_bp = std::max(epsilon, curr_split.M_bp);
        
        M_sp = log(M_sp);
        M_bp = log(M_bp);

        double sum_s = (std::accumulate(M_s.begin(), M_s.end(),0.0)+M_sp)/(M_s.size()+1);
        double sum_b = (std::accumulate(M_b.begin(), M_b.end(),0.0)+M_bp)/(M_b.size()+1);
        
        std::vector<double> W_s_mean = calcMean(curr_split.W_s);
        std::vector<double> W_b_mean = calcMean(curr_split.W_b);
        
        for(unsigned int p=0; p<M_s.size(); ++p){
          update_s[p] = M_s[p] - sum_s - W_s_mean[p];
          update_b[p] = M_b[p] - sum_b - W_b_mean[p];
        }
        
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
        
        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;
        
        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = std::max(epsilon, S); });
        
        double sum_sp = std::max(epsilon, curr_split.M_sp);
        double sum_bp = std::max(epsilon, curr_split.M_bp);  
        
        for(unsigned int p = 0; p<sum_s.size(); ++p){
          update_s[p] = log(sum_s[p] / sum_sp);
          update_b[p] = log(sum_b[p] / sum_bp);
        }
        
        for(int individual: curr_split.leaf_index->individuals){
          for(unsigned int p = 0; p<update_s.size(); ++p){
            if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }else{
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
            }
          }
        }
        
        break;
      }
      case LossType::exponential_2: {
        
         std::vector<double> sum_s = curr_split.M_s;
         std::vector<double> sum_b = curr_split.M_b;
         std::vector<double> sum_s2 = curr_split.M_s;
         std::vector<double> sum_b2 = curr_split.M_b;
         
         std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = std::max(epsilon, S); });
         std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = std::max(epsilon, S); });
         
         std::for_each(sum_s2.begin(), sum_s2.end(), [this](double &S) { S = std::max(epsilon, 1 - S); });
         std::for_each(sum_b2.begin(), sum_b2.end(), [this](double &S) { S = std::max(epsilon, 1 - S); }); 
         
         for(int p = 0; p<value_size; ++p){
           update_s[p] = log(sum_s[p] / sum_s2[p]);
           update_b[p] = log(sum_b[p] / sum_b2[p]);
         }

         for(int individual: curr_split.leaf_index->individuals){
           for(int p = 0; p<value_size; ++p){
             if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
             }else{
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
             }
           }
         }
                  
        break;
      }
      case LossType::exponential_3: {
        
        std::vector<double> sum_s = curr_split.M_s;
        std::vector<double> sum_b = curr_split.M_b;
        
        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = std::max(epsilon, S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = std::max(epsilon, S); });
        
        std::for_each(sum_s.begin(), sum_s.end(), [this](double &S) { S = log(S); });
        std::for_each(sum_b.begin(), sum_b.end(), [this](double &S) { S = log(S); });

        double sum_sp = std::max(epsilon, curr_split.M_sp);
        double sum_bp = std::max(epsilon, curr_split.M_bp);  
        
        sum_sp = log(sum_sp);
        sum_bp = log(sum_bp);

        sum_sp += std::accumulate(sum_s.begin(), sum_s.end(),0.0);
        sum_bp += std::accumulate(sum_b.begin(), sum_b.end(),0.0);
        
        sum_sp = sum_sp/(sum_s.size()+1);
        sum_bp = sum_bp/(sum_b.size()+1);
        
        for(unsigned int p = 0; p<sum_s.size(); ++p){
          update_s[p] = sum_s[p] - sum_sp;
          update_b[p] = sum_b[p] - sum_bp;
        }
        
        for(int individual: curr_split.leaf_index->individuals){
          for(unsigned int p = 0; p<update_s.size(); ++p){
            if(samples_X[individual][curr_split.split_coordinate-1] < curr_split.split_point){
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_s[p]);
            }else{
              weights[individual][p] *= exp(-0.5 * samples_Y[individual][p] * update_b[p]);
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
    initial_leaf.value = std::vector<double>(value_size, 0);
    initial_leaf.individuals = initial_individuals;
    initial_leaf.intervals = initial_intervals;
  }
  std::vector<Leaf> initial_leaves{initial_leaf}; // vector with initial leaf
  
  // initialize tree families
  this->tree_families = std::vector<TreeFamily>(n_trees); 
  
  // iterate over families of trees and modify
  if(parallelize){
    int n_threads = std::thread::hardware_concurrency() - 1;
    for(int n = 0; n<n_trees; n+=n_threads){
      if(n >= (n_trees - n_threads)) n_threads = n_trees % n_threads;
      std::vector<std::thread> threads(n_threads);
      for(size_t t=0; t<n_threads; ++t){
        threads[t] = std::thread(&ClassificationRPF::create_tree_family, this, std::ref(initial_leaves), n + t);
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

/*  retrospectively change parameters of existing class object, 
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data', 
 note that changing cv does not trigger cross validation */
void ClassificationRPF::set_parameters(StringVector keys, NumericVector values){
  if(keys.size() != values.size()) {
    Rcout << "Size of input vectors is not the same. " << std::endl;
    return;
  }
  
  for(unsigned int i=0; i<keys.size(); ++i){
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
        Rcout << "Unkown loss function." << std::endl;
      }
    }else if(keys[i] == "delta"){
      this->delta = values[i];
    }else if(keys[i] == "epsilon"){
      this->epsilon = values[i];
    }else{
      Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}


// ----------------- Rcpp include  -----------------

RCPP_MODULE(mod_rpf) {

    class_<RandomPlantedForest>("RandomPlantedForest")
      .constructor<const NumericMatrix, const NumericMatrix, const NumericVector>()
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
      .constructor<const NumericMatrix, const NumericMatrix, const String, const NumericVector>()
      .method("set_parameters", &ClassificationRPF::set_parameters)
    ;
}
