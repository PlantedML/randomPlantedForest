#include "rpf.hpp"


bool RandomPlantedForest::is_purified()
{
  return purified;
}

void RandomPlantedForest::L2_loss(Split &split)
{

  // new meanq
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();

  split.min_sum = 0;
  for (size_t p = 0; p < value_size; ++p)
  {
    split.min_sum += -2 * split.M_s[p] * split.sum_s[p] + split.I_s.size() * pow(split.M_s[p], 2);
    split.min_sum += -2 * split.M_b[p] * split.sum_b[p] + split.I_b.size() * pow(split.M_b[p], 2);
  }
}

// constructor
RandomPlantedForest::RandomPlantedForest(const NumericMatrix &samples_Y, const NumericMatrix &samples_X,
                                         const NumericVector parameters)
{

  // Ensure correct Rcpp RNG state
  Rcpp::RNGScope scope;

  // initialize class members
  std::vector<double> pars = to_std_vec(parameters);
  if (pars.size() != 9)
  {
    Rcout << "Wrong number of parameters - set to default." << std::endl;
    this->max_interaction = 1;
    this->n_trees = 50;
    this->n_splits = 30;
    this->split_try = 10;
    this->t_try = 0.4;
    this->purify_forest = 0;
    this->deterministic = 0;
    this->nthreads = 1;
    this->cross_validate = 0;
  }
  else
  {
    this->max_interaction = pars[0];
    this->n_trees = pars[1];
    this->n_splits = pars[2];
    this->split_try = pars[3];
    this->t_try = pars[4];
    this->purify_forest = pars[5];
    this->deterministic = pars[6];
    this->nthreads = pars[7];
    this->cross_validate = pars[8];
  }

  // set data and data related members
  this->set_data(samples_Y, samples_X);
}

// determine optimal split
Split RandomPlantedForest::calcOptimalSplit(const std::vector<std::vector<double>> &Y, const std::vector<std::vector<double>> &X,
                                            std::multimap<int, std::shared_ptr<DecisionTree>> &possible_splits, TreeFamily &curr_family)
{

  Split curr_split, min_split;
  curr_split.Y = &Y;
  std::set<int> tree_dims;
  std::vector<double> unique_samples;
  int k;
  unsigned int n = 0;
  double leaf_size, sample_point;

  // sample possible splits
  unsigned int n_candidates = ceil(t_try * possible_splits.size()); // number of candidates that will be considered
  std::vector<int> split_candidates(possible_splits.size());
  std::iota(split_candidates.begin(), split_candidates.end(), 0); // consecutive indices of possible candidates

  if (!deterministic)
  {
    shuffle_vector(split_candidates.begin(), split_candidates.end()); // shuffle for random order
  }

  // consider a fraction of possible splits
  while (n < n_candidates)
  {

    if (possible_splits.empty())
      break;
    if (split_candidates[n] >= 0 && (size_t)split_candidates[n] >= possible_splits.size())
      continue;

    auto candidate = possible_splits.begin();
    std::advance(candidate, split_candidates[n]); // get random split candidate without replacement
    k = candidate->first - 1;                     // split dim of current candidate, converted to index starting at 0
    leaf_size = n_leaves[k];

    // Test if splitting in the current tree w.r.t. the coordinate "k" is an element of candidate tree
    tree_dims = candidate->second->split_dims;
    tree_dims.erase(k + 1);
    tree_dims.erase(0);

    std::vector<std::shared_ptr<DecisionTree>> curr_trees;
    if (tree_dims.size() == 0)
      curr_trees.push_back(curr_family[std::set<int>{0}]);
    if (curr_family.find(tree_dims) != curr_family.end())
      curr_trees.push_back(curr_family[tree_dims]);
    if (curr_family.find(candidate->second->split_dims) != curr_family.end())
      curr_trees.push_back(curr_family[candidate->second->split_dims]);

    // go through all trees in current family
    for (auto &curr_tree : curr_trees)
    {

      // skip if tree has no leaves
      if (curr_tree->leaves.size() == 0)
        continue;

      // go through all leaves of current tree
      for (auto &leaf : curr_tree->leaves)
      {

        std::vector<double> tot_sum(value_size, 0);

        // extract sample points according to individuals from X and Y
        unique_samples = std::vector<double>(leaf.individuals.size());
        for (unsigned int i = 0; i < leaf.individuals.size(); ++i)
        {
          unique_samples[i] = X[leaf.individuals[i]][k];
        }
        std::sort(unique_samples.begin(), unique_samples.end());
        unique_samples.erase(std::unique(unique_samples.begin(), unique_samples.end()), unique_samples.end());

        // check if number of sample points is within limit
        if (unique_samples.size() < 2 * leaf_size)
          continue;

        // consider split_try-number of samples
        std::vector<int> samples;
        if (deterministic)
        { // sequential samples if deterministic
          samples = std::vector<int>(std::min((int)unique_samples.size() - 1, 9));
          std::iota(samples.begin(), samples.end(), 1);
        }
        else
        { // randomly picked samples otherwise
          samples = std::vector<int>(split_try);
          for (size_t i = 0; i < samples.size(); ++i)
            samples[i] = R::runif(leaf_size, unique_samples.size() - leaf_size);
          std::sort(samples.begin(), samples.end());
        }

        // go through samples
        for (size_t sample_pos = 0; sample_pos < samples.size(); ++sample_pos)
        {

          // get samplepoint
          sample_point = unique_samples[samples[sample_pos]];

          // clear current split
          {
            curr_split.I_s.clear();
            curr_split.I_b.clear();
            curr_split.I_s.reserve(leaf.individuals.size());
            curr_split.I_b.reserve(leaf.individuals.size());
            curr_split.M_s = std::vector<double>(value_size, 0);
            curr_split.M_b = std::vector<double>(value_size, 0);
          }

          // get samples greater/smaller than samplepoint
          if (sample_pos == 0)
          {
            curr_split.sum_s = std::vector<double>(value_size, 0);
            curr_split.sum_b = std::vector<double>(value_size, 0);

            for (int individual : leaf.individuals)
            {
              if (X[individual][k] < sample_point)
              {
                curr_split.I_s.push_back(individual);
                curr_split.sum_s += Y[individual];
              }
              else
              {
                curr_split.I_b.push_back(individual);
                curr_split.sum_b += Y[individual];
              }
            }

            tot_sum = curr_split.sum_s + curr_split.sum_b;
          }
          else
          {

            for (int individual : leaf.individuals)
            {
              if (X[individual][k] < sample_point)
              {
                if (X[individual][k] >= unique_samples[samples[sample_pos - 1]])
                {
                  curr_split.sum_s += Y[individual];
                }
                curr_split.I_s.push_back(individual);
              }
              else
              {
                curr_split.I_b.push_back(individual);
              }
            }

            curr_split.sum_b = tot_sum - curr_split.sum_s;
          }

          // accumulate squared mean and get mean
          L2_loss(curr_split);

          // update split if squared sum is smaller
          if (curr_split.min_sum < min_split.min_sum)
          {
            min_split = curr_split;
            min_split.tree_index = curr_tree;
            min_split.leaf_index = &leaf;
            min_split.split_coordinate = k + 1;
            min_split.split_point = sample_point;
          }
        }
      }
    }

    ++n;
  }

  return min_split;
}

void RandomPlantedForest::set_data(const NumericMatrix &samples_Y, const NumericMatrix &samples_X)
{

  this->Y = to_std_vec(samples_Y);
  this->X = to_std_vec(samples_X);

  // Check for correct input
  if (Y.size() == 0)
    throw std::invalid_argument("Y empty - no data provided.");
  if (X.size() == 0)
    throw std::invalid_argument("X empty - no data provided.");
  this->feature_size = X[0].size();
  this->value_size = Y[0].size(); // multiclass
  for (const auto &vec : X)
  {
    if (vec.size() != (size_t)feature_size)
      throw std::invalid_argument("Feature dimensions of X not uniform.");
  }
  if (Y.size() != X.size())
    throw std::invalid_argument("X and Y are not of the same length!");

  this->n_leaves = std::vector<int>(feature_size, 1);
  this->sample_size = X.size();
  this->upper_bounds = std::vector<double>(feature_size);
  this->lower_bounds = std::vector<double>(feature_size);

  // get upper/lower bounds
  double minVal, maxVal, currVal;
  for (int i = 0; i < feature_size; ++i)
  {
    minVal = maxVal = X[0][i];
    for (size_t j = 0; j < sample_size; ++j)
    {
      currVal = X[j][i];
      if (currVal < minVal)
        minVal = currVal;
      if (currVal > maxVal)
        maxVal = currVal;
    }
    this->upper_bounds[i] = maxVal + 2 * eps; // to consider samples at max value
    this->lower_bounds[i] = minVal;
  }

  this->fit();

  if (cross_validate)
  {
    this->cross_validation();
  }
}

void RandomPlantedForest::create_tree_family(std::vector<Leaf> initial_leaves, size_t n)
{

  TreeFamily curr_family;
  curr_family.insert(std::make_pair(std::set<int>{0}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{0}, initial_leaves)))); // save tree with one leaf in the beginning
  // store possible splits in map with splitting variable as key and pointer to resulting tree
  std::multimap<int, std::shared_ptr<DecisionTree>> possible_splits;
  for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
  {
    // add pointer to resulting tree with split dimension as key
    curr_family.insert(std::make_pair(std::set<int>{feature_dim}, std::make_shared<DecisionTree>(DecisionTree(std::set<int>{feature_dim}))));
    possible_splits.insert(std::make_pair(feature_dim, curr_family[std::set<int>{feature_dim}]));
  }

  // sample data points with replacement
  int sample_index;
  std::vector<std::vector<double>> samples_X;
  std::vector<std::vector<double>> samples_Y;

  // deterministic
  if (deterministic)
  {
    samples_X = X;
    samples_Y = Y;
    this->t_try = 1;
  }
  else
  {
    samples_X = std::vector<std::vector<double>>(sample_size);
    samples_Y = std::vector<std::vector<double>>(sample_size);

    for (size_t i = 0; i < sample_size; ++i)
    {

      sample_index = R::runif(0, sample_size - 1);
      samples_Y[i] = Y[sample_index];
      samples_X[i] = X[sample_index];
    }
  }

  // modify existing or add new trees through splitting
  Split curr_split;
  for (int split_count = 0; split_count < n_splits; ++split_count)
  {

    // find optimal split
    curr_split = calcOptimalSplit(samples_Y, samples_X, possible_splits, curr_family);

    // continue only if we get a significant result
    if (!std::isinf(curr_split.min_sum))
    {

      // update possible splits
      for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
      { // consider all possible dimensions

        // create union of split coord, feature dim and dimensions of old tree
        std::set<int> curr_dims = curr_split.tree_index->split_dims;
        curr_dims.insert(curr_split.split_coordinate);
        curr_dims.insert(feature_dim);
        curr_dims.erase(0);

        // skip if possible_split already exists
        if (possibleExists(feature_dim, possible_splits, curr_dims))
          continue;

        // do not exceed maximum level of interaction
        if (max_interaction >= 0 && curr_dims.size() > (size_t)max_interaction)
          continue;

        // check if resulting tree already exists in family
        std::shared_ptr<DecisionTree> found_tree = treeExists(curr_dims, curr_family);

        // update possible_splits if not already existing
        if (found_tree)
        { // if yes add pointer
          possible_splits.insert(std::make_pair(feature_dim, found_tree));
        }
        else
        { // if not create new tree
          curr_family.insert(std::make_pair(curr_dims, std::make_shared<DecisionTree>(DecisionTree(curr_dims))));
          possible_splits.insert(std::make_pair(feature_dim, curr_family[curr_dims]));
        }
      }

      // update values of individuals of split interval with mean
      for (int individual : curr_split.leaf_index->individuals)
      { // todo: loop directly over I_s I_b
        if (samples_X[individual][curr_split.split_coordinate - 1] < curr_split.split_point)
        {
          samples_Y[individual] -= curr_split.M_s;
        }
        else
        {
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
        leaf_s.intervals[curr_split.split_coordinate - 1].second = curr_split.split_point;
        // interval of leaf with bigger individuals has new lower bound in splitting dimension
        leaf_b.intervals[curr_split.split_coordinate - 1].first = curr_split.split_point;
      }

      // construct split_dims of resulting tree when splitting in split_coordinate
      std::set<int> resulting_dims = curr_split.tree_index->split_dims;
      resulting_dims.insert(curr_split.split_coordinate);
      resulting_dims.erase(0);

      // check if resulting tree already exists in family
      std::shared_ptr<DecisionTree> found_tree = treeExists(resulting_dims, curr_family);

      // determine which tree is modified
      if (curr_split.tree_index->split_dims.count(curr_split.split_coordinate))
      { // if split variable is already in tree to be split
        // change values
        {
          leaf_s.value += curr_split.leaf_index->value;
          leaf_b.value += curr_split.leaf_index->value;
        }
        *curr_split.leaf_index = leaf_b;                 // replace old interval
        curr_split.tree_index->leaves.push_back(leaf_s); // add new leaf
      }
      else
      {                                       // otherwise
        found_tree->leaves.push_back(leaf_s); // append new leaves
        found_tree->leaves.push_back(leaf_b);
      }
    }
  }

  // remove empty trees & clear individuals of each tree
  auto keys = getKeys(curr_family);
  for (auto &key : keys)
  {
    if (curr_family[key]->leaves.size() == 0)
    {
      curr_family.erase(key);
      continue;
    }
    for (auto &leaf : curr_family[key]->leaves)
    {
      leaf.individuals.clear();
    }
  }

  tree_families[n] = curr_family;
}

// fit forest to new data
void RandomPlantedForest::fit()
{

  // setup initial set of individuals
  std::vector<int> initial_individuals(sample_size);
  std::iota(initial_individuals.begin(), initial_individuals.end(), 0);

  // initialize intervals with lower and upper bounds
  std::vector<Interval> initial_intervals(feature_size);
  for (int i = 0; i < feature_size; ++i)
    initial_intervals[i] = Interval{lower_bounds[i], upper_bounds[i]};

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

  // Loop over number of tree families and dispatch threads in batches
  // of nhreads at once
  if (nthreads > 1)
  {
    if (nthreads > std::thread::hardware_concurrency())
    {
      Rcout << "Requested " << nthreads << " threads but only " << std::thread::hardware_concurrency() << " available" << std::endl;
    }
    // Create local thread count to not overwrite nthreads,
    // would get reported wrongly by get_parameters()
    unsigned int current_threads = nthreads;
    for (int n = 0; n < n_trees; n += current_threads)
    {
      if (n >= (n_trees - current_threads + 1))
      {
        current_threads = n_trees % current_threads;
      }

      std::vector<std::thread> threads(current_threads);
      for (int t = 0; t < current_threads; ++t)
      {
        // Rcout << "Dispatching thread " << (n + t + 1) << "/" << n_trees << std::endl;
        threads[t] = std::thread(&RandomPlantedForest::create_tree_family, this, std::ref(initial_leaves), n + t);
      }
      for (auto &t : threads)
      {
        if (t.joinable())
          t.join();
      }
    }
  }
  else
  {
    for (int n = 0; n < n_trees; ++n)
    {
      create_tree_family(initial_leaves, n);
    }
  }

  // optionally purify tree
  if (purify_forest)
  {
    this->purify_3();
  }
  else
  {
    purified = false;
  }
}

void RandomPlantedForest::cross_validation(int n_sets, IntegerVector splits, NumericVector t_tries, IntegerVector split_tries)
{

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
std::vector<double> RandomPlantedForest::predict_single(const std::vector<double> &X, std::set<int> component_index)
{

  std::vector<double> total_res = std::vector<double>(value_size, 0);

  if (!purified)
  {
    // consider all components
    if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (auto &dim : tree.first)
            {
              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[std::max(0, dim - 1)] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              {
                valid = false;
              }
            }
            if (valid)
            {

              // Rcout << leaf.value[0] << "\n";
              total_res += leaf.value;
            }
          }
        }
      }
    }
    else
    { // choose components for prediction
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {

          // only consider trees with same dimensions as component_index
          if (tree.first != component_index)
            continue;

          std::vector<int> dims;
          for (auto dim : tree.first)
          {
            dims.push_back(dim);
          }

          for (auto &leaf : tree.second->leaves)
          {
            bool valid = true;
            for (unsigned int i = 0; i < dims.size(); ++i)
            {

              int dim = dims[i];

              if (!((leaf.intervals[std::max(0, dim - 1)].first <= X[i] || leaf.intervals[std::max(0, dim - 1)].first == lower_bounds[std::max(0, dim - 1)]) && (leaf.intervals[std::max(0, dim - 1)].second > X[i] || leaf.intervals[std::max(0, dim - 1)].second == upper_bounds[std::max(0, dim - 1)])))
              {
                valid = false;
              }
            }
            if (valid)
              total_res += leaf.value;
          }
        }
      }
    }
  }
  else
  {
    if (component_index == std::set<int>{-1})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);
          // add value of null tree
          if (tree.first == std::set<int>{0})
          {

            // Rcout << tree.first.size() ;
            leaf_index = std::vector<int>(tree.first.size(), 0);
            total_res += tree.second->GridLeaves.values[leaf_index];
          }
        }
      }
    }
    else if (component_index == std::set<int>{0})
    {
      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {
          std::vector<int> leaf_index(tree.first.size(), -1);

          // add value of null tree
          if (tree.first == std::set<int>{0})
          {

            // Rcout << tree.first.size() ;
            leaf_index = std::vector<int>(tree.first.size(), 0);
          }
          else
          {

            // go through limits of grid
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              // get dim at dim_index
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // transform into index
              }

              auto bounds = tree.second->GridLeaves.lim_list[dim];
              for (double bound : bounds)
              {

                // check if sample in leaf at dimension
                if (X[dim] < bound)
                  break; // changed

                // if no interval smaller, set to end of bounds, otherwise set to leaf index
                leaf_index[dim_index] = std::min(leaf_index[dim_index] + 1, (int)bounds.size() - 2);
              }
            }
          }

          // if interval of first leaf smaller smaller
          for (int &index : leaf_index)
            index = std::max(0, index);

          total_res += tree.second->GridLeaves.values[leaf_index];
        }
      }
    }
    else
    {

      for (auto &tree_family : this->tree_families)
      {
        for (auto &tree : tree_family)
        {

          // only consider trees with same dimensions as component_index
          if (tree.first != component_index)
            continue;

          std::vector<int> leaf_index(tree.first.size(), -1);

          // add value of null tree
          if (tree.first == std::set<int>{0})
          {
            leaf_index = std::vector<int>(tree.first.size(), 0);
          }
          else
          {

            // go through limits of grid
            for (size_t dim_index = 0; dim_index < tree.first.size(); ++dim_index)
            {
              // get dim at dim_index
              int dim = 0;
              {
                auto dim_pnt = tree.first.begin();
                std::advance(dim_pnt, dim_index);
                dim = *dim_pnt;
                --dim; // transform into index
              }

              auto bounds = tree.second->GridLeaves.lim_list[dim];
              for (double bound : bounds)
              {

                // check if sample in leaf at dimension
                if (X[dim_index] < bound)
                  break; // changed

                // if no interval smaller, set to end of bounds, otherwise set to leaf index
                leaf_index[dim_index] = std::min(leaf_index[dim_index] + 1, (int)bounds.size() - 2);
              }
            }
          }

          // if interval of first leaf smaller smaller
          for (int &index : leaf_index)
            index = std::max(0, index);

          total_res += tree.second->GridLeaves.values[leaf_index];
        }
      }
    }
  }

  return total_res / n_trees;
}

// predict multiple feature vectors
Rcpp::NumericMatrix RandomPlantedForest::predict_matrix(const NumericMatrix &X, const NumericVector components)
{
  std::vector<std::vector<double>> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions;

  // todo: sanity check for X
  if (feature_vec.empty())
    throw std::invalid_argument("Feature vector is empty.");
  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec[0].size() != (size_t)this->feature_size)
    throw std::invalid_argument("Feature vector has wrong dimension.");
  if (component_index != std::set<int>{0} && component_index != std::set<int>{-1} && component_index.size() != feature_vec[0].size())
    throw std::invalid_argument("The input X has the wrong dimension in order to calculate f_i(x)");

  for (auto &vec : feature_vec)
  {
    predictions.push_back(predict_single(vec, component_index));
  }

  return from_std_vec(predictions);
}

Rcpp::NumericMatrix RandomPlantedForest::predict_vector(const NumericVector &X, const NumericVector components)
{
  std::vector<double> feature_vec = to_std_vec(X);
  std::set<int> component_index = to_std_set(components);
  std::vector<std::vector<double>> predictions;
  Rcpp::NumericMatrix res;

  // todo: sanity check for X
  if (feature_vec.empty())
  {
    Rcout << "Feature vector is empty." << std::endl;
    return res;
  }

  if (component_index == std::set<int>{0} && this->feature_size >= 0 && feature_vec.size() != (size_t)this->feature_size)
  {
    Rcout << "Feature vector has wrong dimension." << std::endl;
    return res;
  }

  if (component_index == std::set<int>{0})
  {
    predictions.push_back(predict_single(feature_vec, component_index));
  }
  else
  {
    for (auto vec : feature_vec)
    {
      predictions.push_back(predict_single(std::vector<double>{vec}, component_index));
    }
  }

  res = from_std_vec(predictions);
  return res;
}

double RandomPlantedForest::MSE_vec(const NumericVector &Y_predicted, const NumericVector &Y_true)
{
  return sum(Rcpp::pow(Y_true - Y_predicted, 2)) / Y_true.size();
}

double RandomPlantedForest::MSE(const NumericMatrix &Y_predicted, const NumericMatrix &Y_true)
{
  // todo: multiclass
  double sum = 0;
  int Y_size = Y_predicted.size();

  for (int i = 0; i < Y_size; ++i)
  {
    sum += MSE_vec(Y_predicted(i, _), Y_true(i, _));
  }

  return sum / Y_size;
}

void RandomPlantedForest::purify_1()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (auto tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    while (curr_max >= 1)
    {

      // go through split dimensions of all trees
      auto keys = getKeys(curr_family);
      std::vector<std::set<int>>::reverse_iterator key = keys.rbegin();
      while (key != keys.rend())
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (curr_max == 1)
              {
                tree = curr_family[std::set<int>{0}];
              }
              else
              {
                if (!tree)
                {
                  curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                  tree = curr_family[tree_dims];
                }
              }

              // go through leaves of current tree
              int n_leaves = curr_tree->leaves.size();
              for (int l = 0; l < n_leaves; ++l)
              {
                auto &curr_leaf = curr_tree->leaves[l];

                double multiplier = (curr_leaf.intervals[feature_dim - 1].second - curr_leaf.intervals[feature_dim - 1].first) / (upper_bounds[feature_dim - 1] - lower_bounds[feature_dim - 1]);

                // new leaf including intervals and value
                Leaf new_leaf = curr_leaf; // initialize intervals with first leaf
                new_leaf.intervals[feature_dim - 1].first = lower_bounds[feature_dim - 1];
                new_leaf.intervals[feature_dim - 1].second = upper_bounds[feature_dim - 1];
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = -curr_leaf.value[i] * multiplier; // update value of new leaf

                // append new leaf
                if (!leafExists(new_leaf.intervals, curr_tree))
                  curr_tree->leaves.push_back(new_leaf);
                for (size_t i = 0; i < value_size; ++i)
                  new_leaf.value[i] = curr_leaf.value[i] * multiplier; // update value of new leaf
                if (!leafExists(new_leaf.intervals, tree))
                  tree->leaves.push_back(new_leaf);
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

  purified = true;
}

void RandomPlantedForest::purify_2()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);

    // go through all variables of the component
    for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
    {
      std::vector<double> bounds;

      // go through trees of family
      for (const auto &curr_tree : curr_family)
      {

        // consider only relevant trees that have current dimension as variable
        if (!curr_tree.first.count(curr_dim))
          continue;

        // go through leaves of tree
        for (const auto &curr_leaf : curr_tree.second->leaves)
        {
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
        }
      }
      std::sort(bounds.begin(), bounds.end());
      bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
      lim_list[curr_dim - 1] = bounds;
    }

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids(curr_family.size() - 1);
    std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
    std::vector<std::set<int>> variables(curr_family.size() - 1);

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {

      if (curr_tree.first == std::set<int>{0})
        continue; // ignore null tree

      // fill space with dimensions
      std::vector<int> dimensions;
      for (const auto &dim : curr_tree.first)
      {
        dimensions.push_back(lim_list[dim - 1].size() - 1); // size - 1 ?
      }

      // setup grid for leaf indices
      auto grid = grid::NDGrid(dimensions);

      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
      values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
      variables[tree_index] = curr_tree.first;

      // fill grid points with individuals and values
      while (!grid.nextPoint())
      {

        std::vector<int> gridPoint = grid.getPoint();

        bool in_leaf = true;

        // go through sample points to sum up individuals
        for (const auto &feature_vec : X)
        {
          int dim_index = 0;
          in_leaf = true;
          for (const auto &dim : curr_tree.first)
          {
            double val = feature_vec[dim - 1];
            if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // consider individuals only if all in
          if (in_leaf)
            individuals[tree_index][gridPoint] += 1;
        }

        // go through leaves of tree to sum up values
        for (const auto &leaf : curr_tree.second->get_leaves())
        {

          in_leaf = true;
          int dim_index = 0;
          for (const auto &dim : curr_tree.first)
          {
            // consider values only if all in
            if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // sum up values
          if (in_leaf)
            values[tree_index][gridPoint] += leaf.value; // todo: multiclass
        }
      }

      ++tree_index;
    }

    // ------------- create new trees -------------

    // insert null tree
    grids.insert(grids.begin(), grid::NDGrid());
    values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
    variables.insert(variables.begin(), std::set<int>{0});

    // recap maximum number of dimensions of current family
    unsigned int curr_max = 0;
    for (const auto &tree : curr_family)
    {
      if (tree.first.size() > curr_max)
        curr_max = tree.first.size();
    }

    auto keys = getKeys(curr_family);
    while (curr_max > 1)
    {

      // go through split dimensions of all trees
      for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
      {

        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;

        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {

          // go through feature dims
          int dim_index = 0;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {

            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {

              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree

              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (!tree)
              {

                // get index of old and new tree
                auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
                curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));

                // remove matrix dimension of respective variable
                std::vector<int> matrix_dimensions = values[old_tree_index].dims;
                matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);

                // initialize data for new tree
                auto grid = grid::NDGrid(matrix_dimensions);
                grids.insert(grids.begin() + tree_index, grid);
                values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(0, value_size)));
                individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
                variables.insert(variables.begin() + tree_index, tree_dims);

                // fill individuals of new trees
                while (!grid.nextPoint())
                {

                  std::vector<int> gridPoint = grid.getPoint();
                  bool in_leaf = true;

                  // go through sample points to sum up individuals
                  for (const auto &feature_vec : X)
                  {
                    int dim_index = 0;
                    in_leaf = true;
                    for (const auto &dim : tree_dims)
                    {
                      double val = feature_vec[dim - 1];
                      if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
                        in_leaf = false;
                      ++dim_index;
                    }

                    // consider individuals only if all in
                    if (in_leaf)
                      individuals[tree_index][gridPoint] += 1;
                  }
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
    for (TreeFamily::reverse_iterator curr_tree = curr_family.rbegin(); curr_tree != curr_family.rend(); ++curr_tree)
    {
      iter = 0;
      std::set<int> curr_dims = curr_tree->second->get_split_dims();

      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;

      // repeat until tolerance small enough and (?) maximum number of iterations reached
      while ((tol[curr_tree_index] > 0.00000000001) && (iter < 100))
      {

        // go through feature dims
        int curr_dim_index = 0;
        for (const auto &feature_dim : curr_dims)
        {

          // get tree that has same variables as curr_tree minus j-variable
          std::set<int> tree_dims = curr_dims;
          tree_dims.erase(tree_dims.find(feature_dim));
          int tree_index = 0; // if tree not exist, set to null tree
          if (curr_family.find(tree_dims) != curr_family.end())
            tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims)) - 1;

          // update values
          if (grids[curr_tree_index].dimensions.size() == 1)
          { // one dimensional case

            int sum_ind = 0;
            std::vector<double> avg(value_size, 0);

            // get sum of individuals
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              sum_ind += individuals[curr_tree_index][tmp];
            }
            if (sum_ind == 0)
              continue;

            // calc avg
            for (int i = 0; i < individuals[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              avg += (individuals[curr_tree_index][tmp] * values[curr_tree_index][tmp]) / sum_ind;
            }

            // update values of one dimensional and null tree
            for (int i = 0; i < values[curr_tree_index].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              values[curr_tree_index][tmp] -= avg;
            }
            std::vector<int> tmp{0};
            values[tree_index][tmp] += avg;
          }
          else
          { // higher dimensional case

            // setup new grid without dimension j
            std::vector<int> new_dimensions = grids[curr_tree_index].dimensions;
            int j_dim = new_dimensions[curr_dim_index];
            new_dimensions.erase(new_dimensions.begin() + curr_dim_index);
            grid::NDGrid grid = grid::NDGrid(new_dimensions);

            // go through values without dimension j
            while (!grid.nextPoint())
            {
              auto gridPoint = grid.getPoint();
              gridPoint.push_back(0);

              int sum_ind = 0;
              std::vector<double> avg(value_size, 0);

              // go through slice to sum up individuals
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // get sum of individuals
                sum_ind += individuals[curr_tree_index][gridPoint];
              }

              // go through slice to calc avg
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // calc avg
                avg += (individuals[curr_tree_index][gridPoint] * values[curr_tree_index][gridPoint]) / sum_ind;
              }

              // go through slice to update values
              for (int j = 0; j < j_dim; ++j)
              {
                gridPoint.back() = j;

                // update values of current slice
                values[curr_tree_index][gridPoint] -= avg;
              }

              // update lower dimensional tree
              gridPoint.pop_back();
              values[tree_index][gridPoint] += avg;
            }
          }

          ++curr_dim_index;
        }

        // update tolerance
        if (variables[curr_tree_index].size() == 1)
        {
          tol[curr_tree_index] = 1; // todo
        }
        else
        {
          tol[curr_tree_index] = 1;
        }

        ++iter;
      }

      --curr_tree_index;
    }

    // ------------- attach to rpf class -------------

    // fill with new trees
    for (size_t tree_index = 0; tree_index < variables.size(); ++tree_index)
    {
      LeafGrid curr_gridLeaf;
      curr_gridLeaf.grid = grids[tree_index];
      curr_gridLeaf.individuals = individuals[tree_index];
      curr_gridLeaf.lim_list = lim_list;
      curr_gridLeaf.values = values[tree_index];
      curr_family[variables[tree_index]]->GridLeaves = curr_gridLeaf;
    }
  }

  purified = true;
}

void RandomPlantedForest::purify_3()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list(feature_size);

    // go through all variables of the component
    for (int curr_dim = 1; curr_dim <= feature_size; ++curr_dim)
    {
      std::vector<double> bounds;

      // go through trees of family
      for (const auto &curr_tree : curr_family)
      {

        // consider only relevant trees that have current dimension as variable
        if (!curr_tree.first.count(curr_dim))
          continue;

        // go through leaves of tree
        for (const auto &curr_leaf : curr_tree.second->leaves)
        {
          // get interval ends of variable
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].first);
          bounds.push_back(curr_leaf.intervals[curr_dim - 1].second);
        }
      }
      std::sort(bounds.begin(), bounds.end());
      bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
      // int i_last = bounds.size()-1;
      // double bibi = bounds[i_last] + 0.0001;
      // bounds[i_last] = bounds[i_last] + 0.0001;
      lim_list[curr_dim - 1] = bounds;
    }

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids(curr_family.size() - 1);
    std::vector<utils::Matrix<int>> individuals(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values(curr_family.size() - 1);
    std::vector<utils::Matrix<std::vector<double>>> values_old(curr_family.size() - 1);
    std::vector<std::set<int>> variables(curr_family.size() - 1);

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {

      if (curr_tree.first == std::set<int>{0})
      {

        // values[tree_index] = rpf::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
        continue; // ignore null tree
      }

      // fill space with dimensions
      std::vector<int> dimensions;
      for (const auto &dim : curr_tree.first)
      {
        dimensions.push_back(lim_list[dim - 1].size()); // size - 1 ? WICHTIG
      }

      // setup grid for leaf indices
      auto grid = grid::NDGrid(dimensions);

      // initialize data for current tree
      grids[tree_index] = grid;
      individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
      values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));     // changed
      values_old[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0)); // changed
      variables[tree_index] = curr_tree.first;

      // fill grid points with individuals and values
      while (!grid.nextPoint())
      {

        std::vector<int> gridPoint = grid.getPoint();

        bool in_leaf = true;

        // go through sample points to sum up individuals
        for (const auto &feature_vec : X)
        {
          int dim_index = 0;
          in_leaf = true;
          for (const auto &dim : curr_tree.first)
          {
            double val = feature_vec[dim - 1];
            if (!((val >= lim_list[dim - 1][gridPoint[dim_index]]) && (val < lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // consider individuals only if all in
          if (in_leaf)
            individuals[tree_index][gridPoint] += 1;
        }

        // go through leaves of tree to sum up values
        for (const auto &leaf : curr_tree.second->get_leaves())
        {

          in_leaf = true;
          int dim_index = 0;
          for (const auto &dim : curr_tree.first)
          {
            // consider values only if all in
            if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
              in_leaf = false;
            ++dim_index;
          }

          // sum up values
          if (in_leaf)
          {

            values[tree_index][gridPoint] += leaf.value;     // todo: multiclass
            values_old[tree_index][gridPoint] += leaf.value; // todo: multiclass
          }
        }
      }

      ++tree_index;
    }

    // Rcout << variables.size();
    // for(int i = 0; i<variables.size(); ++i){
    //
    //   // Rcout << variables[i].size();
    //
    //   for(auto dim: variables[i]) Rcout << dim << ",";
    //
    //   //  Rcout << variables[i][j] << ",";
    //   //}
    //
    //   Rcout << std::endl;
    // }

    // ------------- create new trees -------------

    // insert null tree
    grids.insert(grids.begin(), grid::NDGrid());
    values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    values_old.insert(values_old.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
    individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
    variables.insert(variables.begin(), std::set<int>{0});

    // recap maximum number of dimensions of current family
    unsigned int curr_max = curr_family.rbegin()->first.size();

    while (curr_max > 1)
    {

      auto keys = getKeys(curr_family);
      // go through split dimensions of all trees
      for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
      {
        auto &curr_tree = curr_family[(*key)];
        std::set<int> curr_dims = curr_tree->split_dims;
        // check if number of dims same as current max_interaction
        if (curr_dims.size() == curr_max)
        {
          // go through feature dims
          int dim_index = 0;
          for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
          {
            // continue only if dim in current tree
            if (curr_tree->split_dims.count(feature_dim) != 0)
            {
              std::set<int> tree_dims = curr_tree->split_dims;
              tree_dims.erase(tree_dims.find(feature_dim)); // remove current feature dim from current tree
              // check if tree with dimensions exists, if not create
              std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
              if (!tree)
              {
                // get index of old and new tree
                auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
                curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
                auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
                // remove matrix dimension of respective variable
                std::vector<int> matrix_dimensions = values[old_tree_index].dims;
                // std::vector<int> matrix_dimensions = values_old[old_tree_index].dims;

                // Rcout << typeof(matrix_dimensions.begin()) << std::endl;

                matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);
                // initialize data for new tree
                auto grid = grid::NDGrid(matrix_dimensions);
                grids.insert(grids.begin() + tree_index, grid);
                values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
                values_old.insert(values_old.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
                individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
                variables.insert(variables.begin() + tree_index, tree_dims);
                // fill individuals of new trees
                while (!grid.nextPoint())
                {
                  std::vector<int> gridPoint = grid.getPoint();
                  bool in_leaf = true;
                  // go through sample points to sum up individuals
                  for (const auto &feature_vec : X)
                  {
                    int dim_index2 = 0;
                    in_leaf = true;
                    for (const auto &dim : tree_dims)
                    {
                      double val = feature_vec[dim - 1];
                      if (!((val >= lim_list[dim - 1][gridPoint[dim_index2]]) && (val < lim_list[dim - 1][gridPoint[dim_index2] + 1])))
                        in_leaf = false;
                      ++dim_index2;
                    }
                    // consider individuals only if all in
                    if (in_leaf)
                      individuals[tree_index][gridPoint] += 1;
                  }
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

    // Rcout << std::endl;
    // Rcout << std::endl;
    // Rcout << std::endl;
    //
    // for(int i = 0; i<variables.size(); ++i){
    //
    //   // Rcout << variables[i].size();
    //
    //   for(auto dim: variables[i]) Rcout << dim << ",";
    //
    //   //  Rcout << variables[i][j] << ",";
    //   //}
    //
    //   Rcout << std::endl;
    // }

    // ------------- purify -------------
    // iterate backwards through tree family
    int tree_index_t = curr_family.size() - 1;
    for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t)
    {
      std::set<int> curr_dims = *tree_t;
      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;
      // Rcout << std::endl << tree_index_t << " - T: ";
      //  Rcout << "tree_t:";
      //  for(auto dim: curr_dims) Rcout << dim << ", ";
      //  Rcout << std::endl;

      auto grid = grids[tree_index_t];
      //     Rcout << "Grid dimensions of T: ";
      //     for(auto dim: grid.dimensions) Rcout << dim << ", ";
      //     Rcout << std::endl;
      // go through subtrees of t
      int tree_index_u = variables.size();
      for (auto tree_u = variables.rbegin(); tree_u != variables.rend(); ++tree_u)
      {
        --tree_index_u;
        // j_dims = dims of t without u
        std::set<int> j_dims = curr_dims;
        if (tree_u->size() > curr_dims.size())
          continue;
        // check if subset
        bool subset = true;
        for (const auto dim : *tree_u)
        {
          if (tree_t->count(dim) == 0)
          {
            subset = false;
            break;
          }
          j_dims.erase(dim);
        }
        if (!subset)
          continue;

        // Rcout << "Hello";
        // Rcout << "   " << tree_index_u << " - U: ";
        // for(auto dim: *tree_u) Rcout << dim << ", ";
        // Rcout << std::endl;
        // Rcout << "   Individuals: ";

        double tot_sum = 0;
        grid = grids[tree_index_u];
        while (!grid.nextPoint())
        {
          auto gridPoint = grid.getPoint();
          //     Rcout << individuals[tree_index_u][gridPoint] << ", ";
          tot_sum += individuals[tree_index_u][gridPoint];
        }
        // Rcout << "Total sum: " << tot_sum << std::endl;
        // Rcout << std::endl;

        grid = grids[tree_index_u];
        //     Rcout << "      Grid dimensions of U: ";
        //     for(auto dim: grid.dimensions) Rcout << dim << ", ";
        //     Rcout << std::endl;

        // Rcout<< "j_dims: "<<j_dims.size() << std::endl;;

        std::vector<double> update(value_size, 0);

        if (j_dims.size() == 0)
        {

          // grid = grids[tree_index_u];
          while (!grid.nextPoint())
          {
            auto gridPoint_i = grid.getPoint();
            //     Rcout << "         " << "i: ";
            //     for(auto p: gridPoint_i) Rcout << p << ", ";
            //     Rcout << std::endl << "         ";
            double curr_sum = individuals[tree_index_u][gridPoint_i];
            //     Rcout << ", Current Sum: " << curr_sum << std::endl;
            //     Rcout << std::endl << "         " << "i, j: ";
            update += (curr_sum / tot_sum) * values_old[tree_index_t][gridPoint_i];
            //     Rcout << std::endl;
          }

          int tree_index_s = variables.size();
          for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
          {

            // Rcout << "tree_s:";
            // for(auto dim: *tree_s) Rcout << dim << ", ";
            // Rcout << std::endl;

            --tree_index_s;
            if (*tree_s == std::set<int>{0})
            {

              auto gridPoint_0 = std::vector<int>{0};
              values[tree_index_s][gridPoint_0] += update;
              //     Rcout << std::endl;
              //}

              /*
               for(auto tree_0: curr_family){

               if(tree_0.first == std::set<int>{0}){

               Rcout << tree_0.first.size();
               std::vector<int> leaf_index(tree_0.first.size(), 0);
               std::vector<int> leaf_index(tree_0.second->GridLeaves.values.size(), 0);

               int Test = tree_0.second->GridLeaves.values.size();
               Rcout << Test;
               tree_0.second->GridLeaves.values[leaf_index] += update;
               }
               }
               */
            }
            else
            {

              // check if S subset of T

              bool subset = true;
              for (const auto dim : *tree_s)
              {
                if (tree_t->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              if (!subset)
                continue;

              // Rcout << pow(-1, (*tree_s).size()) << std::endl;

              auto grid_k = grids[tree_index_s];
              while (!grid_k.nextPoint())
              {
                auto gridPoint_k = grid_k.getPoint();
                //
                //      if((*tree_s).size()>2){
                //      Rcout << std::endl << "            " << "j, k: ";
                //      for(auto p: gridPoint_k) Rcout << p << ", ";
                //      Rcout << std::endl;
                //      }
                //
                //      Rcout << pow(-1, (*tree_s).size()) * update << std::endl;
                values[tree_index_s][gridPoint_k] += pow(-1, (*tree_s).size()) * update;
              }
            }
          }
          // Rcout << std::endl;
        }
        else
        {

          std::vector<int> j_sizes(j_dims.size(), 0);
          for (size_t j = 0; j < j_dims.size(); ++j)
          {
            auto tmp = j_dims.begin();
            std::advance(tmp, j);
            int j_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*tmp));
            j_sizes[j] = grids[tree_index_t].dimensions[j_index];
          }

          // Rcout<<"Hello 1";

          grid::NDGrid grid_j = grid::NDGrid(j_sizes);
          while (!grid_j.nextPoint())
          {

            std::vector<double> update(value_size, 0);
            auto gridPoint_j = grid_j.getPoint();
            //     Rcout << "         " << "j: ";
            //     for(auto p: gridPoint_j) Rcout << p << ", ";
            //     Rcout << std::endl;
            // calc update
            grid = grids[tree_index_u];
            while (!grid.nextPoint())
            {
              auto gridPoint_i = grid.getPoint();
              //     Rcout << "         " << "i: ";
              //     for(auto p: gridPoint_i) Rcout << p << ", ";
              //     Rcout << std::endl << "         ";
              double curr_sum = individuals[tree_index_u][gridPoint_i];
              //     Rcout << ", Current Sum: " << curr_sum << std::endl;
              std::vector<int> gridPoint_ij(tree_t->size(), 0);
              for (size_t j = 0; j < gridPoint_j.size(); ++j)
              {
                auto j_dim = j_dims.begin();
                std::advance(j_dim, j);
                int j_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*j_dim));
                //     Rcout << "         j_dim=" << *j_dim << ", j_index=" << j_index;
                gridPoint_ij[j_index] = gridPoint_j[j];
              }
              for (size_t i = 0; i < gridPoint_i.size(); ++i)
              {
                auto i_dim = tree_u->begin();
                std::advance(i_dim, i);
                int i_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*i_dim));
                //     Rcout << "         i_dim=" << *i_dim << ", i_index=" << i_index;
                gridPoint_ij[i_index] = gridPoint_i[i];
              }
              //     Rcout << std::endl << "         " << "i, j: ";
              //     for(auto p: gridPoint_ij) Rcout << p << ", ";
              //     Rcout << std::endl;
              update += (curr_sum / tot_sum) * values_old[tree_index_t][gridPoint_ij];
              //     Rcout << std::endl;
            }

            // Rcout << "Hello_2";
            // update trees
            int tree_index_s = variables.size();
            for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
            {
              --tree_index_s;
              // check if T\U=j_dims subset of S and S subset of T
              bool subset = true;
              for (const auto dim : j_dims)
              {
                if (tree_s->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              for (const auto dim : *tree_s)
              {
                if (tree_t->count(dim) == 0)
                {
                  subset = false;
                  break;
                }
              }
              if (!subset)
                continue;
              //     Rcout << "         " << "S: ";
              //     for(auto dim: *tree_s) Rcout << dim << ", ";
              //     Rcout << std::endl;
              // S cap U
              std::set<int> k_dims = *tree_s;
              std::set<int> k_dims_h1 = *tree_s;
              std::set<int> k_dims_h2 = *tree_u;
              for (const auto dim : *tree_u)
                k_dims.insert(dim);
              for (const auto dim : *tree_s)
                k_dims_h2.erase(dim);
              for (const auto dim : *tree_u)
                k_dims_h1.erase(dim);
              for (const auto dim : k_dims_h1)
                k_dims.erase(dim);
              for (const auto dim : k_dims_h2)
                k_dims.erase(dim);

              // std::set<int> k_dims = *tree_s;
              // for(const auto dim: *tree_t) k_dims.erase(dim);
              // for(const auto dim: *tree_u) k_dims.insert(dim);

              //     Rcout << "         " << "k_dims: ";
              //     for(auto dim: k_dims) Rcout << dim << ", ";
              //     Rcout << std::endl;

              if (k_dims.size() == 0)
              {

                values[tree_index_s][gridPoint_j] += pow(-1, (*tree_s).size() - j_dims.size()) * update;
              }
              else
              {

                // Rcout <<"k_dims :";
                // for(auto dim: k_dims) Rcout << dim << ", ";
                // Rcout << std::endl;

                std::vector<int> k_sizes(k_dims.size(), 0);
                for (size_t k = 0; k < k_dims.size(); ++k)
                {
                  auto tmp = k_dims.begin();
                  std::advance(tmp, k);
                  int k_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*tmp));
                  k_sizes[k] = grids[tree_index_t].dimensions[k_index];
                }
                // Rcout << "         " << "k_sizes: ";
                // for(auto dim: k_sizes) Rcout << dim << ", ";
                // Rcout << std::endl;
                grid::NDGrid grid_k = grid::NDGrid(k_sizes);
                while (!grid_k.nextPoint())
                {
                  auto gridPoint_k = grid_k.getPoint();
                  // Rcout << "            " << "k: ";
                  // for(auto p: gridPoint_k) Rcout << p << ", ";
                  // Rcout << std::endl << "         ";
                  std::vector<int> gridPoint_jk(tree_s->size(), 0);
                  for (size_t j = 0; j < gridPoint_j.size(); ++j)
                  {
                    auto j_dim = j_dims.begin();
                    std::advance(j_dim, j);
                    int j_index = std::distance(variables[tree_index_s].begin(), variables[tree_index_s].find(*j_dim));
                    // Rcout << "         j_dim=" << *j_dim << ", j_index=" << j_index;
                    gridPoint_jk[j_index] = gridPoint_j[j];
                  }
                  for (size_t k = 0; k < gridPoint_k.size(); ++k)
                  {
                    auto k_dim = k_dims.begin();
                    std::advance(k_dim, k);
                    int k_index = std::distance(variables[tree_index_s].begin(), variables[tree_index_s].find(*k_dim));
                    // Rcout << "         k_dim=" << *k_dim << ", k_index=" << k_index;
                    gridPoint_jk[k_index] = gridPoint_k[k];
                  }
                  // Rcout << std::endl << "            " << "j, k: ";
                  // for(auto p: gridPoint_jk) Rcout << p << ", ";
                  // Rcout << std::endl;

                  // Rcout << pow(-1, (*tree_s).size() - j_dims.size()) * update[0];
                  values[tree_index_s][gridPoint_jk] += pow(-1, (*tree_s).size() - j_dims.size()) * update;
                }
              }
            }
          }
        }
      }
      --tree_index_t;
    }

    // ------------- attach to rpf class -------------

    // fill with new trees
    for (size_t tree_index = 0; tree_index < variables.size(); ++tree_index)
    {
      LeafGrid curr_gridLeaf;
      curr_gridLeaf.grid = grids[tree_index];
      curr_gridLeaf.individuals = individuals[tree_index];
      curr_gridLeaf.lim_list = lim_list;
      curr_gridLeaf.values = values[tree_index];
      curr_family[variables[tree_index]]->GridLeaves = curr_gridLeaf;
    }
  }

  purified = true;
}

void RandomPlantedForest::print()
{
  for (int n = 0; n < n_trees; ++n)
  {
    TreeFamily family = tree_families[n];
    auto keys = getKeys(family);
    for (size_t m = 0; m < keys.size(); ++m)
    {
      DecisionTree tree = *(family[keys[m]]);
      Rcout << m + 1 << " Tree: ";
      Rcout << "Dims=";
      for (const auto &dim : tree.split_dims)
        Rcout << dim << ",";
      Rcout << std::endl
            << "Leaves: (" << tree.leaves.size() << ")" << std::endl;
      for (const auto &leaf : tree.leaves)
      {
        Rcout << "Intervals=";
        for (const auto &interval : leaf.intervals)
        {
          Rcout << interval.first << "," << interval.second << "/";
        }
        Rcout << " Value=";
        for (const auto &val : leaf.value)
          Rcout << val << ", ";
        Rcout << std::endl;
      }
      Rcout << std::endl;
    }
    Rcout << std::endl
          << std::endl;
  }
}

// print parameters of the model to the console
void RandomPlantedForest::get_parameters()
{
  Rcout << "Parameters: n_trees=" << n_trees << ", n_splits=" << n_splits << ", max_interaction=" << max_interaction << ", t_try=" << t_try
        << ", split_try=" << split_try << ", purified=" << purified << ", deterministic=" << deterministic << ", nthreads=" << nthreads
        << ", feature_size=" << feature_size << ", sample_size=" << sample_size << std::endl;
}

/*  retrospectively change parameters of existing class object,
 updates the model, so far only single valued parameters supported,
 for replacing training data use 'set_data',
 note that changing cv does not trigger cross validation */
void RandomPlantedForest::set_parameters(StringVector keys, NumericVector values)
{
  if (keys.size() != values.size())
  {
    Rcout << "Size of input vectors is not the same. " << std::endl;
    return;
  }

  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    if (keys[i] == "deterministic")
    {
      this->deterministic = values[i];
    }
    else if (keys[i] == "nthreads")
    {
      this->nthreads = values[i];
    }
    else if (keys[i] == "purify")
    {
      this->purify_forest = values[i];
    }
    else if (keys[i] == "n_trees")
    {
      this->n_trees = values[i];
    }
    else if (keys[i] == "n_splits")
    {
      this->n_splits = values[i];
    }
    else if (keys[i] == "t_try")
    {
      this->t_try = values[i];
    }
    else if (keys[i] == "split_try")
    {
      this->split_try = values[i];
    }
    else if (keys[i] == "max_interaction")
    {
      this->max_interaction = values[i];
    }
    else if (keys[i] == "cv")
    {
      this->cross_validate = values[i];
    }
    else
    {
      Rcout << "Unkown parameter key  '" << keys[i] << "' ." << std::endl;
    }
  }
  this->fit();
}

List RandomPlantedForest::get_model()
{
  List model;
  for (const auto &family : tree_families)
  {
    List variables, family_values, family_intervals;
    for (const auto &tree : family)
    {
      List tree_values;
      List tree_intervals;
      variables.push_back(from_std_set(tree.first));
      for (const auto &leaf : tree.second->leaves)
      {
        NumericMatrix leaf_values;
        for (const auto &val : leaf.value)
        {
          leaf_values.push_back(val);
        }
        tree_values.push_back(leaf_values);

        NumericVector intervals;
        for (const auto &interval : leaf.intervals)
        {
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
  return (model);
}
