#include "rpf.hpp"

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

// Extracted body to allow multithreading over families (now exposed as purify_3(TreeFamily&))
void RandomPlantedForest::purify_3(TreeFamily &curr_family)
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
      continue; // ignore null tree
    }

    // fill space with dimensions
    std::vector<int> dimensions;
    dimensions.reserve(curr_tree.first.size());
    for (const auto &dim : curr_tree.first)
    {
      dimensions.push_back(lim_list[dim - 1].size());
    }

    // setup grid for leaf indices
    auto grid = grid::NDGrid(dimensions);

    // initialize data for current tree
    grids[tree_index] = grid;
    individuals[tree_index] = utils::Matrix<int>(dimensions, 0);
    values[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));
    values_old[tree_index] = utils::Matrix<std::vector<double>>(dimensions, std::vector<double>(value_size, 0));
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
          {
            in_leaf = false;
            break;
          }
          ++dim_index;
        }
        // consider individuals only if all in
        if (in_leaf)
          individuals[tree_index][gridPoint] += 1;
      }

      // go through leaves of tree to sum up values
      const auto leaves_once = curr_tree.second->get_leaves();
      for (const auto &leaf : leaves_once)
      {
        in_leaf = true;
        int dim_index = 0;
        for (const auto &dim : curr_tree.first)
        {
          // consider values only if all in
          if (!((leaf.intervals[dim - 1].first <= lim_list[dim - 1][gridPoint[dim_index]]) && (leaf.intervals[dim - 1].second >= lim_list[dim - 1][gridPoint[dim_index] + 1])))
          {
            in_leaf = false;
            break;
          }
          ++dim_index;
        }
        // sum up values
        if (in_leaf)
        {
          values[tree_index][gridPoint] += leaf.value;
          values_old[tree_index][gridPoint] += leaf.value;
        }
      }
    }

    ++tree_index;
  }

  // ------------- create new trees -------------
  grids.insert(grids.begin(), grid::NDGrid());
  values.insert(values.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
  values_old.insert(values_old.begin(), utils::Matrix<std::vector<double>>(std::vector<int>{1}, std::vector<double>(value_size, 0)));
  individuals.insert(individuals.begin(), utils::Matrix<int>(std::vector<int>{1}));
  variables.insert(variables.begin(), std::set<int>{0});

  unsigned int curr_max = curr_family.rbegin()->first.size();
  while (curr_max > 1)
  {
    auto keys = getKeys(curr_family);
    for (std::vector<std::set<int>>::reverse_iterator key = keys.rbegin(); key != keys.rend(); ++key)
    {
      auto &curr_tree = curr_family[(*key)];
      std::set<int> curr_dims = curr_tree->split_dims;
      if (curr_dims.size() == curr_max)
      {
        int dim_index = 0;
        for (int feature_dim = 1; feature_dim <= feature_size; ++feature_dim)
        {
          if (curr_tree->split_dims.count(feature_dim) != 0)
          {
            std::set<int> tree_dims = curr_tree->split_dims;
            tree_dims.erase(tree_dims.find(feature_dim));
            std::shared_ptr<DecisionTree> tree = treeExists(tree_dims, curr_family);
            if (!tree)
            {
              auto old_tree_index = std::distance(std::begin(curr_family), curr_family.find(curr_tree->get_split_dims()));
              curr_family.insert(std::make_pair(tree_dims, std::make_shared<DecisionTree>(DecisionTree(tree_dims))));
              auto tree_index = std::distance(std::begin(curr_family), curr_family.find(tree_dims));
              std::vector<int> matrix_dimensions = values[old_tree_index].dims;
              matrix_dimensions.erase(matrix_dimensions.begin() + dim_index);
              auto grid = grid::NDGrid(matrix_dimensions);
              grids.insert(grids.begin() + tree_index, grid);
              values.insert(values.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
              values_old.insert(values_old.begin() + tree_index, utils::Matrix<std::vector<double>>(matrix_dimensions, std::vector<double>(value_size, 0)));
              individuals.insert(individuals.begin() + tree_index, utils::Matrix<int>(matrix_dimensions));
              variables.insert(variables.begin() + tree_index, tree_dims);
              while (!grid.nextPoint())
              {
                std::vector<int> gridPoint = grid.getPoint();
                bool in_leaf = true;
                for (const auto &feature_vec : X)
                {
                  int dim_index2 = 0;
                  in_leaf = true;
                  for (const auto &dim : tree_dims)
                  {
                    double val = feature_vec[dim - 1];
                    if (!((val >= lim_list[dim - 1][gridPoint[dim_index2]]) && (val < lim_list[dim - 1][gridPoint[dim_index2] + 1])))
                    {
                      in_leaf = false;
                      break;
                    }
                    ++dim_index2;
                  }
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
    --curr_max;
  }

  // ------------- purify -------------
  std::vector<std::vector<int>> dim_to_pos(variables.size(), std::vector<int>(feature_size + 1, -1));
  for (size_t idx = 0; idx < variables.size(); ++idx)
  {
    int pos = 0;
    for (const auto dim : variables[idx])
    {
      if (dim >= 0 && dim <= feature_size) dim_to_pos[idx][dim] = pos++;
    }
  }

  std::vector<double> total_individuals(variables.size(), 0.0);
  for (size_t idx = 0; idx < variables.size(); ++idx)
  {
    double tot = 0.0;
    if (variables[idx] == std::set<int>{0})
    {
      std::vector<int> only{0};
      tot += individuals[idx][only];
    }
    else
    {
      auto grid_sum = grids[idx];
      while (!grid_sum.nextPoint())
      {
        auto gp = grid_sum.getPoint();
        tot += individuals[idx][gp];
      }
    }
    total_individuals[idx] = tot;
  }

  int tree_index_t = curr_family.size() - 1;
  for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t)
  {
    std::set<int> curr_dims = *tree_t;
    if (curr_dims == std::set<int>{0})
      continue;

    auto grid = grids[tree_index_t];
    int tree_index_u = variables.size();
    for (auto tree_u = variables.rbegin(); tree_u != variables.rend(); ++tree_u)
    {
      --tree_index_u;
      std::set<int> j_dims = curr_dims;
      if (tree_u->size() > curr_dims.size())
        continue;
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

      double tot_sum = total_individuals[tree_index_u];
      if (tot_sum == 0.0)
        continue;
      const double inv_tot_sum = 1.0 / tot_sum;

      grid = grids[tree_index_u];
      std::vector<double> update(value_size, 0);

      if (j_dims.size() == 0)
      {
        while (!grid.nextPoint())
        {
          auto gridPoint_i = grid.getPoint();
          double curr_sum = individuals[tree_index_u][gridPoint_i];
          update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_i];
        }

        int tree_index_s = variables.size();
        for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
        {
          --tree_index_s;
          if (*tree_s == std::set<int>{0})
          {
            auto gridPoint_0 = std::vector<int>{0};
            values[tree_index_s][gridPoint_0] += update;
          }
          else
          {
            bool subset2 = true;
            for (const auto dim : *tree_s)
            {
              if (tree_t->count(dim) == 0)
              {
                subset2 = false;
                break;
              }
            }
            if (!subset2)
              continue;
            auto grid_k = grids[tree_index_s];
            while (!grid_k.nextPoint())
            {
              auto gridPoint_k = grid_k.getPoint();
              int sign0 = ((*tree_s).size() % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_k] += sign0 * update;
            }
          }
        }
      }
      else
      {
        std::vector<int> j_sizes(j_dims.size(), 0);
        for (size_t j = 0; j < j_dims.size(); ++j)
        {
          auto tmp = j_dims.begin();
          std::advance(tmp, j);
          int j_index = dim_to_pos[tree_index_t][*tmp];
          j_sizes[j] = grids[tree_index_t].dimensions[j_index];
        }
        grid::NDGrid grid_j = grid::NDGrid(j_sizes);
        while (!grid_j.nextPoint())
        {
          std::vector<double> update(value_size, 0);
          auto gridPoint_j = grid_j.getPoint();
          grid = grids[tree_index_u];
          while (!grid.nextPoint())
          {
            auto gridPoint_i = grid.getPoint();
            double curr_sum = individuals[tree_index_u][gridPoint_i];
            std::vector<int> gridPoint_ij(tree_t->size(), 0);
            for (size_t j = 0; j < gridPoint_j.size(); ++j)
            {
              auto j_dim = j_dims.begin();
              std::advance(j_dim, j);
              int j_index = dim_to_pos[tree_index_t][*j_dim];
              gridPoint_ij[j_index] = gridPoint_j[j];
            }
            for (size_t i = 0; i < gridPoint_i.size(); ++i)
            {
              auto i_dim = tree_u->begin();
              std::advance(i_dim, i);
              int i_index = dim_to_pos[tree_index_t][*i_dim];
              gridPoint_ij[i_index] = gridPoint_i[i];
            }
            update += (curr_sum * inv_tot_sum) * values_old[tree_index_t][gridPoint_ij];
          }

          int tree_index_s = variables.size();
          for (auto tree_s = variables.rbegin(); tree_s != variables.rend(); ++tree_s)
          {
            --tree_index_s;
            bool subset2 = true;
            for (const auto dim : j_dims)
            {
              if (tree_s->count(dim) == 0)
              {
                subset2 = false;
                break;
              }
            }
            for (const auto dim : *tree_s)
            {
              if (tree_t->count(dim) == 0)
              {
                subset2 = false;
                break;
              }
            }
            if (!subset2)
              continue;

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

            if (k_dims.size() == 0)
            {
              size_t diff = (*tree_s).size() - j_dims.size();
              int sign = (diff % 2 == 0) ? 1 : -1;
              values[tree_index_s][gridPoint_j] += sign * update;
            }
            else
            {
              std::vector<int> k_sizes(k_dims.size(), 0);
              for (size_t k = 0; k < k_dims.size(); ++k)
              {
                auto tmp = k_dims.begin();
                std::advance(tmp, k);
                int k_index = dim_to_pos[tree_index_t][*tmp];
                k_sizes[k] = grids[tree_index_t].dimensions[k_index];
              }
              grid::NDGrid grid_k = grid::NDGrid(k_sizes);
              while (!grid_k.nextPoint())
              {
                auto gridPoint_k = grid_k.getPoint();
                std::vector<int> gridPoint_jk(tree_s->size(), 0);
                for (size_t j = 0; j < gridPoint_j.size(); ++j)
                {
                  auto j_dim = j_dims.begin();
                  std::advance(j_dim, j);
                  int j_index = dim_to_pos[tree_index_s][*j_dim];
                  gridPoint_jk[j_index] = gridPoint_j[j];
                }
                for (size_t k = 0; k < gridPoint_k.size(); ++k)
                {
                  auto k_dim = k_dims.begin();
                  std::advance(k_dim, k);
                  int k_index = dim_to_pos[tree_index_s][*k_dim];
                  gridPoint_jk[k_index] = gridPoint_k[k];
                }
                size_t diff = (*tree_s).size() - j_dims.size();
                int sign2 = (diff % 2 == 0) ? 1 : -1;
                values[tree_index_s][gridPoint_jk] += sign2 * update;
              }
            }
          }
        }
      }
    }
    --tree_index_t;
  }

  // ------------- attach to rpf class -------------
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
            {
              in_leaf = false;
              break;
            }
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
            {
              in_leaf = false;
              break;
            }
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
  

  // If threading is enabled, parallelize across tree families
  unsigned int threads_to_use = static_cast<unsigned int>(nthreads);
  if (threads_to_use == 0) threads_to_use = 1;
  if (threads_to_use > 1)
  {
    if (threads_to_use > std::thread::hardware_concurrency())
    {
      Rcout << "Requested " << threads_to_use << " threads but only " << std::thread::hardware_concurrency() << " available" << std::endl;
    }
    for (size_t start = 0; start < this->tree_families.size(); start += (size_t)threads_to_use)
    {
      size_t batch = std::min<size_t>((size_t)threads_to_use, this->tree_families.size() - start);
      if (batch == 0) break;
      std::vector<std::thread> threads(batch);
      for (size_t i = 0; i < batch; ++i)
      {
        size_t fam_index = start + i;
        threads[i] = std::thread([this](TreeFamily *fam_ptr){
          this->purify_3(*fam_ptr);
        }, &this->tree_families[fam_index]);
      }
      for (auto &th : threads)
      {
        if (th.joinable()) th.join();
      }
    }
    purified = true;
    return;
  }

  // Serial path: reuse the per-family overload for identical behavior
  for (auto &curr_family : this->tree_families)
  {
    this->purify_3(curr_family);
  }
  purified = true;
  return;
}
