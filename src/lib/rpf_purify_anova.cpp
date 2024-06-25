#include "rpf.hpp"

void RandomPlantedForest::purify_2()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list = get_lim_list(curr_family);

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
                      const double lower = lim_list[dim - 1][gridPoint[dim_index2]];
                      const double upper = lim_list[dim - 1][gridPoint[dim_index2] + 1];
                      if (!((val >= lower) && (val < upper)))
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

    // measure tolerance and number of iterations
    std::vector<double> tol(curr_family.size(), 1);
    int iter;

    // iterate backwards through tree family
    int tree_index_t = curr_family.size() - 1;

    for (auto tree_t = variables.rbegin(); tree_t != variables.rend(); ++tree_t)
    {
      iter = 0;
      std::set<int> curr_dims = *tree_t;

      // do not purify null
      if (curr_dims == std::set<int>{0})
        continue;
      // Rcout << std::endl << tree_index_t << " - T: ";
      //  Rcout << "tree_t:";
      //  for(auto dim: curr_dims) Rcout << dim << ", ";
      //  Rcout << std::endl;

      auto grid = grids[tree_index_t];

      // repeat until tolerance small enough and (?) maximum number of iterations reached
      while ((tol[tree_index_t] > 0.00000000001) && (iter < 100))
      {
        tol[tree_index_t] = 0;

        // go through feature dims
        int curr_dim_index = 0;
        for (const auto j : curr_dims)
        {

          // get tree that has same variables as curr_tree minus j-variable
          std::set<int> tree_dims_minusj = curr_dims;
          tree_dims_minusj.erase(tree_dims_minusj.find(j));
          int tree_index_minusj = 0; // if tree not exist, set to null tree
          if (curr_family.find(tree_dims_minusj) != curr_family.end())
            tree_index_minusj = std::distance(std::begin(curr_family), curr_family.find(tree_dims_minusj));

          // update values
          if (grids[tree_index_t].dimensions.size() == 1)
          { // one dimensional case (tree_t one dimensional)

            int sum_ind = 0;
            std::vector<double> avg(value_size, 0);

            // get sum of individuals
            for (int i = 0; i < individuals[tree_index_t].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              sum_ind += individuals[tree_index_t][tmp];
            }
            if (sum_ind == 0)
              continue;

            // calc avg
            for (int i = 0; i < individuals[tree_index_t].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              avg += (individuals[tree_index_t][tmp] * values[tree_index_t][tmp]) / sum_ind;
            }

            // update values of one dimensional and null tree
            for (int i = 0; i < values[tree_index_t].n_entries; ++i)
            {
              std::vector<int> tmp{i};
              values[tree_index_t][tmp] -= avg;
            }
            std::vector<int> tmp{0};
            values[tree_index_minusj][tmp] += avg;
          }
          else
          { // higher dimensional case (tree_t dimension >1)

            // setup new grid without dimension j
            std::vector<int> new_dimensions = grids[tree_index_t].dimensions;
            int index_j_in_t = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(j));
            std::vector<int> j_dim;
            j_dim.push_back(new_dimensions[index_j_in_t]);
            new_dimensions.erase(new_dimensions.begin() + index_j_in_t);

            grid::NDGrid grid_minusj = grid::NDGrid(new_dimensions);
            grid::NDGrid grid_minusj2 = grid::NDGrid(grids[tree_index_minusj]);

            // go through values without dimension j
            while (!grid_minusj.nextPoint())
            {

              std::vector<double> update(value_size, 0);
              auto gridPoint_minusj = grid_minusj.getPoint();
              //     Rcout << "         " << "j: ";
              //     for(auto p: gridPoint_j) Rcout << p << ", ";
              //     Rcout << std::endl;
              // calc update

              int tree_index_j = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(j));

              grid::NDGrid grid_j = grid::NDGrid(j_dim);

              double sum_ind = 0;
              std::vector<double> avg(value_size, 0);

              while (!grid_j.nextPoint())
              {
                auto gridPoint_j = grid_j.getPoint();
                //     Rcout << "         " << "i: ";
                //     for(auto p: gridPoint_i) Rcout << p << ", ";
                //     Rcout << std::endl << "         ";
                std::vector<int> gridPoint_t(tree_t->size(), 0);

                gridPoint_t[tree_index_j] = gridPoint_j[0];

                for (size_t minusj = 0; minusj < gridPoint_minusj.size(); ++minusj)
                {
                  auto minusj_dim = tree_dims_minusj.begin();
                  std::advance(minusj_dim, minusj);
                  int minusj_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*minusj_dim));
                  //     Rcout << "         i_dim=" << *i_dim << ", i_index=" << i_index;
                  gridPoint_t[minusj_index] = gridPoint_minusj[minusj];
                }
                double sum_curr = individuals[tree_index_t][gridPoint_t];
                sum_ind += sum_curr;
                avg += (sum_curr * values[tree_index_t][gridPoint_t][0]);
              }

              if (sum_ind != 0)
              {
                avg = avg / sum_ind;
              }
              else
              {
                continue;
              }
              double avgtol = avg[0];
              tol[tree_index_t] = std::max(std::fabs(avgtol), tol[tree_index_t]);
              values[tree_index_minusj][gridPoint_minusj] += avg;

              // update values for t-tree
              grid_j = grid::NDGrid(j_dim);

              while (!grid_j.nextPoint())
              {
                auto gridPoint_j = grid_j.getPoint();
                //     Rcout << "         " << "i: ";
                //     for(auto p: gridPoint_i) Rcout << p << ", ";
                //     Rcout << std::endl << "         ";
                std::vector<int> gridPoint_t(tree_t->size(), 0);

                gridPoint_t[tree_index_j] = gridPoint_j[0];

                for (size_t minusj = 0; minusj < gridPoint_minusj.size(); ++minusj)
                {
                  auto minusj_dim = tree_dims_minusj.begin();
                  std::advance(minusj_dim, minusj);
                  int minusj_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*minusj_dim));
                  //     Rcout << "         i_dim=" << *i_dim << ", i_index=" << i_index;
                  gridPoint_t[minusj_index] = gridPoint_minusj[minusj];
                }

                values[tree_index_t][gridPoint_t] -= avg;
              }

            } // end of minusj loop

          } // end of loop for higher dim t case
          ++curr_dim_index;

        } // end of j loop

        ++iter;
      } // finished tol loop

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
      auto &curr_tree = curr_family[variables[tree_index]];
      curr_tree->GridLeaves = curr_gridLeaf;
      
      auto &tree_dims = variables[tree_index];

      if (tree_dims == std::set{0}) {
        curr_tree->leaves[0].value = curr_gridLeaf.values[{0}];
        continue;
      }
      curr_tree->leaves.clear();

      while (!curr_gridLeaf.grid.nextPoint()) {
        auto grid_point = curr_gridLeaf.grid.getPoint();

        bool not_end_leaf = true;
        Leaf new_leaf;
        {
          new_leaf.value = curr_gridLeaf.values[grid_point];
          // initialize interval with split interval
          int dim_idx = 0;
          for (int dim = 1; dim <= feature_size; ++dim) {
            const int gp = grid_point[dim_idx];

            double lower, upper;
            if (tree_dims.count(dim) == 0) {
              lower = lower_bounds[dim-1];
              upper = upper_bounds[dim-1];
            } else {
              if (gp >= lim_list[dim-1].size()-1){
                not_end_leaf = false;
                break;
              }
              lower = lim_list[dim-1][gp];
              upper = lim_list[dim-1][gp + 1];
              ++dim_idx;
            }

            new_leaf.intervals.push_back(Interval{lower, upper});

          }
        }

        if (not_end_leaf)
          curr_tree->leaves.push_back(new_leaf);
      }
    }
  }

  purified = true; // debug
}
