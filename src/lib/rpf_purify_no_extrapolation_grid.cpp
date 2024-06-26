#include "rpf.hpp"

void RandomPlantedForest::purify_no_extrapolation_existing_grid()
{

  // go through all n_trees families
  for (auto &curr_family : this->tree_families)
  {

    // lim_list is a list giving for each variable all interval end-points
    std::vector<std::vector<double>> lim_list = get_lim_list(curr_family);

    // initialize values and individuals for each tree in family
    std::vector<grid::NDGrid> grids;
    std::vector<utils::Matrix<int>> individuals;
    std::vector<utils::Matrix<std::vector<double>>> values;
    std::vector<utils::Matrix<std::vector<double>>> values_old;
    std::vector<std::set<int>> variables;

    //  ------------- setup finer grid  -------------

    int tree_index = 0;
    for (const auto &curr_tree : curr_family)
    {
      grids.push_back(curr_tree.second->GridLeaves.grid);
      individuals.push_back(curr_tree.second->GridLeaves.individuals);
      values.push_back(curr_tree.second->GridLeaves.values);
      values_old.push_back(curr_tree.second->GridLeaves.values);
      variables.push_back(curr_tree.first);
      ++tree_index;
    }



    // recap maximum number of dimensions of current family
    unsigned int curr_max = curr_family.rbegin()->first.size();

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
        // while (!grid.nextPoint())
        // {
        //   auto gridPoint = grid.getPoint();
        //   //     Rcout << individuals[tree_index_u][gridPoint] << ", ";
        //   tot_sum += individuals[tree_index_u][gridPoint];
        // }
        // Rcout << "Total sum: " << tot_sum << std::endl;
        // Rcout << std::endl;

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
            update += curr_sum * values_old[tree_index_t][gridPoint_i];
            tot_sum += curr_sum;
            //     Rcout << std::endl;
          }
          update /= tot_sum;

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
          while (!grid_j.nextPoint()) // looping over every x_{i, T\U} essentially
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
                auto j_dim = j_dims.begin(); // to keep: T\U
                std::advance(j_dim, j);
                int j_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*j_dim));
                //     Rcout << "         j_dim=" << *j_dim << ", j_index=" << j_index;
                gridPoint_ij[j_index] = gridPoint_j[j];
              }
              for (size_t i = 0; i < gridPoint_i.size(); ++i)
              {
                auto i_dim = tree_u->begin(); // to marginalize: U
                std::advance(i_dim, i);
                int i_index = std::distance(variables[tree_index_t].begin(), variables[tree_index_t].find(*i_dim));
                //     Rcout << "         i_dim=" << *i_dim << ", i_index=" << i_index;
                gridPoint_ij[i_index] = gridPoint_i[i];
              }

              if (individuals[tree_index_t][gridPoint_ij] == 0) {
                continue; // Skip non-support areas
              }
              //     Rcout << std::endl << "         " << "i, j: ";
              //     for(auto p: gridPoint_ij) Rcout << p << ", ";
              //     Rcout << std::endl;
              update += curr_sum * values_old[tree_index_t][gridPoint_ij];
              tot_sum += curr_sum;
              //     Rcout << std::endl;
            }
            update /= tot_sum;

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
