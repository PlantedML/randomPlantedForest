# print() summarises a regression fit and its purification state

    Code
      print(fit)
    Output
      -- Regression Random Planted Forest --------------------------------------------
      Formula: `mpg ~ cyl + wt`
      3 tree families with 30 splits each on 2 predictors, main effects only.
      i Forest is not purified.
      
      -- Tree growing 
         split_structure: leaves
               split_try: 10
                   t_try: 0.4
          max_candidates: 50
        split_decay_rate: 0.1
           delete_leaves: TRUE
      
      i Fit using 1 thread, also the default for `predict()` and `purify()`.

---

    Code
      print(fit)
    Output
      -- Regression Random Planted Forest --------------------------------------------
      Formula: `mpg ~ cyl + wt`
      3 tree families with 30 splits each on 2 predictors, main effects only.
      v Forest is purified.
      
      -- Tree growing 
         split_structure: leaves
               split_try: 10
                   t_try: 0.4
          max_candidates: 50
        split_decay_rate: 0.1
           delete_leaves: TRUE
      
      i Fit using 1 thread, also the default for `predict()` and `purify()`.

# print() covers x/y fits, classification losses and deterministic fits

    Code
      print(fit_logit)
    Output
      -- Classification Random Planted Forest ----------------------------------------
      Predictors: `Sepal.Length`, `Sepal.Width`, `Petal.Length`, and `Petal.Width`
      1 tree family with 30 splits each on 4 predictors, main effects only.
      i Forest is not purified.
      ! Fit deterministically.
      
      -- Tree growing 
         split_structure: leaves
               split_try: 10
                   t_try: 0.4
          max_candidates: 50
        split_decay_rate: 0.1
           delete_leaves: TRUE
      
      -- Loss 
           loss: logit
          delta: 0.001
        epsilon: 0.1
      
      i Fit using 1 thread, also the default for `predict()` and `purify()`.

---

    Code
      print(fit_l2)
    Output
      -- Classification Random Planted Forest ----------------------------------------
      Formula: `Species ~ .`
      2 tree families with 30 splits each on 4 predictors, main effects only.
      i Forest is not purified.
      
      -- Tree growing 
         split_structure: leaves
               split_try: 10
                   t_try: 0.4
          max_candidates: 50
        split_decay_rate: 0.1
           delete_leaves: TRUE
      
      -- Loss 
        loss: L2
      
      i Fit using 1 thread, also the default for `predict()` and `purify()`.

# an exported forest prints compactly

    Code
      print(fit$forest)
    Output
      <rpf_forest> of 3 trees
    Code
      str(fit$forest)
    Output
      <rpf_forest> of 3 trees

