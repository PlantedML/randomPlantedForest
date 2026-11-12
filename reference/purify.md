# Purify a Random Planted Forest

Purification makes the components of the forest's functional
decomposition unique, which
[`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
relies on. Unless
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
was called with `purify = TRUE`,
[`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
purifies the forest on first use.

## Usage

``` r
purify(x, ..., maxp_interaction = NULL, mode = 2L, nthreads = NULL)

is_purified(x)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  model.

- ...:

  Reserved for future expansion.

- maxp_interaction:

  `[integer(1) | NULL: NULL]`: Highest interaction order to purify.
  Higher-order components are not computed, but still influence lower
  orders during purification, which are the same as when purifying all
  orders.
  [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  then returns zero for the higher orders and their sum as `remainder`.
  [predict()](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md)
  still returns the whole forest's prediction. `NULL` purifies all
  orders.

- mode:

  `[integer(1): 2]`: Purification algorithm: `2` is the fast exact
  KD-tree based algorithm, which keeps the components as weighted
  rectangles, `1` the original grid-based one, which stores every
  component as a dense grid over the split points and can run out of
  memory from `max_interaction = 3` on. Both give the same components up
  to rounding.

- nthreads:

  `[integer(1) | NULL: NULL]`: Number of threads. `NULL` uses the
  `nthreads` the forest was fitted with, capped at the available cores.

## Value

`purify()` returns `x` invisibly. `is_purified()` returns `TRUE` or
`FALSE`.

## Details

The components are identified by the training data (the marginal
identification of Hiabu, Meyer and Wright, 2023). With `f` the forest's
prediction and `X_1, ..., X_n` the training rows, the component of a set
`S` of predictors is

`m_S(x_S) = sum over the subsets V of S of (-1)^(|S| - |V|) * mean_i f(x_V, X_(i,-V))`,

where `X_(i,-V)` are training row `i`'s values of the predictors outside
`V`. The intercept, the component of the empty set, is the mean
prediction on the training rows. The averages are over the joint
distribution of the training rows, so with dependent predictors a
component need not have mean zero along each of its predictors.

The forest is modified in place: `x` and every copy of it are purified,
whether or not the result is assigned. If the forest is already purified
with all the orders asked for, `purify()` returns it unmodified,
whatever `mode` and `nthreads` are. A forest purified with a
`maxp_interaction` below the order asked for (all orders with the
default `NULL`) is purified again, which needs its training data: a
forest restored by
[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
without it can only be purified up to its earlier cap.

## References

Hiabu, M., Meyer, J. T., & Wright, M. N. (2023). Unifying local and
global model explanations by functional decomposition of low dimensional
structures. Proceedings of The 26th International Conference on
Artificial Intelligence and Statistics (AISTATS), PMLR 206, 7040–7060.
<https://proceedings.mlr.press/v206/hiabu23a.html>

Hiabu, M., Mammen, E., & Meyer, J. T. (2020). Random Planted Forest: a
directly interpretable tree ensemble.
[doi:10.48550/arXiv.2012.14563](https://doi.org/10.48550/arXiv.2012.14563)

## Examples

``` r
rpfit <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, n_families = 10, prune = FALSE)
is_purified(rpfit)
#> [1] FALSE
purify(rpfit)
is_purified(rpfit)
#> [1] TRUE
```
