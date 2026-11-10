# Extract the boxes of a tree family

Every tree in a Random Planted Forest is a set of boxes (leaves) with a
value each, and the prediction is the sum of the values of all boxes
containing the point, averaged over tree families. `rpf_boxes()` returns
the boxes of one tree family as a data frame.

## Usage

``` r
rpf_boxes(x, family = 1L, ...)
```

## Arguments

- x:

  `[rpf]`: A fitted
  [`rpf`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  model.

- family:

  `[integer(1): 1]`: Which tree family to extract, between 1 and
  `ntrees`.

- ...:

  Reserved for future expansion.

## Value

A [base::data.frame](https://rdrr.io/r/base/data.frame.html) with one
row per box and columns

- `tree`: The tree, named by its sorted variables as in
  [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
  e.g. `"x1:x2"`, or `"(Intercept)"`.

- `order`: The number of variables of the tree, `0` for the intercept.

- `box`: The index of the box within its tree.

- `value`: The value of the box. For multiclass classification, one
  column `value_<class>` per outcome level instead.

- `<variable>_lower`, `<variable>_upper` for every predictor: The bounds
  of the box, `NA` for variables the tree does not use.

## Details

Before purification, these are the leaves as grown: boxes in one tree
may overlap, and the trees are not an identified decomposition. After
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
they are the cells of the purified grid, one non-overlapping set per
component of the functional decomposition. As
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
modifies the forest in place, extract the unpurified boxes before
purifying. After `purify(mode = 2)`, the default, the cells are computed
from the training data, so a forest restored by
[`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
needs to have been marshaled with `include_data = TRUE`.

Factor predictors are encoded as integer codes in the order of
`x$factor_levels`, with a box covering the codes `c` with
`lower <= c < upper`.

## See also

[`plot.rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/plot.rpf.md)
to plot the boxes.

## Examples

``` r
rpfit <- rpf(mpg ~ wt + hp, data = mtcars, ntrees = 5, splits = 10)
head(rpf_boxes(rpfit))
#>          tree order box      value wt_lower wt_upper hp_lower hp_upper
#> 1 (Intercept)     0   1  0.0000000       NA       NA       NA       NA
#> 2          hp     1   1 24.0800000       NA       NA       52   123.00
#> 3          hp     1   2 14.4859561       NA       NA      230   335.02
#> 4          hp     1   3 17.0200000       NA       NA      123   205.00
#> 5          hp     1   4  2.9575758       NA       NA       52    66.00
#> 6          hp     1   5 -0.3059561       NA       NA       66   335.02
```
