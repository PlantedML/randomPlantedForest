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
purifying.

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
#> 2          wt     1   1  8.2166667    2.200    2.320       NA       NA
#> 3          wt     1   2  0.1649621    3.190    5.444       NA       NA
#> 4          wt     1   3  1.8166667    1.513    2.200       NA       NA
#> 5          wt     1   4 -1.7078125    2.320    3.190       NA       NA
#> 6          hp     1   1 23.2078125       NA       NA       97      105
```
