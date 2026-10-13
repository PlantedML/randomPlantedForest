# Plot a tree family

Illustrates how a Random Planted Forest predicts by drawing the trees of
one tree family. Each tree belongs to a set of variables, and its leaves
are boxes in those variables with a value each. The prediction for a
point is the root value plus the values of all leaves containing it,
averaged over tree families.

## Usage

``` r
# S3 method for class 'rpf'
plot(
  x,
  family = 1L,
  ...,
  type = c("tree", "boxes"),
  max_interaction = 2L,
  trees = NULL,
  class = NULL,
  factor_levels = TRUE
)
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

- type:

  `[character(1): "tree"]`: `"tree"` for a split diagram, `"boxes"` for
  the leaves in variable space.

- max_interaction:

  `[integer(1): 2]`: Draw trees on up to this many variables, at most 2
  for `type = "boxes"`. The remaining trees are summarized.

- trees:

  `[character | NULL: NULL]`: Names of the trees to draw, as in the
  `tree` column of
  [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md),
  e.g. `c("x1", "x1:x2")`. `NULL` draws all trees up to
  `max_interaction`.

- class:

  `[character(1) | NULL: NULL]`: For multiclass classification, the
  outcome level whose values to plot. `NULL` uses the first level.

- factor_levels:

  `[logical(1): TRUE]`: Label factors with their levels instead of
  integer codes. Levels are ordered by their association with the
  outcome, as used in the fit.

## Value

A `ggplot` for `type = "tree"`, a `patchwork` for `type = "boxes"`.

## Details

- `type = "tree"` draws a split diagram in the style of Figure 1b of
  Hiabu et al. (2023), with the trees in rows by number of variables,
  each leaf filled by its value. The forest stores only its final
  leaves, not the order of splits or which leaf a tree grew from, so the
  splits within a tree are reconstructed: the diagram shows a split
  structure consistent with the leaves, and dashed lines link trees to
  the trees they could have grown from. Leaves that overlap other leaves
  of the same tree hang from dotted lines. Readable for small fits, e.g.
  `ntrees = 1` and `splits = 10`.

- `type = "boxes"` draws one panel per tree in variable space: leaves of
  trees on one variable as bars from zero to their value, leaves of
  trees on two variables as rectangles filled by their value.
  Overlapping leaves add up, which the translucent drawing shows.

Trees on more than `max_interaction` variables, and trees not in
`trees`, are not drawn but summarized as a remainder, similar to the
remainder term of
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md).

After
[`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md),
the trees are the components of the functional decomposition and only
`type = "boxes"` applies. To plot components averaged over all tree
families, use
[`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md),
e.g. with the plotting functions in `glex`.

## References

Hiabu, M., Mammen, E., Meyer, J. T. (2023). Random Planted Forest: a
directly interpretable tree ensemble. arXiv:2012.14563.

## Examples

``` r
rpfit <- rpf(mpg ~ wt + hp + cyl, data = mtcars, ntrees = 1, splits = 10)
plot(rpfit)

plot(rpfit, type = "boxes")

purify(rpfit)
plot(rpfit, type = "boxes")
```
