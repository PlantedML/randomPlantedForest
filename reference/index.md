# Package index

## Fitting and prediction

Fit a model, predict with it, and print it.

- [`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
  : Random Planted Forest
- [`predict(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/reference/predict.rpf_pruned.md)
  : Predictions of a pruned Random Planted Forest
- [`predict(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/reference/predict.rpf.md)
  : Predict with a Random Planted Forest
- [`print(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/reference/print.rpf_pruned.md)
  : Print a pruned Random Planted Forest
- [`print(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/reference/print.rpf.md)
  : Print a Random Planted Forest

## Reading the model

The components that add up to the prediction, and the purification that
makes them unique.

- [`predict_components()`](https://plantedml.com/randomPlantedForest/reference/predict_components.md)
  : The components of a Random Planted Forest's prediction
- [`purify()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  [`is_purified()`](https://plantedml.com/randomPlantedForest/reference/purify.md)
  : Purify a Random Planted Forest

## Pruning

Leave out the small components, by a rule chosen out of bag or by
cross-validation.

- [`prune_components()`](https://plantedml.com/randomPlantedForest/reference/prune_components.md)
  : Prune the components of a Random Planted Forest

## The number of splits and tuning

The out-of-bag error along the splits, and the search over the
greediness and the leaf size.

- [`oob_path_checkpoints()`](https://plantedml.com/randomPlantedForest/reference/oob_path_checkpoints.md)
  : Checkpoints of the out-of-bag path
- [`tune_control()`](https://plantedml.com/randomPlantedForest/reference/tune_control.md)
  : Settings of the out-of-bag tuning
- [`greediness_ladder()`](https://plantedml.com/randomPlantedForest/reference/greediness_ladder.md)
  : The greediness ladder of the out-of-bag tuning

## Looking inside a tree family

Draw the trees of one family, or get its leaves as boxes.

- [`plot(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/reference/plot.rpf.md)
  : Plot a tree family
- [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/reference/rpf_boxes.md)
  : Extract the boxes of a tree family

## Saving and restoring

Models hold a C++ object, which needs marshalling to survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) or a worker process.

- [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/reference/rpf_marshal.md)
  : Serialize and restore Random Planted Forests
- [`rpf_is_valid()`](https://plantedml.com/randomPlantedForest/reference/rpf_is_valid.md)
  : Check whether a model's C++ forest is still alive
- [`bundle(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/reference/bundle.rpf.md)
  [`bundle(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/reference/bundle.rpf.md)
  : Bundle a Random Planted Forest

## Utilities

Compact printing of the flattened forest (`export_forest = TRUE`).

- [`print(`*`<rpf_forest>`*`)`](https://plantedml.com/randomPlantedForest/reference/print.rpf_forest.md)
  [`str(`*`<rpf_forest>`*`)`](https://plantedml.com/randomPlantedForest/reference/print.rpf_forest.md)
  : Compact printing of forest structures
