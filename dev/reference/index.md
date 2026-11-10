# Package index

## Model fitting and prediction

Modeling fundamentals for creation and prediction from Random Planted
Forest models.

- [`rpf()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf.md)
  : Random Planted Forest
- [`predict(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf.md)
  : Random Planted Forest Predictions
- [`oob_path_checkpoints()`](https://plantedml.com/randomPlantedForest/dev/reference/oob_path_checkpoints.md)
  : Checkpoints of the out-of-bag path

## Tuning

Tuning the greediness and the leaf size out of bag

- [`tune_control()`](https://plantedml.com/randomPlantedForest/dev/reference/tune_control.md)
  : Settings of the out-of-bag tuning
- [`greediness_ladder()`](https://plantedml.com/randomPlantedForest/dev/reference/greediness_ladder.md)
  : The greediness ladder of the out-of-bag tuning

## Functional decomposition

Enabling and utilizing functional decomposition

- [`purify()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  [`is_purified()`](https://plantedml.com/randomPlantedForest/dev/reference/purify.md)
  : Purify a Random Planted Forest
- [`predict_components()`](https://plantedml.com/randomPlantedForest/dev/reference/predict_components.md)
  : Extract predicted components from a Random Planted Forest

## Pruning

Dropping small components, chosen by cross-validation or out of bag

- [`prune_components()`](https://plantedml.com/randomPlantedForest/dev/reference/prune_components.md)
  : Prune the components of a Random Planted Forest
- [`predict(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/predict.rpf_pruned.md)
  : Predictions of a pruned Random Planted Forest

## Visualization

Inspecting the trees of a forest

- [`plot(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/plot.rpf.md)
  : Plot a tree family
- [`rpf_boxes()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_boxes.md)
  : Extract the boxes of a tree family

## Utility functions

- [`print(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf.md)
  : Print an rpf
- [`print(`*`<rpf_forest>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_forest.md)
  [`str(`*`<rpf_forest>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_forest.md)
  : Compact printing of forest structures
- [`print(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/print.rpf_pruned.md)
  : Print a pruned Random Planted Forest
- [`rpf_marshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  [`rpf_unmarshal()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_marshal.md)
  : Serialize and restore Random Planted Forests
- [`rpf_is_valid()`](https://plantedml.com/randomPlantedForest/dev/reference/rpf_is_valid.md)
  : Check whether an rpf object's C++ forest is still alive
- [`bundle(`*`<rpf>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/bundle.rpf.md)
  [`bundle(`*`<rpf_pruned>`*`)`](https://plantedml.com/randomPlantedForest/dev/reference/bundle.rpf.md)
  : Bundle an rpf model
