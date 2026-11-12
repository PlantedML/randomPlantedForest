# Articles

### Get started

- [Getting Started with
  randomPlantedForest](https://plantedml.com/randomPlantedForest/articles/randomPlantedForest.md):

  A first model, its predictions and components, classification, and the
  settings you may want.

### Understanding the model

What a fitted model is made of, and how the forest behind it is grown.

- [Reading the
  model](https://plantedml.com/randomPlantedForest/articles/reading-the-model.md):

  The components of a fitted model, what the printed summary says line
  by line, importance, plots with glex, and how a component of three
  predictors is shown.

- [How the forest is
  grown](https://plantedml.com/randomPlantedForest/articles/growing.md):

  Tree families, max_interaction, how a split grows a family, the split
  structures, candidates and aging, shown on one small family and
  checked on data where the truth is known.

### Fitting

The choices
[`rpf()`](https://plantedml.com/randomPlantedForest/reference/rpf.md)
makes for you, how it makes them, and how to change them.

- [The number of splits and further
  tuning](https://plantedml.com/randomPlantedForest/articles/number-of-splits-and-tuning.md):

  How rpf() chooses the number of splits out of bag while the forest
  grows, when it stops, how to fix the number yourself, and what
  rpf(tune = ) searches.

- [Pruning to a sparse
  model](https://plantedml.com/randomPlantedForest/articles/pruning.md):

  Why rpf() prunes the forest’s components, what it leaves out, the
  rules and the threshold, choosing them out of bag or by
  cross-validation, the display of three-way components, and prune =
  FALSE.

- [Classification](https://plantedml.com/randomPlantedForest/articles/classification.md):

  Binary and multiclass outcomes, probabilities and classes, the
  components on the scale of the scores, and the losses.

### In practice

Large data, saving models, and other frameworks.

- [Speed, memory and saving
  models](https://plantedml.com/randomPlantedForest/articles/speed-and-memory.md):

  Threads, what the time and memory grow with, saving and restoring
  models, and using rpf() with tidymodels and mlr3.
