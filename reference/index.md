# Package index

## All functions

- [`mcRigor-package`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor-package.md)
  [`mcRigor`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor-package.md)
  : mcRigor: a statistical tool to enhance the rigor of metacell
  partitioning
- [`mcRigorTS_Step1()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigorTS_Step1.md)
  : Step 1 of mcRigor two-step: Identifying the single cells that
  constitute dubious metacells and need to be re-partitioned.
- [`mcRigorTS_Step2()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigorTS_Step2.md)
  : Step 2 of mcRigor two-step: Re-partition single cells that
  previously formed dubious metacells into smaller metacells.
- [`mcRigor_DETECT()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_DETECT.md)
  : Main functionality 1: To detect dubious metacells for a given
  metacell partition
- [`mcRigor_OPTIMIZE()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_OPTIMIZE.md)
  : Main functionality 2: To select the optimal hyperparameter for
  metacell partitioning
- [`mcRigor_buildmc()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_buildmc.md)
  : To build an metacell object from the given metacell partitioning
- [`mcRigor_covariate()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_covariate.md)
  : To assist metacell object building in determining the metacell
  covariates (metadata)
- [`mcRigor_projection()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_projection.md)
  : To visualize a given metacell partitioning by projecting metacells
  onto the single cell embedding space.
- [`mcRigor_purity()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_purity.md)
  : To assist metacell object building in compute the purity of
  metacells with respect to each covariate
- [`mcRigor_threshold()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_threshold.md)
  : A building block of the main functions. To derive the thresholds for
  detecting dubious metacells based on the output permutation results
  (TabMC)
- [`mcRigor_tradeoff()`](https://jsb-ucla.github.io/mcRigor/reference/mcRigor_tradeoff.md)
  : A building block of the main functions. To evaluate each metacell
  partition and optimize metacell partitioning based on the output
  permutation results (TabMC) and thresholds (threshold)
