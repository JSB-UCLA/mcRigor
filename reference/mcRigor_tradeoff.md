# A building block of the main functions. To evaluate each metacell partition and optimize metacell partitioning based on the output permutation results (TabMC) and thresholds (threshold)

A building block of the main functions. To evaluate each metacell
partition and optimize metacell partitioning based on the output
permutation results (TabMC) and thresholds (threshold)

## Usage

``` r
mcRigor_tradeoff(
  TabMC,
  threshold = NULL,
  test_cutoff = 0.01,
  D_bw = 10,
  optim_method = c("tradeoff", "dub_rate_large", "dub_rate_small"),
  thre_smooth = T,
  thre_bw = 1/6,
  dub_rate = 0.1,
  weight = 0.5,
  draw = T
)
```

## Arguments

- TabMC:

  A dataframe containing the permutation results. Saved in the previous
  steps

- threshold:

  A dataframe containing the dubious metacell detection thresholds given
  by mcRigor_threshold

- test_cutoff:

  The test size for dubious metacell detection testing

- D_bw:

  A boolean indicating whether to smooth the dubious rate with respect
  to metacell size

- optim_method:

  The method used for granularity level optimization. Default is trading
  off between sparsity and dubious rate

- thre_smooth:

  A boolean indicating whether to smooth the threshold function

- thre_bw:

  If thre_smooth is True, thre_bw specifies the bandwidth for smoothing.

- dub_rate:

  If tradeoff is not used for optimization, what is highest acceptable
  dubious rate

- weight:

  The weight for dubious rate in the tradeoff.

- draw:

  A boolean indicating whether to visualize the mcRigor results

## Value

A list containing the following fields:

- optimized:

  The optimization results, containing the optimal gamma and its
  corresponding Sore

- scores:

  A data frame containing the evaluation scores for each gamma

- optim_plot:

  The line plot to visualize the tradeoff for hyperparameter
  opimization.
