# A building block of the main functions. To derive the thresholds for detecting dubious metacells based on the output permutation results (TabMC)

A building block of the main functions. To derive the thresholds for
detecting dubious metacells based on the output permutation results
(TabMC)

## Usage

``` r
mcRigor_threshold(
  TabMC,
  test_cutoff = 0.01,
  thre_smooth = T,
  thre_bw = 1/6,
  draw = F,
  palpha = 1,
  org_color = c("red", "orange", "yellow", "lightblue"),
  null_color = "#666666",
  pur_metric = NULL
)
```

## Arguments

- TabMC:

  A dataframe containing the permutation results. Saved in the previous
  steps

- test_cutoff:

  The test size for dubious metacell detection testing

- thre_smooth:

  A boolean indicating whether to smooth the threshold function

- thre_bw:

  If thre_smooth is True, what is the bandwidth for smoothing

- draw:

  A boolean indicating whether to visualize the mcRigor results

- palpha:

  Point alpha value for transparency in drawing.

- org_color:

  The colors indicating metacell purities or other interested factors

- null_color:

  The color for the null.

- pur_metric:

  Name of the covariate that we want to compute purity on. Can be NULL
  or a metadata variable name, ex. cell type.

## Value

A list containing the following fields:

- threshold:

  The thresholds for dubious metacell detection

- TabMC:

  A dataframe containing the permutation results and the testing results
  given by mcRigor

- test_plot:

  The scatter plots demonstrating the mcDiv values and the obtained
  thresholds for dubious metacell detection

- purity_plot:

  A violin plot showing the purity distribution of the pur_metric
  covariate in dubious metacells and trustworthy metacells
