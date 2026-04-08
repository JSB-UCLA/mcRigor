# Main functionality 2: To select the optimal hyperparameter for metacell partitioning

Main functionality 2: To select the optimal hyperparameter for metacell
partitioning

## Usage

``` r
mcRigor_OPTIMIZE(
  obj_singlecell,
  cell_membership = NULL,
  assay_type = c("RNA", "ATAC"),
  Gammas = NULL,
  aggregate_method = c("mean", "sum", "geom"),
  output_file = NULL,
  Nrep = 1,
  gene_filter = 0.1,
  feature_use = 2000,
  cor_method = c("pearson", "spearman"),
  prePro = T,
  test_cutoff = 0.01,
  thre_smooth = T,
  thre_bw = 1/6,
  D_bw = 10,
  optim_method = c("tradeoff", "dub_rate_large", "dub_rate_small"),
  weight = 0.5,
  dub_rate = 0.1,
  draw = T,
  pur_metric = NULL,
  check_purity = T,
  fields = NULL,
  step_save = T
)
```

## Arguments

- obj_singlecell:

  A Seurat object of the original single cells.

- cell_membership:

  A dataframe, each column of which contains the metacell membership of
  single cells under a given gamma. The column names should be the
  corresponding gamma's. The row names should be the single cell id's.

- assay_type:

  The type of data assay yuo are using, depending on which different
  normalization would be used.

- Gammas:

  The candidate pool of granularity levels to consider in optimization

- aggregate_method:

  The method to aggregate single cell profiles into metacell profiles

- output_file:

  The directory for output saving

- Nrep:

  The number of permutation repetitions we use when deriving the null.

- gene_filter:

  A proportion. Genes expressed lower than this proportion will be
  filtered out.

- feature_use:

  The number of genes to use in metacell testing.

- cor_method:

  The method for gene correlation calculation description

- prePro:

  A boolean indicating whether to normalize obj_singlecell for
  preprocessing.

- test_cutoff:

  The test size for dubious metacell detection testing

- thre_smooth:

  A boolean indicating whether to smooth the threshold function

- thre_bw:

  If thre_smooth is True, thre_bw specifies the bandwidth for smoothing.

- D_bw:

  A boolean indicating whether to smooth the dubious rate with respect
  to metacell size

- optim_method:

  The method used for granularity level optimization. Default is trading
  off between sparsity and dubious rate

- weight:

  The weight for dubious rate in the tradeoff.

- dub_rate:

  If tradeoff is not used for optimization, what is highest acceptable
  dubious rate

- draw:

  A boolean indicating whether to visualize the mcRigor results

- pur_metric:

  Can be NULL or a metadata variable name, ex. cell type.

- check_purity:

  A boolean indicating whether to calculate the metacell purity of
  specific fields or not.

- fields:

  A vector of the fields of interest, ex. celltype. It should be a
  subset of obj_singlecell's meta.data.

- step_save:

  A boolean indicating whether to save the outputs step by step

## Value

A list containing the following fields:

- best_granularity_level:

  The optimal granularity level selected by mcRigor

- best_Score:

  The evaluation score for the metacell partition given by the optimal
  granularity level selected by mcRigor

- opt_metacell:

  The metacell object build under the optimal granularity level

- scores:

  A data frame containing the evaluation scores for each gamma

- optim_plot:

  The line plot to visualize the tradeoff for hyperparameter
  opimization.

- thre:

  The thresholds for dubious metacell detection

- TabMC:

  A dataframe containing the permutation results, elements to calculate
  the test statistics mcDiv and mcDiv null
