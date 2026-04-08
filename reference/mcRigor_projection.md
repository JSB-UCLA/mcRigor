# To visualize a given metacell partitioning by projecting metacells onto the single cell embedding space.

To visualize a given metacell partitioning by projecting metacells onto
the single cell embedding space.

## Usage

``` r
mcRigor_projection(
  obj_singlecell,
  sc_membership = NULL,
  sc.reduction = c("umap", "tsne", "pca"),
  dims = c(1, 2),
  metric = "size",
  add_testres = F,
  test_stats = NULL,
  Thre = NULL,
  test_cutoff = 0.01,
  color_field = NULL,
  cpalette = c("#13678A", "#45C4B0", "#9AEBA3", "#BF847E", "#F2C12E", "#7E57C2",
    "#FACFCE", "#9E9D24", "#DAFDBA", "#86ABD4", "#42A5F5", "#546E7A", "#D4E157",
    "#76FF03", "#6D4C41", "#004D40", "#AB47BC", "#D81B60"),
  sc.alpha = 0.5,
  mc.alpha = 1,
  pt_size = 1,
  axis_lab = F,
  max_mcsize = 1000,
  continuous_metric = T,
  dub_mc.label = F,
  dub_mc_test.label = F,
  label_text = F
)
```

## Arguments

- obj_singlecell:

  the Seurat object of single cells

- sc_membership:

  A named vector (or dataframe) of the metacell membership of single
  cells

- sc.reduction:

  The single cell reduction method to produce single cell embeddings

- dims:

  The dimensions to use in the single cell embeddings

- metric:

  The variable that determines the sizes of dots representing the
  metacells

- add_testres:

  A bool indicating whether to add the mcRigor results (dubious or
  trustworthy) as part of obj_metacell's metadata

- test_stats:

  If add_testres = True, this argument is needed. Usually should be
  TabMC, an output from previous steps.

- Thre:

  The threshold for dubious metacell detection. If not inputed, it will
  be computed based on test_stats.

- test_cutoff:

  The test size for dubious metacell detection testing

- color_field:

  The variable based on which to color the cells.

- cpalette:

  The color palette for visualization

- sc.alpha:

  The level of single cell point transparency for visualization

- mc.alpha:

  The level of metacell circle transparency for visualization

- pt_size:

  The single cell point size for visualization

- axis_lab:

  A bool indicating whether to add axis labels in visualization

- max_mcsize:

  The maximum metacell size used for visualization scaling

- continuous_metric:

  A bool indicating whether the field "metric" contains continuous
  values

- dub_mc.label:

  A bool indicating whether to mark the ground truth dubious metacells
  using red color or not

- dub_mc_test.label:

  A bool indicating whether to mark the detected dubious metacells using
  red color or not

- label_text:

  A bool indicating whether to label the names of metacells or not

## Value

a scatter plot of metacells projected to the single cell 2D embedding
space.
