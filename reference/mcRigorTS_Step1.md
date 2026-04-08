# Step 1 of mcRigor two-step: Identifying the single cells that constitute dubious metacells and need to be re-partitioned.

Step 1 of mcRigor two-step: Identifying the single cells that constitute
dubious metacells and need to be re-partitioned.

## Usage

``` r
mcRigorTS_Step1(
  obj_singlecell,
  sc_membership,
  TabMC = NULL,
  method = c("seacells", "mc1", "mc2", "supercell", "metaq")
)
```

## Arguments

- obj_singlecell:

  A Seurat object of the original single cells.

- sc_membership:

  A named vector (or dataframe) of the metacell membership of single
  cells.

- TabMC:

  A dataframe containing the permutation results and mcRigor results in
  dubious metacell detection. Saved in the previous steps through
  function mcRigor_DETECT.

- method:

  The metacell partitioning method used to build metacells. Available:
  "seacells" (default, SEACells), "mc1" (MetaCell), "mc2" (MetaCell2),
  "supercell" (SuperCell), "metaq" (MetaQ).

## Value

A list containing the following fields:

- obj_sc_dub:

  A Seurat object of single cells that are marked as dubious and will be
  re-partitioned.

- target_res:

  A list containing the permutation results, derived thresholds, and
  mcRigor results in dubious metacell detection, given by function
  mcRigor_threshold.

- obj_metacell_step1:

  The metacell object built under the membership sc_membership, with
  more dubious metacells marked.
