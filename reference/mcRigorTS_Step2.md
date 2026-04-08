# Step 2 of mcRigor two-step: Re-partition single cells that previously formed dubious metacells into smaller metacells.

Step 2 of mcRigor two-step: Re-partition single cells that previously
formed dubious metacells into smaller metacells.

## Usage

``` r
mcRigorTS_Step2(
  step1_res,
  obj_singlecell,
  cell_membership_twostep,
  twostep_gamma = NULL,
  method = c("seacells", "mc1", "mc2", "supercell", "metaq"),
  custom_thre = F,
  color_field = NULL
)
```

## Arguments

- step1_res:

  A list of results obtained form Step 1 of mcRigor two-step, given by
  function mcRigorTS_Step1.

- obj_singlecell:

  A Seurat object of the original single cells.

- cell_membership_twostep:

  A dataframe, each column of which contains the metacell membership of
  single cells that need to be re-partitioned under a given gamma
  (granularity level). The column names should be the corresponding
  gamma's (in the character type). The row names should be the single
  cell id's.

- twostep_gamma:

  The granularity level value for the second step. Default is NULL.

- method:

  The metacell partitioning method used to build metacells. Available:
  "seacells" (default, SEACells), "mc1" (MetaCell), "mc2" (MetaCell2),
  "supercell" (SuperCell), "metaq" (MetaQ).

- custom_thre:

  A boolean indicating whether to re-compute the threshold values.

- color_field:

  The variable based on which to color the cells.

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
