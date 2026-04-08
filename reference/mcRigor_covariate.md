# To assist metacell object building in determining the metacell covariates (metadata)

To assist metacell object building in determining the metacell
covariates (metadata)

## Usage

``` r
mcRigor_covariate(
  sc_covariate,
  sc_membership,
  method = c("jaccard", "relative", "absolute")
)
```

## Arguments

- sc_covariate:

  a vector (or dataframe) of single cell covariates

- sc_membership:

  a vector (or dataframe) of metacell memberships of single cells

- method:

  method to define the most abundant cell covariate within metacells.
  Available: "jaccard" (default), "relative", "absolute".

  - jaccard - assign metacell to covariate with the maximum jaccard
    coefficient (recommended)

  - relative - assign metacell to covariate with the maximum relative
    abundance (normalized by cluster size), may result in assignment of
    metacells to poorly represented (small) covariate due to
    normalization

  - absolute - assign metacell to covariate with the maximum absolute
    abundance within metacell, may result in disappearance of poorly
    represented (small) clusters

## Value

a vector of assigned metacell covariates
