# To assist metacell object building in compute the purity of metacells with respect to each covariate

To assist metacell object building in compute the purity of metacells
with respect to each covariate

## Usage

``` r
mcRigor_purity(
  sc_covariate,
  sc_membership,
  method = c("max_proportion", "entropy")[1]
)
```

## Arguments

- sc_covariate:

  vector of single cell covariates

- sc_membership:

  vector of assignment of single-cell data to metacells

- method:

  method to compute metacell purity. `"max_proportion"` if the purity is
  defined as a proportion of the most abundant covariate (cell type)
  within super-cell or `"entropy"` if the purity is defined as the
  Shanon entropy of the covariates metacell consists of.

## Value

a vector of metacell purity, which is defined as: - proportion of the
most abundant covariate within metacell for `method = "max_proportion"`
or - Shanon entropy for `method = "entropy"`. With 1 meaning that
metacell consists of single cells from one covariate (reference
assignment)
