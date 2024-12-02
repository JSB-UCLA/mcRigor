% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcRigor_function.R
\name{mcRigor_buildmc}
\alias{mcRigor_buildmc}
\title{To build an metacell object from the given metacell partitioning}
\usage{
mcRigor_buildmc(
  obj_singlecell,
  sc_membership = NULL,
  assay_type = c("RNA", "ATAC"),
  doNorm = T,
  aggregate_method = c("mean", "sum", "geom"),
  doAssign = T,
  fields = NULL,
  covariate_method = "absolute",
  purity_method = "max_proportion",
  add_testres = F,
  test_stats = NULL,
  Thre = NULL,
  test_cutoff = 0.01,
  prePro = F,
  feature_use = 2000,
  gene_filter = 0.1,
  cor_method = c("pearson", "spearman")
)
}
\arguments{
\item{obj_singlecell}{the Seurat object of single cells}

\item{sc_membership}{A named vector (or dataframe) of the metacell membership of single cells}

\item{assay_type}{The type of data assay yuo are using, depending on which different normalization would be used.}

\item{doNorm}{A bool indicating whether to perform normalization for the metacell object (obj_metacell)}

\item{aggregate_method}{The method to aggregate single cell profiles into metacell profiles}

\item{doAssign}{A bool indicating whether to assign covariates to metacells or not}

\item{fields}{A vector of covariate names to assign}

\item{covariate_method}{Method to define the most abundant cell covariate within metacells. Available: "jaccard", "relative", "absolute" (default).
\itemize{
  \item jaccard - assign metacell to covariate with the maximum jaccard coefficient (recommended)
  \item relative - assign metacell to covariate with the maximum relative abundance (normalized by cluster size), may result in assignment of metacells to poorly represented (small) covariate due to normalization
  \item absolute - assign metacell to covariate with the maximum absolute abundance within metacell, may result in disappearance of poorly represented (small) clusters
}}

\item{purity_method}{method to compute metacell purity.
\code{"max_proportion"} if the purity is defined as a proportion of the most abundant covariate (cell type) within super-cell or
\code{"entropy"} if the purity is defined as the Shanon entropy of the covariates metacell consists of.}

\item{add_testres}{A bool indicating whether to add the mcRigor results (dubious or trustworthy) as part of obj_metacell's metadata}

\item{test_stats}{If add_testres = True, this argument is needed. Usually should be TabMC, an output from previous steps.}

\item{Thre}{The threshold for dubious metacell detection. If not inputed, it will be computed based on test_stats.}

\item{test_cutoff}{The test size for dubious metacell detection testing}

\item{prePro}{A boolean indicating whether to normalize obj_singlecell for preprocessing.}

\item{feature_use}{The number of genes to use in metacell testing.}

\item{gene_filter}{A proportion. Genes expressed lower than this proportion will be filtered out.}

\item{cor_method}{The method for gene correlation calculation description}
}
\value{
a Seurat object of the metacells
}
\description{
To build an metacell object from the given metacell partitioning (sc_membership)
}
