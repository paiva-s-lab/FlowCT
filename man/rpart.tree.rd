% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop.cutoff.R
\name{rpart.tree}
\alias{rpart.tree}
\title{Recursive partitioning tree}
\usage{
rpart.tree(fcs.SCE, cell.clusters, variables, value = "percentage", time.var, event.var, xerror, return.data = F)
}
\arguments{
\item{fcs.SCE}{A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.}

\item{cell.clusters}{Name of column containing clusters identified through \code{\link[FlowCT.v2:clustering.flow]{FlowCT.v2::clustering.flow()}}.}

\item{variables}{Vector with variables for calculating the cutoff. If nothing is detailed (\code{NULL}, default), all immune populations from \code{cell.clusters} will be considered.}

\item{value}{String specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").}

\item{time.var}{Survival time variable.}

\item{event.var}{Variable with event censoring.}

\item{xerror}{Numeric value with cross-validation error for pruning the tree. If missing, the minimal error will be automatically selected.}

\item{return.data}{Original, pruned and used metadata should be returned?. Default = \code{FALSE}.}

}
\description{
This function performs a recursive partitioning tree using the percentage (or raw counts) of identified cell populations. Depending if \code{time.var} is present or absent, final tree were a cassification or regression tree, respectively.
}
\examples{
\dontrun{
tr <- rpart.tree(fcs.SCE = fcs, cell.clusters = "clusters_named", time.var = "PFS", event.var = "PFS_c", return.data = t)
}
}
\keyword{recursive}
\keyword{partitioning}
\keyword{CART}
\keyword{classification}
\keyword{regression}
\keyword{tree}