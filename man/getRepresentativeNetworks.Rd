% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRepresentativeNetworks.R
\name{getRepresentativeNetworks}
\alias{getRepresentativeNetworks}
\title{Representative Networks}
\usage{
getRepresentativeNetworks(sampleIDs, dim_subset, SubpopSizeFilter,
  num_networkEdge)
}
\arguments{
\item{sampleIDs}{sampleID vector}

\item{dim_subset}{a string vector of string names to subset the data columns for PAC; set to NULL to use all columns}

\item{SubpopSizeFilter}{the cutoff for small subpopulations. Smaller subpopulations have unstable covariance structure, so no network structure is calculated}

\item{num_networkEdge}{the number of edges to draw for each subpopulation mutual information network}
}
\description{
Outputs representative networks for clades/subpopulations larger than a size filter (very small subpopulations are not considered in downstream analyses)
}

