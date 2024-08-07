% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_clustering.R
\name{Clustering}
\alias{Clustering}
\title{Clustering function}
\usage{
Clustering(
  dataset,
  nb_clust_possible,
  n_iter = 5,
  Dmin = 0.28,
  SampleName = NULL
)
}
\arguments{
\item{dataset}{dataset with Contrast and SigStren for each individuals (as SampleName) and each markers (as MarkerName)}

\item{nb_clust_possible}{number of cluster possible (ploidy+1)}

\item{n_iter}{number of iterations to perform for clustering}

\item{Dmin}{minimal distance between two clusters}

\item{SampleName}{vector with all SampleName (important when missing genotype)}
}
\value{
list of results of clustering
}
\description{
Clustering function to run clustering with no parallelization process nor auto save
}
\examples{
data(GenoTriplo_to_clust)
ploidy=3
res = Clustering(dataset=GenoTriplo_to_clust,
                 nb_clust_possible=ploidy+1,n_iter=5)

}
