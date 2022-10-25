\name{resClus}
\docType{data}
\alias{resClus}
\title{An example of output of \code{clustering} function in SGCP pipeline.}
\description{
This is an example of \code{clustering} function output. Firstly, the adjacency matrix is produced using \code{adjacencyMatrix} function in SGCP over cheng dataset. The matrix is then used in~\code{clustering} function to produce the clustering result.}
\format{An object of clas \code{list} containing the clustering information.}

\usage{data(resClus)}

\details{\code{resClus} is a list containing the following clustering information. \code{dropped.indices}: a vector of dropped gene indices, \code{geneID}:a vector of geneIDs, \code{method}:indicates the selected method for number of cluster,\code{k}:selected number of clusters, \code{Y}: transformed matrix with 2*k columns, \code{X}:eigenvalues correspond to 2*k columns in Y, \code{cluster}:an object of class \code{kmeans}, \code{clusterLabels}:a vector containing the cluster label per gene, there is a 1-to-1 correspondence between geneID and clusterLabes, \code{conductance}:a list containing mean and median, and individual cluster conductance index for clusters per method. Index in `clusterConductance` field denotes the \code{cvGOdf}:a dataframe used for gene ontology validation. For each method, it returns the gene ontology enrichment result on the cluster with minimum conductance index, \code{cv}: an string indicates the validation method for number of cluster, "cvGO" means gene ontology validation used, \code{clusterNumberPlot}:an object of class \code{ggplot2} for  relativeGap"", "secondOrderGap", and "additiveGap", \code{silhouette}:a dataframe that indicates the silhouette for genes,\code{original}:a list with matrix transformation and corresponding eigenvalues and n_egvec, where n_egvec top columns of transformation is kept.}

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{adjacencyMatrix}}
    \code{\link{clustering}}
}
\examples{
library(SGCP)
data(resClus)
summary(resClus)
resClus
}
