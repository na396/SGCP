\name{resClus}
\docType{data}
\alias{resClus}
\title{An example of the output from \code{clustering} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{clustering} function, representing the network clustering step in the SGCP pipeline. Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}}


\format{An object of clas \code{list} containing the clustering information.}

\usage{data(resClus)}

\details{\code{resClus} is a list containing the following clustering information:

\itemize{
        \item{\code{dropped.indices}: A vector of dropped gene indices.}
        \item{\code{geneID}: A vector of gene IDs.}
        \item{\code{method}: Indicates the selected method for determining the number of clusters.}
        \item{\code{k}: The selected number of clusters.}
        \item{\code{Y}: Transformed matrix with 2*k columns.}
        \item{\code{X}: Eigenvalues corresponding to the 2*k columns in Y.}
        \item{\code{cluster}: An object of class \code{kmeans}.}
        \item{\code{clusterLabels}: A vector containing the cluster label for each gene. There is a 1-to-1 correspondence between \code{geneID} and \code{clusterLabels}}.
        \item{\code{conductance}: A list containing the mean, median, and individual cluster conductance index for clusters per method. The index in the \code{clusterConductance} field denotes the method.}
        \item{\code{cvGOdf}: A dataframe used for gene ontology validation. For each method, it returns the gene ontology enrichment result on the cluster with the minimum conductance index.}
        \item{\code{cv}: A string indicating the validation method for the number of clusters; "cvGO" means gene ontology validation was used.}
        \item{\code{clusterNumberPlot}: An object of class \code{ggplot2} for \code{relativeGap}, \code{secondOrderGap}, and \code{additiveGap}.}
        \item{\code{silhouette}: A dataframe indicating the silhouette values for genes.}
        \item{\code{original}: A list with matrix transformation, corresponding eigenvalues, and \code{n_egvec}, where the top}
        \item{\code{n_egvec} columns of the transformation are kept.}

    }
}

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
