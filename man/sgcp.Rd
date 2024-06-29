\name{sgcp}
\docType{data}
\alias{sgcp}
\title{An example of the output of \code{ezSGCP} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{ezSGCP} function, representing the entire SGCP pipeline.Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}. The clustering result, \code{resClus}, is subsequently utilized in the \code{geneOntology} function to derive \code{resInitialGO}, which captures the initial gene ontology (GO) enrichment results. The \code{resInitialGO} output is then processed through the \code{semiLabeling} function to produce \code{resSemiLabel}, indicating the semi-labeled genes based on their clustering characteristics. This semi-labeled information is further employed in the \code{semiSupervised} function, yielding \code{resSemiSupervised}, which includes the final supervised classification outcomes for the unremarkable genes. Finally, the results from \code{resSemiSupervised} are fed into the \code{geneOntology} function once more to generate \code{resFinalGO}, which represents the final GO enrichment analysis.}

\usage{data(sgcp)}
\format{An object of class \code{list} containing the \code{ezSGCP} function information.}


\details{
\code{sgcp} contains a list with the following fields:
\itemize{
  \code{clustering}: List of clustering
  \itemize{
  \item{\code{dropped.indices}:} Dropped gene indices.
    \item{\code{geneID}:} Vector of geneIDs.
    \item{\code{method}:} Selected method for determining the number of clusters.
    \item{\code{k}:} Selected number of clusters.
    \item{\code{Y}:} Transformed matrix with 2*k columns.
    \item{\code{X}:} Eigenvalues corresponding to the 2*k columns in Y.
    \item{\code{cluster}:} Object of class \code{kmeans}.
    \item{\code{clusterLabels}:} Vector containing cluster labels for each gene.
    \item{\code{conductance}:} List containing mean, median, and individual cluster conductance indices. Each method's \code{clusterConductance} field denotes the cluster label with its corresponding conductance index.
  \item{\code{cvGOdf}:} DataFrame used for gene ontology validation. For each method, it shows gene ontology enrichment on the cluster with the smallest conductance index.
  \item{\code{cv}:} String indicating the validation method for the number of clusters (\code{"cvGO"} for gene ontology validation).
  \item{\code{clusterNumberPlot}:} Object of class \code{ggplot2} for displaying relativeGap, secondOrderGap, and additiveGap.
  \item{\code{silhouette}:} DataFrame indicating silhouette values for genes.
  \item{\code{original}:} List with matrix transformation, eigenvalues, and \code{n_egvec}, retaining the top columns of transformation.

  }

  \code{initial.GO}: List of GO term analysis results for initial clusters
  \itemize{
    \item{\code{GOresults}: DataFrame summarizing GO term information.}
    \item{\code{FinalGOTermGenes}: List containing geneIDs of each GO term per cluster.}
  }

  \code{semiLabeling}: List of semi-labeling results
  \itemize{
    \item{\code{cutoff}}: Numeric indicating selected cutoff.
    \item{\code{geneLabel}}: DataFrame with geneID and corresponding cluster label (or NA if unremarkable).
  }

  \code{semiSupervised}: List of semi-supervised learning results
  \itemize{
  \item{\code{semiSupervised}: Object of classification result.}
  \item{\code{prediction}: Vector of predicted labels for unremarkable genes.}
  \item{\code{FinalLabeling}: DataFrame of geneID with corresponding semi-label and final label.}
  }

  \code{final.GO}}: List of GO term analysis results for final modules
  \itemize{
  \item{\code{GOresults}: DataFrame summarizing GO term information.}
  \item{\code{FinalGOTermGenes}: List containing geneIDs of each GO term per cluster.}
  }


}

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{ezSGCP}}
}
\examples{
library(SGCP)
data(sgcp)
summary(sgcp)

# clustering step
summary(sgcp$clustering)

# intial gene ontology step
summary(sgcp$initial.GO)

# semilabeling step
summary(sgcp$semiLabeling)

# semi-supervised step
summary(sgcp$semiSupervised)

# final gene ontology step
summary(sgcp$final.GO)

}
