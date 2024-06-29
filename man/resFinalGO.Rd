\name{resFinalGO}
\docType{data}
\alias{resFinalGO}
\title{An example of the output from \code{geneOntololgy} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{geneOntology} function, representing the final step in the SGCP pipeline. Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}. The clustering result, \code{resClus}, is subsequently utilized in the \code{geneOntology} function to derive \code{resInitialGO}, which captures the initial gene ontology (GO) enrichment results. The \code{resInitialGO} output is then processed through the \code{semiLabeling} function to produce \code{resSemiLabel}, indicating the semi-labeled genes based on their clustering characteristics. This semi-labeled information is further employed in the \code{semiSupervised} function, yielding \code{resSemiSupervised}, which includes the final supervised classification outcomes for the unremarkable genes. Finally, the results from \code{resSemiSupervised} are fed into the \code{geneOntology} function once more to generate \code{resFinalGO}, which represents the final GO enrichment analysis}
\format{An object of class \code{list} containing the gene ontology information for final gene ontology.}

\usage{data(resFinalGO)}

\details{code{resFinalGO} is a list containing the following information:
    \itemize{
    \item{\code{GOresults}: A dataframe of significant gene ontology terms and their corresponding test statistics.}
    \item{\code{FinalGOTermGenes}: A list of genes belonging to significant gene ontology terms per cluster.. }}
    }

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{geneOntology}}
}
\examples{
library(SGCP)
data(resFinalGO)
summary(resFinalGO)

# dataframe of significant gene ontology terms
head(resFinalGO$GOresults)

# a list of genes belong to significant gene ontology term for cluster 1
head(resFinalGO$FinalGOTermGenes$Cluster1_GOTermGenes)

# a list of genes belong to significant gene ontology term for cluster 2
head(resFinalGO$FinalGOTermGenes$Cluster2_GOTermGenes)
}

