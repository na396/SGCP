\name{resInitialGO}
\docType{data}
\alias{resInitialGO}
\title{An example of the output from the \code{geneOntololgy} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{geneOntology} function, representing the third step in the SGCP pipeline. Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}. The clustering result, \code{resClus}, is subsequently utilized in the \code{geneOntology} function to derive \code{resInitialGO}, which captures the initial gene ontology (GO) enrichment results.}


\format{An object of class \code{list} containing the gene ontology information for initial gene ontology.}

\usage{data(resInitialGO)}

\details{ \code{resInitialGO} is a list containing the following information.
    \itemize{
    \item{\code{GOresults}: a dataframe of significant gene ontology terms and their corresponding test statistics information.}
    \item{\code{FinalGOTermGenes}: a list of the genes belong to significant gene ontology terms per cluster. }}
    }


\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{geneOntology}}
}
\examples{
    library(SGCP)
    data(resInitialGO)
    summary(resInitialGO)

    # dataframe of significant gene ontology terms
    head(resInitialGO$GOresults)

    # a list of genes belong to significant gene ontology term for cluster 1
    head(resInitialGO$FinalGOTermGenes$Cluster1_GOTermGenes)

    # a list of genes belong to significant gene ontology term for cluster 2
    head(resInitialGO$FinalGOTermGenes$Cluster2_GOTermGenes)
}

