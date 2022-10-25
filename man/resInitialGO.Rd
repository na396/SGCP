\name{resInitialGO}
\docType{data}
\alias{resInitialGO}
\title{An example of output of \code{geneOntololgy} function in SGCP pipeline.}
\description{
    This is an example of \code{geneOntology} function output as the third step in SGCP pipeline. Firstly, the adjacency matrix is produced using \code{adjacencyMatrix} function in SGCP over cheng dataset. The matrix is then used in~\code{clustering} function to produce the clustering result \code{resClus}. resClus is then used in \code{geneOnology} to produce~\code{resInitialGO}. }
\format{An object of class \code{list} containing the gene ontology information for final gene ontology.}

\usage{data(resInitialGO)}

\details{\code{resInitialGO} is a list containin the following information. \code{GOresults}: a dataframe of significant gene ontology terms and their corresponding test statistics information.  \code{FinalGOTermGenes}: a list of the genes belong to ignificant gene ontology terms per cluster. }

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

