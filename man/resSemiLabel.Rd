\name{resSemiLabel}
\docType{data}
\alias{resSemiLabel}
\title{An example of the output from \code{semiLabeling} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{semiLabeling} function, representing the semi-label step in the SGCP pipeline. Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}. The clustering result, \code{resClus}, is subsequently utilized in the \code{geneOntology} function to derive \code{resInitialGO}, which captures the initial gene ontology (GO) enrichment results. The \code{resInitialGO} output is then processed through the \code{semiLabeling} function to produce \code{resSemiLabel}, indicating the semi-labeled genes based on their clustering characteristics.}

\format{An object of class \code{list} containing the semi-labeling information.}

\usage{data(resSemiLabel)}

\details{\code{resSemiLabel} is a list containing the following information.
\itemize{
    \item{\code{cutoff}: a numeric in (0,1) that shows the base line for identifying remarkable genes.}
    \item{\code{geneLabel}: a dataframe of geneIDs and its corresponding label, NA labels means that correpsonding genes are unremarkable.}}
}

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{semiLabeling}}
}
\examples{
    library(SGCP)
    data(resSemiLabel)
    summary(resSemiLabel)

    # cutoff value
    head(resSemiLabel$cutoff)

    # gene semi-label
    head(resSemiLabel$geneLabel)

}

