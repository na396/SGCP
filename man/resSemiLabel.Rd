\name{resSemiLabel}
\docType{data}
\alias{resSemiLabel}
\title{An example of output of \code{semiLabeling} function in SGCP pipeline.}
\description{
    This is an example of \code{geneOntology} function output as the last step in SGCP pipeline. Firstly, the adjacency matrix is produced using \code{adjacencyMatrix} function in SGCP over cheng dataset. The matrix is then used in~\code{clustering} function to produce the clustering result \code{resClus}. resClus is then used in \code{geneOnology} to produce~\code{resInitialGO}. This result is fed to \code{semiLabeling} to produce~\code{resSemiLabel}.}

\format{An object of class \code{list} containing the semi-labeling information.}

\usage{data(resSemiLabel)}

\details{\code{resSemiLabel} is a list containin the following information. \code{cutoff}: a numeric in (0,1) that shows the base line for identifying remarkable genes.  \code{geneLabel}: a dataframe of geneIDs and its corresponding label, NA labels means that correpsonding genes are unremarkable.}

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

