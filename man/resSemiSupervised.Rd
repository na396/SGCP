\name{resSemiSupervised}
\docType{data}
\alias{resSemiSupervised}
\title{An example of output of \code{semiSupervised} function in SGCP pipeline.}
\description{
    This is an example of \code{geneOntology} function output as the last step in SGCP pipeline. Firstly, the adjacency matrix is produced using \code{adjacencyMatrix} function in SGCP over cheng dataset. The matrix is then used in~\code{clustering} function to produce the clustering result \code{resClus}. resClus is then used in \code{geneOnology} to produce~\code{resInitialGO}. This result is fed to \code{semiLabeling} to produce~\code{resSemiLabel}. This is result is used as input to~\code{semiSupervised} function to produce \code{resSemiSupervised}.}

\format{An object of class \code{list} containing the semi-supervised information.}

\usage{data(resSemiSupervised)}

\details{\code{resSemiSupervised} is a list containin the following information. \code{semiSupervised}: an object of \code{caret} for the training model.  \code{prediction}: A vector of predicted labels for unremakable genes. \code{FinalLabeling}: a dataframe gene semil-label and final predicted labels.}

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{semiLabeling}}
}
\examples{
    library(SGCP)
    data(resSemiSupervised)

    # supervised model information
    summary(resSemiSupervised$semiSupervised)

    # predicted label for unremarkable genes
    head(resSemiSupervised$prediction)

    # gene semi and final labeling
    head(resSemiSupervised$FinalLabeling)

}
