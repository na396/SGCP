\name{resSemiSupervised}
\docType{data}
\alias{resSemiSupervised}
\title{An example of the output from \code{semiSupervised} function in the SGCP pipeline}
\description{
This is an example of the output from the \code{semiSupervised} function, representing the semi-supervised step in the SGCP pipeline. Initially, the adjacency matrix is generated using the \code{adjacencyMatrix} function within the SGCP framework applied to the \code{cheng} dataset. This adjacency matrix serves as input to the \code{clustering} function, resulting in the clustering outcome stored in \code{resClus}. The clustering result, \code{resClus}, is subsequently utilized in the \code{geneOntology} function to derive \code{resInitialGO}, which captures the initial gene ontology (GO) enrichment results. The \code{resInitialGO} output is then processed through the \code{semiLabeling} function to produce \code{resSemiLabel}, indicating the semi-labeled genes based on their clustering characteristics. This semi-labeled information is further employed in the \code{semiSupervised} function, yielding \code{resSemiSupervised}, which includes the final supervised classification outcomes for the unremarkable genes. }

\format{An object of class \code{list} containing the semi-supervised information.}

\usage{data(resSemiSupervised)}

\details{\code{resSemiSupervised} is a list containin the following information.
\itemize{
    \item{\code{semiSupervised}: an object of \code{caret} for the training model.}
    \item{\code{prediction}: A vector of predicted labels for unremakable genes.}
    \item{\code{FinalLabeling}: a dataframe gene semil-label and final predicted labels.}}
}

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
