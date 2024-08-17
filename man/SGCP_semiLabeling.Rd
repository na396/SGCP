% semi-labeling
\name{semiLabeling}
\alias{semiLabeling}
\title{Performs gene semi-labeling step in the SGCP pipeline}
\description{
    Performs the Semi-labeling step in the SGCP pipeline to identify remarkable and unremarkable genes. This step involves collecting all gene ontology (GO) terms from all clusters and selecting terms in the top 0.1 percent. Genes associated with these terms are considered remarkable, while the remaining genes are categorized as unremarkable.
}
\usage{
    semiLabeling(geneID, df_GO, GOgenes, cutoff = NULL,
                percent = 0.10, stp = 0.01)
}
\arguments{
    \item{geneID}{A vector containing gene IDs, where n is the number of genes.}
    \item{df_GO}{The \code{GOresults} dataframe returned by the \code{geneOntology} function, containing information on GO terms in the clusters.}
    \item{GOgenes}{The \code{FinalGOTermGenes} list returned by the \code{geneOntology} function, listing genes associated with GO terms for each cluster.}
    \item{cutoff}{A numeric value in (0, 1) (default: NULL), serving as a baseline for GO term significance.}
    \item{percent}{A numeric value in (0, 1) (default: 0.1), indicating the percentile for selecting top GO terms.}
    \item{stp}{A numeric value in (0, 1) (default: 0.01), increment added to the \code{percent} parameter for stepwise selection of top GO terms.}
    }

\details{
Genes associated with GO terms more significant than the cutoff value are considered remarkable. If the cutoff value is \code{NULL}, SGCP determines the cutoff based on the significance level of the GO terms. Otherwise, SGCP selects the top percent (default: 0.1) of GO terms from all clusters combined, considering genes associated with these terms as remarkable. If all remarkable genes originate from a single cluster, SGCP incrementally increases the percent parameter by 0.01 to identify both remarkable and unremarkable genes. This process continues until remarkable genes originate from at least two clusters.
}

\value{
    \item{cutoff}{a numeric in (0,1) which indicates the selected cutoff.}
    \item{geneLabel}{a dataframe containing the information of geneID and
                its corresponding cluster label if is remarkable otherwise NA.}
}



\references{
    \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
}

\seealso{
    \code{\link{geneOntology}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}

\examples{

library(SGCP)
# load the output of clustering, gene ontology function

data(resClus)
data(resInitialGO)

# call the function

res <- semiLabeling(geneID = resClus$geneID, df_GO = resInitialGO$GOresults,
                GOgenes = resInitialGO$FinalGOTermGenes)
# cutoff value
res$cutoff

# gene semi-labeling information
head(res$geneLabel)

}
