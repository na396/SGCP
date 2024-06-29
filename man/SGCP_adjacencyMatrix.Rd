% adjacencyMatrix description
\name{adjacencyMatrix}
\alias{adjacencyMatrix}
\title{Performs netwrok construction step in the SGCP pipeline}
\description{
    It creates the adjacency matrix of the gene co-expression network in the SGCP pipeline. Users can specify steps in the following order: calibration, norm, Gaussian kernel, and tom. If \code{calibration} is set to \code{TRUE}, SGCP performs calibration as the first step (refer to the manuscript for details). If \code{norm} is \code{TRUE}, each gene is normalized by its L2 norm. The Gaussian kernel metric is then calculated as a mandatory step to determine pairwise gene similarity values. If \code{tom} is \code{TRUE}, SGCP incorporates second-order node neighborhood information into the network. The pipeline concludes by returning a symmetric adjacency matrix \code{adja} of size m**n, where n is the number of genes. All values in the adjacency matrix range from 0 to 1, with 1 indicating maximum similarity. The diagonal elements of the matrix are set to zero.
}
\usage{
    adjacencyMatrix(expData, calibration = FALSE, norm = TRUE,
                    tom = TRUE, saveAdja = FALSE,
                    adjaNameFile = "adjacency.RData",
                    hm = "adjaHeatMap.png")
}
\arguments{
  \item{expData}{A dataframe or matrix containing the expression data, where rows correspond to genes and columns to samples.}
  \item{calibration}{Logical, default FALSE. If TRUE, performs calibration step.}
  \item{norm}{Logical, default TRUE. If TRUE, divides each gene (row) by its norm2.}
  \item{tom}{Logical, default TRUE. If TRUE, adds TOM to the network.}
  \item{saveAdja}{Logical, default FALSE. If TRUE, saves the adjacency matrix.}
  \item{adjaNameFile}{String indicating the name of the file for saving the adjacency matrix.}
  \item{hm}{String indicating the name of the file for saving the adjacency matrix heatmap.}
}
\value{
  \item{adja}{A symmetric matrix of dimension n * n representing the adjacency matrix, where n is the number of genes. Values range in (0, 1) with a zero diagonal.}
}
\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \href{https://www.frontiersin.org/articles/10.3389/fbinf.2021.704817/pdf}{calibration step information}
}

\examples{
## create an adjcency matrix
GeneExpression <- matrix(runif(1000, 0,1), nrow = 200, ncol = 5)
diag(GeneExpression) <- 0

## call the function
adja <- adjacencyMatrix(GeneExpression, hm= NULL)
head(adja)
}


\keyword{graphs}
