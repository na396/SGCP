% adjacencyMatrix description
\name{adjacencyMatrix}
\alias{adjacencyMatrix}
\title{Performs Netwrok Construction step In SGC Pipeline}
\description{
    It creates the adjanecy matrix of gene co-expression network.
}
\usage{
    adjacencyMatrix(expData, calibration = FALSE, norm = TRUE,
                    tom = TRUE, saveAdja = FALSE,
                    adjaNameFile = "adjacency.RData",
                    hm = "adjaHeatMap.png")
}
\arguments{
    \item{expData}{a dataframe or matrix containing the expression data,
                    rows correspond to genes and columns to samples.}
    \item{calibration}{boolean, default FALSE,
                    if TRUE it performs calibration step.}
    \item{norm}{boolean, default TRUE,
                    if TRUE will divide each genes (rows) by its norm2.}
    \item{tom}{boolean, default TRUE,
                if TRUE it adds TOM to the network.}
    \item{saveAdja}{boolean, default FALSE,
                    if TRUE, the adjancency matrix will be saved.}
    \item{adjaNameFile}{string indicates the name of file for saving adjancency matrix.}
    \item{hm}{string indicates the name of file for saving adjancency matrix heat map.}
    }
\value{
    \item{adja}{a matrix of dimension n * n of the adjacency matrix
                where n is the number of the genes.
                This matrix is symmetric and entry values are in (0,1) with 0 diagonal}
}
\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \href{https://github.com/na396/SGC}{SGC Toturial}
}

\examples{
## create an adjcency matrix
GeneExpression = matrix(runif(0,1, 1000), nrow = 200, ncol = 5)
diag(GeneExpression) = 0

## call the function
adja = adjacencyMatrix(GeneExpression)
head(adja)
}


\keyword{graphs}
