% adjacencyMatrix description
\name{adjacencyMatrix}
\alias{adjacencyMatrix}
\title{Performs Netwrok Construction step In SGCP Pipeline}
\description{
    It creates the adjanecy matrix of gene co-expression network in SGCP pipeline. Here, users can select steps. The order of steps are calibration, norm, Gaussian kernel, and tom. If calibration is TRUE, then SGCP will perform calibration at first (please see the manuscript for more information). Next, if norm is TRUE, SGCP will divide each gene by its norm2. Then, SGCP will calculate Gaussian kernel metric as the similarity function to calculate the pairwise gene similarity values, this step is mandatory and user cannot change it. If tom is TRUE, SGCP will add the information of the second order of the node neighborhood to the network. At the end, SGCP returns a symmetric adjacency matrix adja of size n*n where n is the number of genes. All values in the adjancency matrix range from 0 to 1, where 1 is the most similar. The diagonal of the returning matrix is zero.
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
adja <- adjacencyMatrix(GeneExpression)
head(adja)
}


\keyword{graphs}
