% plot pca
\name{SGCP_plot_pca}
\alias{SGCP_plot_pca}
\title{PCA visualization in the SGCP Pipeline}
\description{
    Generates Principal Component Analysis (PCA) of the gene expression and transformed gene expression; comparision with and without labels.
}
\usage{
    SGCP_plot_pca(m, clusLabs, tit = "PCA plot", ps = .5)
}
\arguments{
    \item{m}{A numeric matrix of size n*m.}
    \item{clusLabs}{Either \code{NULL} or a vector of size n indicating cluster labels. There is a 1-to-1 correspondence between the rows in \code{m} and \code{clusLabs}.}
    \item{tit}{Plot title (default: "PCA plot")}
    \item{ps}{Point size (default: 0.5)}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}
\examples{

library(SGCP)
GeneExpression <- matrix(runif(200, 0,1), nrow = 40, ncol = 5)
diag(GeneExpression) <- 0

## call the function
plt <- SGCP_plot_pca(m = GeneExpression, clusLabs = NULL)

print(plt)

}
