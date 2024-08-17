%  plot heatMap
\name{SGCP_plot_heatMap}
\alias{SGCP_plot_heatMap}
\title{Adjacency matrix heatmap in the SGCP pipeline}
\description{
    Generates a Heatmap of the the adjacency matrix (network) in the SGCP pipeline.
}
\usage{
    SGCP_plot_heatMap(m, tit = "Adjacency Heatmap",
        xname = "genes", yname = "genes")
}
\arguments{
    \item{m}{An adjacency matrix returned by the \code{adjacencyMatrix} function in the SGCP pipeline, or any symmetric matrix with values in (0, 1) excluding the diagonal.}
    \item{tit}{Plot title (default: "Adjacency Heatmap")}
    \item{xname}{X-axis title (default: "genes")}
    \item{yname}{Y-axis title (default: "genes")}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
}

\seealso{
    \code{\link{adjacencyMatrix}}
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}
\examples{

library(SGCP)
GeneExpression <- matrix(runif(200, 0,1), nrow = 40, ncol = 5)
diag(GeneExpression) <- 0

## call the function
adja <- adjacencyMatrix(GeneExpression)

plt <- SGCP_plot_heatMap(m =  adja)
print(plt)

}
