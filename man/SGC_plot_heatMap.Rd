%  plot heatMap
\name{SGC_plot_heatMap}
\alias{SGC_plot_heatMap}
\title{Plots Adjacency Matrix HeatMap In SGC Pipeline}
\description{
    Plots the HeatMap of Adjacency Matrix (Network) in SGC pipeline.
}
\usage{
    SGC_plot_heatMap(m, tit = "Adjacency Heatmap",
        xname = "genes", yname = "genes")
}
\arguments{
    \item{m}{an adjacency matrix returned by adjacencyMatrix function in SGC pipeline, or must be a matrix, symmetric, with values in (0, 1) and zero diagonal}
    \item{tit}{plot title, default "Adjacency Heatmap"}
    \item{xname}{x-axis title, default "genes"}
    \item{yname}{y-axis title, default "genes"}
}


\value{
    returns the plot, an object of ggplot class.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{adjacencyMatrix}}
    \code{\link{SGC_ezPLOT}}
    \href{https://github.com/na396/SGC}{SGC Toturial}
}
\examples{

## lat mat be symmetric matrix of n*n with values in (0,1) and zero diagonal
## heatmap of mat
pl = SGC_plot_heatMap(m = mat)

}
