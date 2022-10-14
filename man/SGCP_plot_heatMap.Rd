%  plot heatMap
\name{SGCP_plot_heatMap}
\alias{SGCP_plot_heatMap}
\title{Plots Adjacency Matrix HeatMap In SGCP Pipeline}
\description{
    Plots the HeatMap of Adjacency Matrix (Network) in SGCP pipeline.
}
\usage{
    SGCP_plot_heatMap(m, tit = "Adjacency Heatmap",
        xname = "genes", yname = "genes")
}
\arguments{
    \item{m}{an adjacency matrix returned by adjacencyMatrix function in SGCP pipeline, or must be a matrix, symmetric, with values in (0, 1) and zero diagonal}
    \item{tit}{plot title, default "Adjacency Heatmap"}
    \item{xname}{x-axis title, default "genes"}
    \item{yname}{y-axis title, default "genes"}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{adjacencyMatrix}}
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}
