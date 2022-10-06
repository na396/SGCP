%  plot conductance
\name{SGC_plot_conductance}
\alias{SGC_plot_conductance}
\title{Plots Cluster Conductance Index In SGC Pipeline}
\description{
    Plots the cluster conductance index per cluster in SGC pipeline.
}
\usage{
    SGC_plot_conductance(conduct, tit = "Clustering Conductance Index",
                        xname = "cluster", yname = "conductance")
}
\arguments{
    \item{conduct}{conductance field returned by clustering function in SGC pipeline}
    \item{tit}{plot title, default "Clustering Conductance Index"}
    \item{xname}{x-axis title, default "cluster"}
    \item{yname}{y-axis title, default "conductance"}
}



\value{
    returns the plot, an object of ggplot class.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{clustering}}
    \code{\link{SGC_ezPLOT}}
    \href{https://github.com/na396/SGC}{SGC Toturial}
}
\examples{

## let conductance be the conductance field in output of clustering function
# to plot the conductance value
pl = SGC_plot_conductance(conduct = conductance)

}

