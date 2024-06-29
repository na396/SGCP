%  plot conductance
\name{SGCP_plot_conductance}
\alias{SGCP_plot_conductance}
\title{Cluster conductance index bar chart in the SGCP Pipeline}
\description{
    Generates a bar chart displaying the cluster conductance index in the SGCP pipeline.
}
\usage{
    SGCP_plot_conductance(conduct, tit = "Clustering Conductance Index",
                        xname = "cluster", yname = "conductance")
}
\arguments{
    \item{conduct}{The conductance field returned by the \code{clustering} function in the SGCP pipeline.}
    \item{tit}{Plot title (default: "Clustering Conductance Index")}
    \item{xname}{X-axis title (default: "cluster")}
    \item{yname}{Y-axis title (default: "conductance")}
}



\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{clustering}}
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}

\examples{

library(SGCP)
# load the output of geneOntology function
data(resClus)

# call the function

plt <- SGCP_plot_conductance(conduct = resClus$conductance)
print(plt)
}
