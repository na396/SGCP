% plot silhouette
\name{SGC_plot_silhouette}
\alias{SGC_plot_silhouette}
\title{Plots Gene Silhouette Index In SGC Pipeline}
\description{
    Plots the Silhouette index of genes in SGC pipeline.
}
\usage{
    SGC_plot_silhouette(df, tit = "Gene Silhouette Index",
                        xname = "genes", yname = "silhouette index")
}
\arguments{
    \item{df}{the silhouette dataframe returned by clustering function in SGC pipeline}
    \item{tit}{plot title, default "Gene Silhouette Index"}
    \item{xname}{x-axis title, default "genes"}
    \item{yname}{y-axis title, default "silhouette index"}
}

\details{
In order to plot silhouette index, `sil` parameter in clustering function must be set to TRUE.
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
