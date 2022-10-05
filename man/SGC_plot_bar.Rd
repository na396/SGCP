%  plot mean
\name{SGC_plot_bar}
\alias{SGC_plot_bar}
\title{Mean Over Gene Ontology Enrichment p-values In SGC Pipeline}
\description{
    Plots the mean over gene ontology enrchment p-values in SGC pipeline.
}
\usage{

    SGC_plot_bar(df, tit = "mean -log10 p-values",
                  xname = "module", yname = "-log10 p-value")
}
\arguments{
    \item{df}{the `GOresults` dataframe returned by geneOntology function in SGC pipeline }
    \item{tit}{plot title, default"mean -log10 p-values"}
    \item{xname}{x-axis title, default "module"}
    \item{yname}{y-axis title, default either "-log10 p-value"}
}


\value{
    returns the plot, an object of ggplot class.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{geneOntology}}
    \code{\link{SGC_ezPLOT}}
    \href{https://github.com/na396/SGC}{SGC Toturial}
}
