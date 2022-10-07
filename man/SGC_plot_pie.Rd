%  plot geneOntology
\name{SGC_plot_pie}
\alias{SGC_plot_pie}
\title{Pie Chart of Gene Ontology Terms In SGC Pipeline}
\description{
    Plots the test direction plus ontology of gene ontology terms in SGC pipeline.
}
\usage{
    SGC_plot_pie(df, tit = "GO Analysis",
                    xname = "module", yname = "count", posx = 1.9)

}
\arguments{
    \item{df}{the `GOresults` dataframe returned by geneOntology function in SGC pipeline }
    \item{tit}{plot title, default "GO Analysis"}
    \item{xname}{x-axis title, default "module"}
    \item{yname}{y-axis title, default "count"}
    \item{posx}{a numeric for label position in pie chart, the higher the number the futhur the label will be from pie chart.}
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

