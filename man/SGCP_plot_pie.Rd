%  plot geneOntology
\name{SGCP_plot_pie}
\alias{SGCP_plot_pie}
\title{Pie Chart of Gene Ontology Terms In SGCP Pipeline}
\description{
    Plots the test direction plus ontology of gene ontology terms in SGCP pipeline.
}
\usage{
    SGCP_plot_pie(df, tit = "GO Analysis",
                    xname = "module", yname = "count", posx = 1.9)

}
\arguments{
    \item{df}{the `GOresults` dataframe returned by geneOntology function in SGCP pipeline }
    \item{tit}{plot title, default "GO Analysis"}
    \item{xname}{x-axis title, default "module"}
    \item{yname}{y-axis title, default "count"}
    \item{posx}{a numeric for label position in pie chart, the higher the number the futhur the label will be from pie chart.}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{geneOntology}}
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}
\examples{

library(SGCP)
# load the output of geneOntology function
data(resInitialGO)

# call the function

plt <- SGCP_plot_pie(df = resInitialGO$GOresults)
print(plt)
}

