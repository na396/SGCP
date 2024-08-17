%  plot geneOntology
\name{SGCP_plot_pie}
\alias{SGCP_plot_pie}
\title{Gene ontology analysis pie chart in the SGCP pipeline}
\description{
    Generate a pie chart illustrating the ontology and test direction of gene ontology terms across the SGCP pipeline
}
\usage{
    SGCP_plot_pie(df, tit = "GO Analysis",
                    xname = "module", yname = "count", posx = 1.9)

}
\arguments{
    \item{df}{The \code{GOresults} dataframe returned by the \code{geneOntology} function in the SGCP pipeline.}
    \item{tit}{Plot title (default: "GO Analysis")}
    \item{xname}{X-axis title (default: "module")}
    \item{yname}{Y-axis title (default: "count")}
    \item{posx}{Numeric value for label position in the pie chart. A higher number places labels further from the pie chart.}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
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

