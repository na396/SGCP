%  plot jitter
\name{SGCP_plot_jitter}
\alias{SGCP_plot_jitter}
\title{P-value jitter chart for gene ontology enrichment in the SGCP pipeline.}
\description{
    Generates a jitter chart illustrating the cluster gene ontology enrichment in the SGCP pipeline.
}
\usage{
    SGCP_plot_jitter(df, tit = "p-values Distribution",
                    xname = "module", yname = "-log10 p-value", ps = 3)
}
\arguments{
    \item{df}{The \code{GOresults} dataframe returned by the \code{geneOntology} function in the SGCP pipeline.}
    \item{tit}{Plot title (default: "p-values Distribution")}
    \item{xname}{X-axis title (default: "module")}
    \item{yname}{Y-axis title (default: "-log10 p-value")}
    \item{ps}{Numeric value for point size (default: 3)}
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

plt <- SGCP_plot_jitter(df = resInitialGO$GOresults)
print(plt)
}
