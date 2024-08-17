%  plot density
\name{SGCP_plot_density}
\alias{SGCP_plot_density}
\title{Visualization of gene ontology term p-value distribution in the SGCP pipeline}
\description{
    Generates a density chart displaying p-values of gene ontology terms in the SGCP pipeline.
}
\usage{
    SGCP_plot_density(df, tit = "p-values Density",
                    xname = "module", yname = "-log10 p-value")

}
\arguments{
    \item{df}{The \code{GOresults} dataframe returned by the \code{geneOntology} function in the SGCP pipeline.}
    \item{tit}{Plot title (default: "p-values Density")}
    \item{xname}{X-axis title (default: "module")}
    \item{yname}{Y-axis title (default: "-log10 p-value")}
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

plt <- SGCP_plot_density(df = resInitialGO$GOresults)
print(plt)
}
