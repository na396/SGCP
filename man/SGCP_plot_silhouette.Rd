% plot silhouette
\name{SGCP_plot_silhouette}
\alias{SGCP_plot_silhouette}
\title{Cluster silhouette index chart in the SGCP Pipeline}
\description{
     Generates a chart displaying the cluster silhouette index in the SGCP pipeline.
}
\usage{
    SGCP_plot_silhouette(df, tit = "Gene Silhouette Index",
                        xname = "genes", yname = "silhouette index")
}
\arguments{
    \item{df}{The silhouette dataframe returned by the \code{clustering} function in the SGCP pipeline.}
    \item{tit}{Plot title (default: "Gene Silhouette Index")}
    \item{xname}{X-axis title (default: "genes")}
    \item{yname}{Y-axis title (default: "silhouette index")}
}

\details{
In order to plot silhouette index, \code{sil} argument in the \code{clustering} function must be set to \code{TRUE}.
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
}

\seealso{
    \code{\link{clustering}}
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}

\examples{

library(SGCP)
data(resClus)

## call the function
plt <- SGCP_plot_silhouette(df = resClus$silhouette)

print(plt)

}


