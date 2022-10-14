% plot pca
\name{SGCP_plot_pca}
\alias{SGCP_plot_pca}
\title{Plots PCA Of The Data In SGCP Pipeline}
\description{
    Plots PCA with and without labels in SGCP pipeline.
}
\usage{
    SGCP_plot_pca(m, clusLabs, tit = "PCA plot", ps = .5)
}
\arguments{
    \item{m}{a numeric matrix of n*m}
    \item{clusLabs}{NULL or a vector of size n showing cluster labels,
            there 1-to-1 correspondance between the rows in m and clusLabs}
    \item{tit}{plot title, default "PCA plot"}
    \item{ps}{point size, default .5}
}


\value{
    returns the plot, an object of class \code{ggplot2}.
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{SGCP_ezPLOT}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}
