\name{cheng}
\docType{data}
\alias{cheng}
\title{Normalized gene expression Of ischemic cardiomyopathy (ICM) from a publication by Cheng et al. }
\description{
This is a normalized gene expression data of 1000 genes * 5 samples. This data is a subset of a larger gene expression of ischemic cardiomyopathy (ICM) with 5000 genes and 57 samples. Gene expression data is normalized using the DESeq method, based on median ratio of gene counts.
}
\format{An object of class \code{SummarizedExperiment}.}
\source{
    \url{https://www.sciencedirect.com/science/article/pii/S0010482520303061?via\%3Dihub}
}

\details{
assays contains the gene expression data and rowData field contains the corresponding gene Entrez IDs. Sample names also are avialable in colData.
}
\usage{cheng}

\keyword{datasets}
