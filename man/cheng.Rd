\name{cheng}
\docType{data}
\alias{cheng}
\title{Normalized gene expression data from Cheng et al.'s publication on ischemic cardiomyopathy (ICM).}
\description{
This dataset contains normalized gene expression data for 1500 genes across 5 samples. It is a subset of a larger dataset related to ischemic cardiomyopathy (ICM), which includes 5000 genes and 57 samples. The normalization was performed using the DESeq method, which utilizes the median ratio of gene counts to achieve normalization.
}
\format{An object of class \code{SummarizedExperiment}.}
\source{
    \url{https://www.sciencedirect.com/science/article/pii/S0010482520303061?via\%3Dihub}
}

\usage{data(cheng)}

\details{
assays contains the gene expression data and rowData field includes the corresponding gene Entrez IDs. Sample names also are avialable in colData.
}

\examples{
## load cheng dataset
library(SGCP)
library(SummarizedExperiment)

data(cheng)
expData <- assay(cheng)
geneID <- rowData(cheng)
geneID <- geneID$ENTREZID

}


\keyword{datasets}
