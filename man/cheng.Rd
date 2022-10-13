\name{cheng}
\docType{data}
\alias{cheng}
\title{Normalized gene expression Of ischemic cardiomyopathy (ICM) from a publication by Cheng et al. }
\description{
This data is a normalized gene expression of 1000 genes * 5 samples of ischemic cardiomyopathy (ICM) In Homo sapiens. 
ICM samples were obtained from the left ventricle of patients undergoing heart cardiac transplantation, whilst patients with no major cardiac history served as NF group.
gene expression were normalized using the DESeq method, based on median ratio of gene counts.
}
\format{An object of class \code{SummarizedExperiment}.}
\source{
    \url{https://www.sciencedirect.com/science/article/pii/S0010482520303061?via\%3Dihub}
}

\details{
assays contains the gene expression data and rowData field contains the corresponding gene Entrez IDs.
}
\usage{cheng}

\keyword{datasets}
