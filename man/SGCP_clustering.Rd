% networkClustering description
\name{clustering}
\alias{clustering}
\title{Perform network clustering step in the SGCP pipeline}
\description{
    It performs clustering on the adjacency network of gene co-expression network in SGCP pipeline. Initially, it transforms the n*n adjacency matrix into a new dimension Y. Subsequently, it determines the number of clusters k using three methods: "relativeGap", "secondOrderGap", and "additiveGap". For each method, k-means clustering is applied to Y with the determined k as input. Conductance indices are computed for the clusters within each method, and the cluster with the smallest conductance index is selected for further analysis. Following this, gene ontology enrichment analysis is performed on the selected clusters to finalize the optimal k. The pipeline concludes by returning the result of k-means clustering based on the selected method, along with the transformed matrix Y and additional information. This step produces the initial clusters.
}
\usage{
    clustering(adjaMat, geneID , annotation_db ,
                kopt = NULL, method = NULL,
                func.GO = sum, func.conduct = min,
                maxIter = 1e8, numStart = 1000, eff.egs = TRUE,
                saveOrig = TRUE, n_egvec = 200, sil = FALSE)
}
\arguments{
    \item{adjaMat}{A squared symmetric matrix of size n*n with values in (0, 1) and 0 diagonal. This is the output of the \code{adjacencyMatrix} function in SGCP.}
    \item{geneID}{A vector containing gene IDs of size n, where n is the number of genes.}
    \item{annotation_db}{A string indicating the genomic-wide annotation database.}
    \item{kopt}{An integer denoting the optimal number of clusters \( k \) chosen by the user (default: NULL).}
    \item{method}{Method for identifying the number of clusters \( k \) (default: NULL). Options include "relativeGap", "secondOrderGap", "additiveGap", or \code{NULL}.}
    \item{func.GO}{A function for gene ontology validation (default: sum).}
    \item{func.conduct}{A function for conductance validation (default: min).}
    \item{maxIter}{An integer specifying the maximum number of iterations for k-means clustering.}
    \item{numStart}{An integer indicating the number of starts for k-means clustering.}
    \item{eff.egs}{Boolean (default: TRUE). If TRUE, uses \code{eigs_sym} to calculate eigenvalues and eigenvectors, which is more efficient than R's default function.}
    \item{saveOrig}{Boolean (default: TRUE). If TRUE, keeps the transformation matrix.}
    \item{n_egvec}{An integer (default: 200) specifying the number of columns of the transformation matrix to retain. Should be less than 200.}
    \item{sil}{Boolean (default: FALSE). If TRUE, calculates silhouette index for each cluster.}


    }

\details{
If \code{kopt} is not null, SGCP will determine clusters based on the specified \code{kopt}. Otherwise, if \code{method} is not \code{NULL}, SGCP will select k using the specified method. If both \code{geneID} and \code{annotation_db} are \code{NULL}, SGCP will determine the optimal method and its corresponding number of clusters based on conductance validation. It selects a method where the conductance, evaluated by \code{func.conduct}, is minimized. Alternatively, SGCP defaults to gene ontology validation to find the optimal method and its corresponding clusters. It performs gene ontology enrichment on clusters, selecting the method where the cluster with the minimum conductance index yeilds the highest \code{func.GO} over log10 of the p-values.
}

\value{
    \itemize{
        \item{\code{dropped.indices}}{ A vector of dropped gene indices.}
        \item{\code{geneID}}{ A vector of gene IDs.}
        \item{\code{method}}{ Indicates the selected method for number of clusters.}
        \item{\code{k}}{ Selected number of clusters.}
        \item{\code{Y}}{ Transformed matrix with 2*k columns.}
        \item{\code{X}}{ Eigenvalues corresponding to 2*k columns in Y.}
        \item{\code{cluster}}{ An object of class \code{kmeans}.}
        \item{\code{clusterLabels}}{ A vector containing the cluster label per gene, with a 1-to-1 correspondence to geneID.}
        \item{\code{conductance}}{ A list containing mean and median conductance indices for clusters per method. The index in \code{clusterConductance} field denotes the cluster label and the value shows the conductance index.}
        \item{\code{cvGOdf}}{ A dataframe used for gene ontology validation. For each method, it returns the gene ontology enrichment result on the cluster with the minimum conductance index.}
        \item{\code{cv}}{ A string indicating the validation method for number of clusters}:
        \itemize{
            \item{"cvGO": Gene ontology validation used.}
            \item{"cvConductance": Conductance validation used.}
            \item{"userMethod": User-defined method.}
            \item{"userkopt": User-defined \code{kopt}.}
        }

        \item{\code{clusterNumberPlot}}{ An object of class \code{ggplot2} for relativeGap, secondOrderGap, and additiveGap.}
        \item{\code{silhouette}} {A dataframe indicating the silhouette index for genes.}
        \item{\code{original}}{ A list with matrix transformation, corresponding eigenvalues, and \code{n_egvec}, where the top \code{n_egvec} columns of the transformation are retained.}

    }
}

\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{adjacencyMatrix}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}

\examples{
## load cheng dataset
library(SGCP)
library(SummarizedExperiment)

data(cheng)
expData <- assay(cheng)
geneID <- rowData(cheng)
geneID <- geneID$ENTREZID


# to create the adjacency matrix un comment the following
## resAdja <- adjacencyMatrix(expData = expData, hm = NULL)
## resAdja[0:10, 0:5]

# to perform clustering
## library(org.Hs.eg.db)
annotation_db = "org.Hs.eg.db"
## resClus = clustering(adjaMat = resAdja, geneID = geneID,
##              annotation_db = annotation_db)


}

\keyword{clustering}
