% networkClustering description
\name{clustering}
\alias{clustering}
\title{Perform Network Clustering step In SGCP Pipeline}
\description{
    It performs clustering on the adjacency network of gene co-expression network in SGCP pipeline. To this end, it firstly transforms the adjacency matrix of size n*n into new dimension Y of n*n using gaph Laplacian. It then, calculates number of cluster k based on three methods "relativeGap", "secondOrderGap", and "additiveGap". For each method, it performs kmeans on Y and the corresponding k as the input. Then it calculates conductance index for the clusters in the methods, and for each method, it picks the cluster that has smallest conductance index. Finally it performs gene ontology enrichment on those selected clusters to finalize k. At the end, it returns the result of the kmeans clustering on the selected method along with the transformed matrix Y, and some additional information. At this step, initial clusters are produced.
}
\usage{
    clustering(adjaMat, geneID , annotation_db ,
                kopt = NULL, method = NULL,
                func.GO = sum, func.conduct = min,
                maxIter = 1e8, numStart = 1000, eff.egs = TRUE,
                saveOrig = TRUE, n_egvec = 200, sil = FALSE)
}
\arguments{
    \item{adjaMat}{adjacency matrix of n*n where n is the number of the genes,
                    this matrix is squared symmetric with values in (0, 1)
                    and 0 diagonal.
                    It is the output of adjacencyMatrix function of SGCP.}
    \item{geneID}{a vector containing the genes IDs of size n
                where n is the number of genes.}
    \item{annotation_db}{a string indicating the genomic wide annotation database.}
    \item{kopt}{an integer denotes the optimal number of clusters k by the user,
                default is NULL.}
    \item{method}{method for identifying the number of clusters k, default NULL,
                either "relativeGap", "secondOrderGap", "additiveGap", or NULL.}
    \item{func.GO}{a function for gene ontology validation,
                    default is sum.}
    \item{func.conduct}{a function for conductance validation,
                        default is min.}
    \item{maxIter}{an integer,
                identifies the maximum number of iteration for kmeans.}
    \item{numStart}{an integer,
                    identifies the number of start for kmeans.}
    \item{eff.egs}{boolean, default TRUE,
                    if TRUE, used \code{eigs_sym} ti calculate the eigenvalues and eigenvectors, more efficient than R default                      function }
    \item{saveOrig}{boolean, default TRUE,
                    if TRUE, keeps the transformation matrix.}
    \item{n_egvec}{an integer less than 200, default = 200,
                indicates the number of columns of transformation matrix to be kept}
    \item{sil}{boolean, default FALSE,
            if TRUE, calculates silhouette index for each cluster.}

    }

\details{
If \code{kopt} is not null, SGCP will find clusters based on \code{kopt}. Otherwise,
if method is not null, SGCP will pick k based on the selected method. Otherwise,
if \code{geneID} or \code{annotation_db} is null, SGCP will determine the optimal method and its
corresponding number of cluster based on condutance validation. It picks a method
that \code{func.conduct} on its cluster is minimum. Otherwise,
SGCP will use gene ontology validation (by default) to find the optimal method and
its corresponding number of clusters. To this end, it will perform gene ontology
enrichment on the cluster with minimum conductanc index per method and pick the
one that has the maximum \code{func.GO} over -log10 of p-values.

}

\value{
a list containing some of the following depending on the initial call
    \item{dropped.indices}{a vector of dropped gene indices.}
    \item{geneID}{a vector of geneIDs.}
    \item{method}{indicates the selected method for number of cluster.}
    \item{k}{selected number of clusters.}
    \item{Y}{transformed matrix with 2*k columns.}
    \item{X}{eigenvalues correspond to 2*k columns in Y.}
    \item{cluster}{an object of class \code{kmeans}.}
    \item{clusterLabels}{a vector containing the cluster label per gene, there is a 1-to-1 correspondence                        between geneID and clusterLabes}
    \item{conductance}{a list containing mean and median, and individual cluster conductance index for                              clusters per method. Index in `clusterConductance` field denotes the cluster label                          and the value shows the conductance index.}
    \item{cvGOdf}{a dataframe used for gene ontology validation. For each method, it returns the gene                           ontology enrichment result on the cluster with minimum conductance index.}
    \item{cv}{an string indicates the validation method for number of cluster,
                "cvGO": if gene ontology validation used,
                "cvConductance": if conductance validation used,
                "userMethod": if user defined the method,
                "userkopt": if user defines the kopt.}
    \item{clusterNumberPlot}{an object of class \code{ggplot2} for
            relativeGap"", "secondOrderGap", and "additiveGap".}
    \item{silhouette}{a dataframe that indicates the silhouette for genes.}
    \item{original}{a list with matrix transformation and
                    corresponding eigenvalues and n_egvec,
                    where n_egvec top columns of transformation is kept.}

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
