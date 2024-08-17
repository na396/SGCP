% ezSGCP
\name{ezSGCP}
\alias{ezSGCP}
\title{Integrated execution of the SGCP pipeline}
\description{
The SGCP pipeline for gene co-expression network construction and analysis integrates multiple steps into a single function. It begins with network construction, where gene expression data and gene IDs are utilized alongside an annotation database to build an adjacency matrix. Next, network clustering identifies initial clusters. Gene ontology enrichment distinguishes genes into remarkable and unremarkable sets, enabling semi-labeling to convert the problem into semi-supervised learning. Remarkable genes serve as the training set for a supervised model, predicting labels for unremarkable genes and producing final modules. Finally, another gene ontology step evaluates module enrichment.
}
\usage{
    ezSGCP(expData, geneID, annotation_db, semilabel = TRUE,
        calib = FALSE, norm = TRUE, tom = TRUE,
        saveAdja = FALSE, adjaNameFile = "adjacency.Rdata",
        hm = "adjaHeatMap.png",
        kopt = NULL, method_k = NULL, f.GO = sum, f.conduct = min,
        maxIteration = 1e8, numberStart = 1000, eff.egs = TRUE,
        saveOrig = TRUE, n_egvec = 100, sil = FALSE,
        dir = c("over", "under"), onto = c("BP", "CC", "MF"),
        hgCut = NULL, condTest = TRUE,
        cutoff = NULL, percent = 0.10, stp = 0.01,
        model = "knn", kn = NULL)
}
\arguments{
  \item{expData}{A dataframe or matrix containing the expression data, where rows correspond to genes and columns to samples.}
  \item{geneID}{A vector containing the gene IDs of size n, where n is the number of genes.}
  \item{annotation_db}{A string indicating the genomic-wide annotation database.}
  \item{semilabel}{Logical, default \code{TRUE}. If \code{TRUE}, performs semilabeling step.}
  \item{calib}{Logical, default \code{FALSE}. If \code{TRUE}, performs calibration step.}
  \item{norm}{Logical, default \code{TRUE}. If \code{TRUE}, divides each gene (row) by its norm2.}
  \item{tom}{Logical, default \code{TRUE}. If \code{TRUE}, adds TOM to the network.}
  \item{saveAdja}{Logical, default \code{FALSE}. If \code{TRUE}, saves the adjacency matrix.}
  \item{adjaNameFile}{String indicating the name of the file for saving the adjacency matrix.}
  \item{hm}{String indicating the name of the file for saving the adjacency matrix heatmap.}
  \item{kopt}{An integer indicating the optimal number of clusters k chosen by the user, default is \code{NULL}.}
  \item{method_k}{Method for identifying the number of clusters k, default \code{NULL}. Options are "relativeGap", "secondOrderGap", "additiveGap", or \code{NULL}.}
  \item{f.GO}{A function for gene ontology validation, default is sum.}
  \item{f.conduct}{A function for conductance validation, default is min.}
  \item{maxIteration}{An integer indicating the maximum number of iterations for kmeans.}
  \item{numberStart}{An integer indicating the number of starts for kmeans.}
  \item{eff.egs}{Boolean, default \code{TRUE}. If \code{TRUE}, uses \code{eigs_sym} to calculate eigenvalues and eigenvectors, which is more efficient than R's default function.}
  \item{saveOrig}{Boolean, default \code{TRUE}. If \code{TRUE}, keeps the transformation matrix.}
  \item{n_egvec}{Either "all" or an integer indicating the number of columns of the transformation matrix to keep, default is 100.}
  \item{sil}{Logical, default \code{FALSE}. If \code{TRUE}, calculates silhouette index for each cluster.}
  \item{dir}{Test direction for GO terms, default c("over", "under").}
  \item{onto}{GO ontologies to consider, default c("BP", "CC", "MF").}
  \item{hgCut}{Numeric value in (0,1) as the p-value cutoff for GO terms, default 0.05.}
  \item{condTest}{Logical, default \code{TRUE}. If \code{TRUE}, performs conditional hypergeometric test.}
  \item{cutoff}{Numeric in (0, 1) default \code{NULL}, baseline for GO term significance.}
  \item{percent}{Numeric in (0,1) default 0.1, percentile for finding top GO terms.}
  \item{stp}{Numeric in (0,1) default 0.01, increasing value for percent parameter.}
  \item{model}{Type of classification model, either "knn" (k nearest neighbors) or "lr" (logistic regression).}
  \item{kn}{Integer indicating the number of neighbors in knn, default \code{NULL}.}
}
\details{
For clustering step; If \code{kopt} is not \code{NULL}, SGCP finds clusters based on \code{kopt}. If \code{method_k} is not \code{NULL}, SGCP picks \code{k} based on the selected method (\code{"relativeGap"}, \code{"secondOrderGap"}, \code{"additiveGap"}). If \code{geneID} or \code{annotation_db} is \code{NULL}, SGCP determines the optimal method and corresponding number of clusters based on conductance validation. It selects the method where \code{func.conduct} on its clusters is minimized.
Otherwise, SGCP uses gene ontology validation (by default) to find the optimal method and its corresponding number of clusters. It performs gene ontology enrichment on the cluster with the minimum conductance index per method and selects the method that maximizes \code{func.GO} over -log10 of p-values.


For semilabeling step; Genes associated with GO terms more significant than \code{cutoff} value are considered remarkable. If \code{cutoff} value is \code{NULL}, SGCP determines the cutoff based on the significance level of GO terms.
SGCP selects the top \code{percent} (default 0.1) GO terms from all clusters collectively and considers genes associated with those as remarkable. If all remarkable genes come from a single cluster, SGCP increases the percent by 0.01 until remarkable genes come from at least two clusters.


For semi-supervise step; Remarkable clusters are those that have at least one remarkable gene. SGCP performs semi-supervised classification using the transformed matrix from clustering and gene semilabels from semilabeling function. It uses remarkable genes as the training set to train either a "k nearest neighbor" (knn) or "logistic regression" (lr) model and makes predictions for unremarkable genes to produce the final modules.}



\value{
    Returns a list with the following fields, depending on the initial call:
    \itemize{
      \item{\code{semilabel}}{ Boolean indicating if semilabeling step was performed.}
      \item{\code{clusterLabels}}{ DataFrame with geneID and its corresponding initial and final labels.}
      \item{\code{clustering}}{ List containing clustering information:}
      \itemize{
        \item{\code{dropped.indices}}{ Vector of dropped gene indices.}
        \item{\code{geneID}}{ Vector of geneIDs.}
        \item{\code{method}}{ Method selected for number of clusters.}
        \item{\code{k}}{ Selected number of clusters.}
        \item{\code{Y}}{ Transformed matrix with 2*k columns.}
        \item{\code{X}}{ Eigenvalues corresponding to 2*k columns in Y.}
        \item{\code{cluster}}{ Object of class \code{kmeans}.}
        \item{\code{clusterLabels}} {Vector containing cluster labels for each gene.}
        \item{\code{conductance}}{ List with mean, median, and individual cluster conductance indices. \code{clusterConductance} field denotes the cluster label and its conductance index.}
        \item{\code{cvGOdf}}{ DataFrame used for gene ontology validation. For each method, shows GO enrichment on the cluster with smallest conductance index.}
        \item{\code{cv}}{ String indicating validation method for number of clusters: "cvGO", "cvConductance", "userMethod", or "userkopt".}
        \item{\code{clusterNumberPlot}}{ Object of class \code{ggplot2} for visualizing cluster number selection methods.}
        \item{\code{silhouette}}{ DataFrame indicating silhouette indices for genes.}
        \item{\code{original}}{ List with matrix transformation, corresponding eigenvalues, and \code{n_egvec} top columns of transformation matrix kept.} }

       \item{\code{initial.GO}}{ List containing initial gene ontology (GO) information:}
    \itemize{
        \item{\code{GOresults}}{ DataFrame summarizing GO term information. Includes clusterNum, GOtype, GOID, Pvalue, OddsRatio, ExpCount, Count, Size, and Term.}
        \item{\code{FinalGOTermGenes}}{ List of geneIDs associated with each GO term per cluster.}}

    \item{\code{semiLabeling}}{ List containing semilabeling information:}
    \itemize{
        \item{\code{cutoff}}{ Numeric (0,1) indicating selected cutoff for significant GO terms.}
        \item{\code{geneLabel}}{ DataFrame with geneID and its corresponding cluster label if remarkable, otherwise NA.}}

    \item{\code{semiSupervised}}{ List containing semi-supervised classification information:}
    \itemize{
        \item{\code{semiSupervised}}{ Object of classification result.}
        \item{\code{prediction}}{ Vector of predicted labels for unremarkable genes.}
        \item{\code{FinalLabeling}}{ DataFrame of geneID with its corresponding semilabel and final label.}}

    \item{\code{final.GO}}{ List containing final gene ontology (GO) information:}
    \itemize{
        \item{\code{GOresults}}{ DataFrame summarizing GO term information. Includes \code{clusterNum}, \code{GOtype}, \code{GOID}, \code{Pvalue}, \code{OddsRatio}, \code{ExpCount}, \code{Count}, \code{Size}, and \code{Term}.}
        \item{\code{ FinalGOTermGenes}}{List of geneIDs associated with each GO term per cluster.}}

  }
}



\references{
   \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
}



\seealso{
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


library(org.Hs.eg.db)

# to call the function uncomment the following
## res <- ezSGCP(expData = expData, geneID = geneID, annotation_db = "org.Hs.eg.db")
## summary(res)
## summary(res$clustering)
## summary(res$initial.GO)
## summary(res$semiLabeling)
## summary(res$semiSupervised)
## summary(res$final.GO)





}
