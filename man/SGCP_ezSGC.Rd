% ezSGCP
\name{ezSGCP}
\alias{ezSGCP}
\title{Performs All SGCP pipeline In One Step}
\description{
    On step SGCP pipeline for gene co-expression network construction and analysis. It takes the gene expression and gene IDs, along with annotation_db and performs all steps of SGCP in a single function. It firstly perform network construction step and build the adjacency matrix. It then apply network clustering step using to find the intial clusters. Using gene ontology enrichent step, it finds and  divides the genes into set of remarkable and unremarkable and use this for the next step to semi-label the data and convert the problem into smei-supervised. It used the remarkable genes as the training set to train a supervised model and make prediction for unremarakable gene and produce the final modules. Finaly, it performs one more gene ontology step to see the module enrichment.
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
    \item{expData}{a dataframe or matrix containing the expression data,
        rows correspond to genes and columns to samples.}
    \item{geneID}{a vector containing the genes IDs of size n
        where n is the number of genes.}
    \item{annotation_db}{a string indicating the genomic wide annotation database.}
    \item{semilabel}{Boolean, default TRUE,
                    if TRUE, semilabeling step will be performed.}
    \item{calib}{boolean, default FALSE,
        if TRUE it performs calibration step.}
    \item{norm}{boolean, default TRUE,
        if TRUE will divide each genes (rows) by its norm2.}
    \item{tom}{boolean, default TRUE,
        if TRUE it adds TOM to the network.}
    \item{saveAdja}{boolean, default FALSE,
        if TRUE, the adjancency matrix will be saved.}
    \item{adjaNameFile}{string indicates the name of file for saving adjancency matrix.}
    \item{hm}{string indicates the name of file for saving adjancency matrix heat map.}
    \item{kopt}{an integer denotes the optimal number of clusters k by the user,
        default is NULL.}
    \item{method_k}{method for identifying the number of clusters k, default NULL,
        either "relativeGap", "secondOrderGap", "additiveGap", or NULL.}
    \item{f.GO}{a function for gene ontology validation,
                    default is sum.}
    \item{f.conduct}{a function for conductance validation,
                        default is min.}
    \item{maxIteration}{an integer,
                    identifies the maximum number of iteration for kmeans.}
    \item{numberStart}{an integer,
                    identifies the number of start for kmeans.}
    \item{eff.egs}{a boolean, default TRUE,
                    if TRUE it uses \code{eigs_sym} to calculate the eigenvalues and eigenvectors, more efficient than R default                      function}
    \item{saveOrig}{boolean, default TRUE,
                    if TRUE, keeps the transformation matrix.}
    \item{n_egvec}{either "all" or an integer, default = 100,
                indicates the number of columns of transformation matrix to be kept}
    \item{sil}{boolean, default FALSE,
            if TRUE, calculates silhouette index for each cluster.}

    \item{dir}{test direction, default c("over", "under"),
        for over-represented, or under-represented GO terms.}
    \item{onto}{GO ontologies, default c("BP", "CC", "MF"),
        BP: Biological Process,
        CC: Cellular Component,
        MF: Molecular Function.}
    \item{hgCut}{a numeric value in (0,1) as the p-value cutoff, default 0.05,
        GO terms smaller than hgCutoff value are kept.}
    \item{condTest}{Boolean, default TRUE,
        if TRUE conditional hypergeometric test is performed.}
    \item{cutoff}{a numeric in (0, 1) default NULL,
        is a base line for GO term significancy.}
    \item{percent}{a number in (0,1) default 0.1,
        indicate the percentile for finding top GO terms.}
    \item{stp}{a number in (0,1) default 0.01,
        indicates increasing value to be added to percent parameter.}
    \item{model}{either "knn" or "lr" for classification model,
        knn: k nearest neighbors,
        lr: logistic regression.}
    \item{kn}{an integer default NULL indicating the number of neighbors in knn,
        if kn is NULL, then
        kn = 20 : (20 + 2 * k) if 2 * k < 30 otherwise 20 : 30, where
        k is the numbber of remarkable cluster}
}

\details{
For clustering step, if kopt is not null, SGCP will find clusters based on kopt. Otherwise,
if method is not null, SGCP will pick k based on the selected method. Otherwise,
if geneID or annotation_db is null, SGCP will determine the optimal method and its
corresponding number of cluster based on condutance validation. It picks a method
that func.conduct on its cluster is minimum. Otherwise,
SGCP will use gene ontology validation (by default) to find the optimal method and
its corresponding number of clusters. To this end, it will perform gene ontology
enrichment on the cluster with minimum conductanc index per method and pick the
one that has the maximum func.GO over -log10 of p-values. In semilabeling step,
gene associated to the GO terms more significant than cutoff value are remarkable.
If cutoff value is NULL, SGCP will find the cutoff depend on the GO terms significant
level. Otherwise, SGCP picks the top percent (by default 0.1) GO terms from all
clusters collectively, and consider the genes associated to those as remakable.
If all remarkable genes come from a single cluster,
then SGCP will increase the precent by 0.01 to find the remarkable and
unremarkable genes. It repeats this process until all remarkable genes come from
at least two clusters. In semi-supervise step, remarkable clusters are the
clusters that have at least one remarkable gene.
}

\value{
It returns a list of clustering, initial.GO, semiLabeling, semiSupervised,
final.GO fields, which contains the infromation of corresponding step.
    \item{semilabel}{Boolean, indicates if semilabeling step is performed.}
    \item{clusterLabels}{a dataframe with geneID and its corresponding initial
    and final labels.}
    clustering field, a list of
    \item{dropped.indices}{dropped gene indices.}
    \item{geneID}{a vector of geneIDs.}
    \item{method}{indicates the selected method for number of cluster.}
    \item{k}{selected number of clusters.}
    \item{Y}{transformed matrix with 2*k columns.}
    \item{X}{eigenvalues correspond to 2*k columns in Y.}
    \item{cluster}{object rof class \code{kmeans}.}
    \item{clusterLabels}{a vector containing the cluster label, for each gene,
                        there is a 1-to-1 correspondance between geneID and clusterLabes.}
    \item{conductance}{a list containing mean and median, and
                    individual cluster conductance index for clusters. In
                    each method, the clusterConductance field denote the cluster label
                    with its corresponding conductance index.}
    \item{cvGOdf}{a dataframe used for gene ontology validation,
                for each method, it shows the gene ontology enrichment
                on the cluster with smallest conductance index.}
    \item{cv}{an string indicates the validation method for number of cluster,
                "cvGO": if gene ontology validation used,
                "cvConductance": if conductance validation used,
                "userMethod": if user defined the method,
                "userkopt": if user defines the kopt}
    \item{clusterNumberPlot}{an objet of class \code{ggplot2} for
            relativeGap"", "secondOrderGap", and "additiveGap".}
    \item{silhouette}{a dataframe that indicates the silhouette for genes.}
    \item{original}{a list with matrix transformation and
                    corresponding eigenvalues and n_egvec,
                    where n_egvec top columns of tranformation is kept.}
    initial.GO field, a list of
    \item{GOresults}{a dataframe containing the summary of the information of GOTerms,
                clusterNum: indicates the cluster label,
                GOtype: indicates the test directions plut ontology,
                GOID: unique GO term id,
                Pvalue: the p-value of hypergeometric test for the GO term,
                OddsRatio: the odds ratio of the GO term,
                ExpCount: expected count value for genes associated the GO term,
                Count: actual count of the genes associated to the GO term in the cluster,
                Size: actual size of the genes associated to the GO term in the entire geneIDs,
                Term: description of the GO term.}

    \item{FinalGOTermGenes}{a list containing the geneIDs of each GOTerms per cluster.}
    semiLabeling field, a list of
    \item{cutoff}{a numeric in (0,1) which indicates the selected cutoff.}
    \item{geneLabel}{a dataframe containing the information of geneID and
                its corresponding cluster label if is remarkable otherwise NA.}
    semiSupervised field which is a list of
    \item{semiSupervised}{an object of classification result.}
    \item{prediction}{a vector of predicted labels for unremarkable genes.}
    \item{FinalLabeling}{a dataframe of geneID with its corresponding semilabel
        and final label.}
    final.GO field, a list of
    \item{GOresults}{a dataframe containing the summary of the information of GOTerms,
                clusterNum: indicates the cluster label,
                GOtype: indicates the test directions plut ontology,
                GOID: unique GO term id,
                Pvalue: the p-value of hypergeometric test for the GO term,
                OddsRatio: the odds ratio of the GO term,
                ExpCount: expected count value for genes associated the GO term,
                Count: actual count of the genes associated to the GO term in the cluster,
                Size: actual size of the genes associated to the GO term in the entire geneIDs,
                Term: description of the GO term.}

    \item{FinalGOTermGenes}{a list containing the geneIDs of each GOTerms per cluster.}


}

\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
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
