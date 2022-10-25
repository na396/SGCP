\name{sgcp}
\docType{data}
\alias{sgcp}
\title{An example of output of \code{zSGCP} function in SGCP pipeline.}
\description{
This is an example of \code{ezSGPC} function output. This function is the automatic SGCP.}

\usage{data(sgcp)}
\format{An object of clas \code{list} containing the \code{ezSGCP} function information.}

\details{\code{sgcp} contains a list of clustering, initial.GO, semiLabeling, semiSupervised, final.GO fields, which contains the infromation of corresponding step.
\code{semilabel}:Boolean, indicates if semilabeling step is performed,
\code{clusterLabels}:a dataframe with geneID and its corresponding initial and final labels. In the clustering field,         following information is provided; \code{dropped.indices}:dropped gene indices, \code{geneID}:a vector of geneIDs,            \code{method}:indicates the selected method for number of cluster, \code{k}:selected number of clusters, \code{Y}:            transformed matrix with 2*k columns, \code{X}:eigenvalues correspond to 2*k columns in Y, \code{cluster}:object rof class     \code{kmeans}, \code{clusterLabels}:a vector containing the cluster label, for each gene,
there is a 1-to-1 correspondance between geneID and clusterLabes, \code{conductance}:a list containing mean and median,       andindividual cluster conductance index for clusters. In  each method, the clusterConductance field denote the cluster        label with its corresponding conductance index, \code{cvGOdf}:a dataframe used for gene ontology validation, for each        method, it shows the gene ontology enrichment on the cluster with smallest conductance index, \code{cv}:an string indicates the validation method for number of cluster, "cvGO" shows that gene ontology validation used, \code{clusterNumberPlot}:an     objet of class \code{ggplot2} for relativeGap"", "secondOrderGap", and "additiveGap", \code{silhouette}:a dataframe that     indicates the silhouette for genes, \code{original}:a list with matrix transformation and corresponding eigenvalues and n_egvec, where n_egvec top columns of tranformation is kept, initial.GO field, a list of \code{GOresults}:a dataframe containing the summary of the information of GOTerms, \code{FinalGOTermGenes}:a list containing the geneIDs of each GOTerms per cluster, semiLabeling field, a list of \code{cutoff}:a numeric in (0,1) which indicates the selected cutoff, \code{geneLabel}:a dataframe containing the information of geneID and its corresponding cluster label if is remarkable otherwise NA, semiSupervised field which is a list of \code{semiSupervised}:an object of classification result, \code{prediction}:a vector of predicted labels for unremarkable genes, \code{FinalLabeling}:a dataframe of geneID with its corresponding semilabel and final label final.GO field, a list of \code{GOresults}:a dataframe containing the summary of the information of GOTerms, \code{FinalGOTermGenes}:a list containing the geneIDs of each GOTerms per cluster.
}


\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
    \code{\link{ezSGCP}}
}
\examples{
library(SGCP)
data(sgcp)
summary(sgcp)

# clustering step
summary(sgcp$clustering)

# intial gene ontology step
summary(sgcp$initial.GO)

# semilabeling step
summary(sgcp$semiLabeling)

# semi-supervised step
summary(sgcp$semiSupervised)

# final gene ontology step
summary(sgcp$final.GO)

}
