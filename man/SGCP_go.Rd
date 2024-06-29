% gene ontology
\name{geneOntology}
\alias{geneOntology}
\title{Performs gene ontology enrichment step in the SGCP pipeline.}
\description{
    It performs gene ontology enrichment step GOstat package in SGCP pipeline. It takes the entire genes in the input with their labels, along with annotation_db to perform gene ontology enrichment for each set of genes that have similar label.
}
\usage{
    geneOntology(geneUniv, clusLab, annotation_db,
                direction = c("over", "under"),
                ontology = c("BP", "CC", "MF"), hgCutoff = NULL,
                cond = TRUE)
}
\arguments{
    \item{geneUniv}{a vector of all the geneIDs in the expression dataset.}
    \item{clusLab}{a vector of cluster label for each geneID.}
    \item{annotation_db}{a string indicating the genomic wide annotation database.}
    \item{direction}{test direction, default c("over", "under"),
                    for over-represented, or under-represented GO terms.}
    \item{ontology}{GO ontologies, default c("BP", "CC", "MF"),
                    BP: Biological Process,
                    CC: Cellular Component,
                    MF: Molecular Function.}
    \item{hgCutoff}{a numeric value in (0,1) as the p-value cutoff, default 0.05,
                    GO terms smaller than hgCutoff value are kept.}
    \item{cond}{Boolean, default TRUE,
                if TRUE conditional hypergeometric test is performed.}

    }


\value{
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
    \href{https://bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf}{GOstat Toturial}
}

\examples{

library(SGCP)
# load the output of clustering function
data(resClus)

# call the function
library(org.Hs.eg.db)

# to call the geneOntology uncomment the following
## res <- geneOntology(geneUniv = resClus$geneID, clusLab = resClus$clusterLabels,
##                     annotation_db = "org.Hs.eg.db")
## summary(res$GOresults)
## summary(res$FinalGOTermGenes)

}


