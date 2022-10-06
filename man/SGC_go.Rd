% gene ontology
\name{geneOntology}
\alias{geneOntology}
\title{Performs Gene Ontology Enrichment step In SGC Pipeline}
\description{
    It performs gene ontology enrichment for the each cluster using GOstat package.
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
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \href{https://github.com/na396/SGC}{SGC Toturial}
    \href{https://bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf}{GOstat Toturial}
}

\examples{

## let geneID is the vector of all gene ids, and geneClusLab is a vector of cluster label for each gene, and annotation_db be the genome wide annotation  database

# to call
resGO = geneOntology(geneUniv = geneID, clusLab = geneClusLab, annotation_db = annotation_db)

}

