% semi-labeling
\name{semiLabeling}
\alias{semiLabeling}
\title{Performs Gene Semi-labeling step In SGC Pipeline}
\description{
    Indentifies remarkable and unremarkable genes
}
\usage{
    semiLabeling(geneID, df_GO, GOgenes, cutoff = NULL,
    percent = 0.10, stp = 0.01)
}
\arguments{
    \item{geneID}{a vector containing the genes IDs of size n,
                  where n is the number of genes.}
    \item{df_GO}{GOresults dataframe returned by geneOntology function,
                 consists the information of GO terms the clusters.}
    \item{GOgenes}{FinalGOTermGenes list returedn by geneOntology function,
                   is a list of genes associated to the GO Terms per each cluster.}
    \item{cutoff}{a numeric in (0, 1) default NULL,
                  is a base line for GO term significancy.}
    \item{percent}{a number in (0,1) default 0.1,
                   indicate the percentile for finding top GO terms.}
    \item{stp}{a number in (0,1) default 0.01,
                indicates increasing value to be added to percent parameter.}
    }

\details{
Gene associated to the GO terms more significant than cutoff value are remarkable.
If cutoff value is NULL, SGC will find the cutoff depend on the GO terms significant
level. Otherwise, SGC picks the top percent (by default 0.1) GO terms from all
clusters collectively, and consider the genes associated to those as remakable.
If all remarkable genes come from a single cluster,
then SGC will increase the precent by 0.01 to find the remarkable and
unremarkable genes. It repeats this process until all remarkable genes come from
at least two clusters.
}

\value{
    \item{cutoff}{a numeric in (0,1) which indicates the selected cutoff.}
    \item{geneLabel}{a dataframe containing the information of geneID and
                its corresponding cluster label if is remarkable otherwise NA.}
}



\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{geneOntology}}
    \href{https://github.com/na396/SGC}{SGC Toturial}
}

\examples{

## let GODF is the dataframe of GOresults from geneOntology function
## let geneID is a vector of all gene ids
## let GOgenesDF is the dataframe of FinalGOTermGenes from geneOntology function

# to call semiLabeling
resSemiLabel = semiLabeling(geneID = geneID, df_GO = GODF, GOgenes = GOgenesDF)

}
