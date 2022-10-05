% semi-supervised
\name{semiSupervised}
\alias{semiSupervised}
\title{Performs Semi-supervised Step In SGC Pipeline}
\description{
    Performs semi-supervised classification.
}
\usage{
    semiSupervised(specExp, geneLab, model = "knn", kn = NULL, seed = NULL)
}
\arguments{
    \item{specExp}{matrix or dataframe with genes in rows and features in columns,
                   this is Y matrix from clustering function output.}
    \item{geneLab}{a dataframe returned by semiLabeling function,
                   contains the geneID and its corresponding label if is remarkable
                   otherwise NA.}
    \item{model}{either "knn" or "lr" for classification model,
                knn: k nearest neighbors,
                lr: logistic regression.}
    \item{kn}{an integer default NULL indicating the number of neighbors in knn,
              if kn is NULL, then
              kn = 20 : (20 + 2 * k) if 2 * k < 30 otherwise 20 : 30, where
              k is the numbber of remarkable cluster}
    \item{seed}{an integer for random state, default NULL.}
}

\details{
remarkable clusters are the clusters that have at least one remarkable gene.
}


\value{
    \item{semiSupervised}{an object of caret train class.}
    \item{prediction}{a vector of predicted labels for unremarkable genes.}
    \item{FinalLabeling}{a dataframe of geneID with its corresponding semilabel
        and final label.}

}

\references{
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGC: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{clustering}}
    \code{\link{semiLabeling}}
    \href{https://github.com/na396/SGC}{SGC Toturial}

}

\keyword{classification}
