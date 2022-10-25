% semi-supervised
\name{semiSupervised}
\alias{semiSupervised}
\title{Performs Semi-supervised step In SGCP Pipeline}
\description{
    Performs semi-supervised classification step in SGCP pipline. It takes the transformed matrix from clustering function along with gene semi-labels from semiLabeling function, and use the labeled (remarkable) genes as the training set to trian either "k nearest neighbor" or "logistic regression" model and make prediction for unlabeled (unremarkable genes). At this step, final modules are produced.
}
\usage{
    semiSupervised(specExp, geneLab, model = "knn", kn = NULL)
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
    \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
}

\seealso{
    \code{\link{clustering}}
    \code{\link{semiLabeling}}
    \href{https://github.com/na396/SGCP}{SGCP Toturial}

}
\examples{

library(SGCP)
# load the output of clustering, gene ontology function

data(resClus)
data(resSemiLabel)

# call the function

res <- semiSupervised(specExp = resClus$Y, geneLab = resSemiLabel$geneLabel)

# model summary
summary(res$semiSupervised)

# prediction label for unremarkable genes
head(res$prediction)

# semi and final gene labels
head(res$FinalLabeling)

}


\keyword{classification}
