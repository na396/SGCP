% semi-supervised
\name{semiSupervised}
\alias{semiSupervised}
\title{Performs the semi-supervised step in the SGCP pipeline}
\description{
    Performs the semi-supervised classification step in the SGCP pipeline. It utilizes the transformed matrix from the \code{clustering} function along with gene semi-labels from the \code{semiLabeling} function. The labeled (remarkable) genes serve as the training set to train either a "k-nearest neighbor" or "logistic regression" model. The trained model then predicts labels for unlabeled (unremarkable) genes, resulting in the final modules.
}
\usage{
    semiSupervised(specExp, geneLab, model = "knn", kn = NULL)
}
\arguments{
    \item{specExp}{Matrix or dataframe where genes are in rows and features are in columns, representing the Y matrix from \code{clustering} function output.}
    \item{geneLab}{A dataframe returned by the \code{semiLabeling} function, containing geneIDs and their corresponding labels (remarkable or NA).}
    \item{model}{Classification model type: "knn" for k-nearest neighbors or "lr" for logistic regression.}
    \item{kn}{An integer (default: NULL) indicating the number of neighbors in k-nearest neighbors (knn) model. If kn is \code{NULL}, the default value is determined by: kn = 20 if 2 * k < 30, otherwise kn = 20 : 30, where k is the number of remarkable clusters.}
    }

\details{
Remarkable clusters are defined as clusters that contain at least one remarkable gene.
}


\value{
    \\item{semiSupervised}{An object of the caret train class representing the semi-supervised classification model.}
    \item{prediction}{A vector containing predicted labels for unremarkable genes.}
    \item{FinalLabeling}{A dataframe containing geneIDs along with their corresponding semi-labels and final labels.}

}

\references{
    \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05848-w}{Aghaieabiane, N and Koutis, I (2024) SGCP: a spectral self-learning method for clustering genes in co-expression networks}
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
