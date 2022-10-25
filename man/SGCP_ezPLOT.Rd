% SGCP_ezPLOT
\name{SGCP_ezPLOT}
\alias{SGCP_ezPLOT}
\title{Performs All SGCP Plots In One Step}
\description{
    On step plotting function for ezSGCP result. It takes the result from ezSGCP along with the expression data, and plot PCA of transformed and expression data, cluster conductance, gene silhouette index, method for number of clusters, distribution of gene ontology terms, density of gene ontology terms, cluster performance for both intial clusters and final modules.
}
\usage{
    SGCP_ezPLOT(sgcp, expreData, keep = FALSE,
            pdf.file = TRUE, pdfname = "ezSGCP.pdf",
            excel.file = TRUE, xlsxname = "ezSGCP.xlsx",
            w = 6, h = 6, sr = 2, sc = 2, ftype = "png",  uni = "in",
            expressionPCA = TRUE, pointSize1 = .5,
            exprePCATitle0 = "Expression Data PCA Without Labels",
            exprePCATitle1 = "Expression Data PCA With Initial Labels",
            exprePCATitle2 = "Expression Data PCA With Final Labels",
            transformedPCA = TRUE, pointSize2 = 0.5,
            transformedTitle0 = "Transformed Data PCA Without Labels",
            transformedTitle1 = "Transformed Data PCA Initial Labels",
            transformedTitle2 = "Transformed Data PCA Final Labels",
            conduct = TRUE,
            conductanceTitle = "Cluster Conductance Index",
            conductx = "clusterLabel", conducty = "conductance index",
            clus_num = TRUE,
            silhouette_index = FALSE,
            silTitle = "Gene Silhouette Index",
            silx = "genes", sily = "silhouette index",
            jitt1 = TRUE,
            jittTitle1 = "Initial GO p-values", jps1 = 3,
            jittx1 = "cluster", jitty1 = "-log10 p-value",
            jitt2 = TRUE,
            jittTitle2 = "Final GO p-values", jps2 = 3,
            jittx2 = "module", jitty2 = "-log10 p-value",
            density1 = TRUE,
            densTitle1 = "Initial GO p-values Density",
            densx1 = "cluster", densy1 = "-log10 p-value",
            density2 = TRUE,
            densTitle2 = "Final GO p-values Density",
            densx2 = "module", densy2 = "-log10 p-value",
            mean1 = TRUE,
            meanTitle1 = "Cluster Performance",
            meanx1 = "cluster", meany1 = "mean -log10 p-value",
            mean2 = TRUE,
            meanTitle2 = "Module Performance",
            meanx2 = "module", meany2 = "mean -log10 p-value",
            pie1 = TRUE, pieTitle1 = "Initial GO Analysis",
            piex1 = "cluster", piey1 = "count", posx1 = 1.8,
            pie2 = TRUE, pieTitle2 = "Final GO Analysis",
            piex2 = "module", piey2 = "count", posx2 = 1.8)
}
\arguments{

    \item{sgcp}{a returning result from ezSGCP function.}
    \item{expreData}{a matrix of initial gene expression dataset.}
    \item{keep}{Boolean, default FALSE. If TRUE plotting objects are kept.}
    \item{pdf.file}{Boolean, default TRUE, if TRUE it stores the plots in a pdf file.}
    \item{pdfname}{name of pdf file, default "ezSGCP.pdf".}
    \item{excel.file}{Boolean, default TRUE, if TRUE it stores the plots in a excel file.}
    \item{xlsxname}{name of pdf file, default "ezSGCP.xlsx".}
    \item{w}{width of plot images in excel, default 6.}
    \item{h}{height of plot images in excel, default 6.}
    \item{sr}{starting row in an excel sheet, default 2.}
    \item{sc}{starting column in an excel sheet, default 2.}
    \item{ftype}{plot image type, default "png".}
    \item{uni}{plot image units , default "in" for inch.}

    \item{expressionPCA}{Boolean, default TRUE, if TRUE PCA of gene expression
        data is plotted.}
    \item{pointSize1}{point size in for expression PCA, default 0.5.}

    \item{exprePCATitle0}{a string for expression PCA plot title without labels,
    default "Expression Data PCA Without Labels".}

    \item{exprePCATitle1}{a string for expression PCA plot title with initial cluster labels,
    default "Expression Data PCA With Initial Labels".}

    \item{exprePCATitle2}{a string for expression PCA plot title with final module labels,
    default "Expression Data PCA With Final Labels".}

    \item{transformedPCA}{Boolean, default TRUE, if TRUE PCA of transformed
        data is plotted.}
    \item{pointSize2}{point size in for transformed PCA, default 0.5.}

    \item{transformedTitle0}{a string for PCA plot title without labels on transformed data,
    default "Transformed Data PCA Without Labels".}

    \item{transformedTitle1}{a string for PCA plot title with initial cluster labels on transformed data,
    default "Transformed Data PCA Initial Labels".}

    \item{transformedTitle2}{a string for PCA plot title with final labels on transformed data,
    default "Transformed Data PCA Final Labels".}

    \item{conduct}{Boolean, default TRUE, if TRUE conductance indices for
    clusters are plotted.}
    \item{conductanceTitle}{a string for conductance indices plot title,
        default "Cluster Conductance Index".}
    \item{conductx}{a string for x-axis title in conductance indices plot,
        default "clusterLabel".}
    \item{conducty}{a string for y-axis title in conductance indices plot,
        default "conductance index".}

    \item{clus_num}{Boolean, default TRUE, if TRUE cluster numbers method are plotted.}


    \item{silhouette_index}{Boolean, default FALSE, if TRUE silhouette indices
    for genes are plotted.}
    \item{silTitle}{a string for silhouette indices plot title,
        default Gene Silhouette Index".}
    \item{silx}{a string for x-axis title in silhouette plot ,
        default "genes".}
    \item{sily}{a string for y-axis title in silhouette indices plot,
        default "silhouette index".}


    \item{jitt1}{Boolean, default TRUE, if TRUE jitter plot of p-values of GO terms
    in intial clusters are plotted.}
    \item{jps1}{point size in jiiter plot for initial clusters, default 3.}
    \item{jittTitle1}{a string for jiiter plot title for initial clusters,
        default "Initial GO p-values".}
    \item{jittx1}{a string for jiiter plot for initial clusters legend ,
        default "cluster".}
    \item{jitty1}{string for y-axis title in jiiter plot for initial clusters,
        default "-log10 p-value".}

    \item{jitt2}{Boolean, default TRUE, if TRUE jitter plot of p-values of GO terms
    in final modules are plotted.}
    \item{jps2}{point size in jiiter plot for final modules, default 3.}
    \item{jittTitle2}{a string for jiiter plot title for final modules,
        default "Final GO p-values".}
    \item{jittx2}{a string for jiiter plot for final modules legend ,
        default "module".}
    \item{jitty2}{string for y-axis title in jiiter plot for final modules,
        default "-log10 p-value".}


    \item{density1}{Boolean, default TRUE, if TRUE density plot of p-values of
    GO terms in initial clusters are plotted.}
    \item{densTitle1}{a string for density plot title for initial clusters,
        default "Initial GO p-values Density".}
    \item{densx1}{a string for density plot for initial clusters legend ,
        default "cluster".}
    \item{densy1}{string for y-axis title in desnity plot for initial clusters,
        default "-log10 p-value".}

    \item{density2}{Boolean, default TRUE, if TRUE density plot of p-values of GO terms
    in final modules are plotted.}
    \item{densTitle2}{a string for density plot title for final modules,
        default "Final GO p-values Density".}
    \item{densx2}{a string for density plot for final modules legend ,
        default "module".}
    \item{densy2}{string for y-axis title in desnity plot for final modules,
        default "-log10 p-value".}

    \item{mean1}{Boolean, default TRUE, if TRUE mean over p-values of GO terms
    in initial clusters are plotted.}
    \item{meanTitle1}{a string for mean plot title for initial clusters,
        default "Cluster Performance".}
    \item{meanx1}{a string for mean plot for initial clusters legend ,
        default "cluster".}
    \item{meany1}{string for y-axis title in mean plot for initial clusters,
        default "mean -log10 p-value".}

    \item{mean2}{Boolean, default TRUE, if TRUE mean over p-values of GO terms
    in final modules are plotted.}
    \item{meanTitle2}{a string for mean plot title for final modules,
        default "Module Performance".}
    \item{meanx2}{a string for mean plot for initial clusters legend ,
        default "module".}
    \item{meany2}{string for y-axis title in mean plot for final modules,
        default "mean -log10 p-value".}

    \item{pie1}{Boolean, default TRUE, if TRUE pie chart of direction and ontology
                of GO terms for initial clusters are plotted.}
    \item{pieTitle1}{a string for pie plot title for initial clusters,
        default "Initial GO Analysis".}
    \item{piex1}{a string for pie plot x-axis title for initial clusters,
        default "cluster".}
    \item{piey1}{string for y-axis title in pie plot for initial clusters,
        default "count".}
    \item{posx1}{a numeric, default 1.8,
                position of label of -log10 p-value of the most significant term.}

    \item{pie2}{Boolean, default TRUE, if TRUE pie chart of direction and ontology
                of GO terms for final modules are plotted.}
    \item{pieTitle2}{a string for pie plot title for initial clusters,
        default Final GO Analysis".}
    \item{piex2}{a string for pie plot x-axis title for final modules,
        default "module".}
    \item{piey2}{string for y-axis title in pie plot for final modules,
        default "count".}
    \item{posx2}{a numeric, default 1.8,
                position of label of -log10 p-value of the most significant term.}


}

\value{
    Returns the plotting object for each plot, if keep is TRUE.
    }

\references{
        \href{https://arxiv.org/abs/2209.10545}{Aghaieabiane, N and Koutis, I (2022) SGCP: A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks}
    }

\seealso{
    \href{https://github.com/na396/SGCP}{SGCP Toturial}
}

\examples{

library(SGCP)
library(SummarizedExperiment)

# load the result of ezSGCP function
data(sgcp)

# load the input expression dataset
data(cheng)
expData <- assay(cheng)

# to call the function uncomment the following
## plt <- SGCP_ezPLOT(sgcp = sgcp, expreData = cheng, keep = TRUE)

## print(plt)

}
