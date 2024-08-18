# SGCP: a spectral self-learning method for clustering genes in co-expression networks, [link](https://link.springer.com/article/10.1186/s12859-024-05848-w)


## SGCP Introduction
The Self-training Gene Clustering Pipeline (`SGCP`) is an innovative framework for constructing and analyzing gene co-expression networks. Its primary objective is to group genes with similar expression patterns into cohesive clusters, often referred to as modules. SGCP introduces several novel steps that enable the computation of highly enriched gene modules in an unsupervised manner. What sets SGCP apart from existing frameworks is its integration of a semi-supervised clustering approach, which leverages Gene Ontology (GO) information. This unique step significantly enhances the quality of the resulting modules, producing highly enriched and biologically relevant clusters.

## SGCP Publication
`SGCP` is available at [BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-024-05848-w). 

## SGCP Installation
For detailed instructions and steps, please refer to the `SGCP` manual on
[Bioconductor page](https://bioconductor.org/packages/release/bioc/html/SGCP.html). To install the latest version of `SGCP`, you can access the GitHub repository using the following command:
```{r}
#install.packages("devtools")
#devtools::install_github("na396/SGCP")
```
## SGCP license
GPL-3

## SGCP encoding
UTF-8


## SGCP Input

`SGCP` requires three main inputs; __expData__ , __geneID__, and __annotation_db__.
*    __expData__: This is a matrix or dataframe of size `m*n` where `m` represents the number of genes and `n` represents the number of samples. It can contain data from either DNA-microarray or RNA-seq experiments . Note that `SGCP` assumes that pre-processing steps, such as normalization and batch effect corection, have already been performed, as these are not handled by the pipeline.
*    __geneID__: A vector of gene identifier corresponding to the rows in __expData__.
*    __anotation_db__: The name of a genome-wide annotation package for the organism of interest, used in the gene ontology (GO) enrichment step. The `annotation_db` package must be installed by user prior to using `SGCP`.

Below are some commonly used `annotation_db` packages along with their corresponding gene identifiers for different organisms.

|organism                     | annotation_db  | gene identifier         |
|:----------------------------|:--------------:|:---------------------   | 
|Homo sapiens (Hs)            | org.Hs.eg.db   | Entrez Gene identifiers |
|Drosophila melanogaster (Dm) | org.Dm.eg.db   | Entrez Gene identifiers |
|Rattus norvegicus (Rn)       | org.Rn.eg.db   | Entrez Gene identifiers |
|Mus musculus (Mm)            | org.Mm.eg.db   | Entrez Gene identifiers |
|Arabidopsis thaliana (At)    | org.At.tair.db | TAIR identifiers        |

Gene expression datasets for your analysis can be obtained from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/), a public repository of high-throughput gene expression data.


### SGCP Input Cleaning
In `SGCP`, the following assumptions are made about the input genes:

* Genes must have expression values available across all samples, with no missing values.
* Genes must exhibit non-zero variance in expression across all samples.
* ach gene must have exactly one unique identifier, specified by __geneID__.
* Genes must be annotated with Gene Ontology (GO) terms.


## SGCP Input Example
Here, we give a brief example of the `SGCP` input. For this documentation, we use the gene expression `GSE181225`. For more information visit its [Gene Expression Omnibus page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181225)). 

Throughout this section, several Bioconductor packages will be required. Make sure to install and load them as needed to follow the example.

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "GEOquery", "AnnotationDbi"))
```

First, set the directory
```{r}
# Display the current working directory
print(getwd())

# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "."
setwd(workingDir)
```

First, we need to download the gene expression file. The R package `GEOquery` is used to obtain gene expression data from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). For detailed information on how to use `GEOquer`y, refer to the [GEOquery guide](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html). 


To download the expression file for `GSE181225`, visit its [Gene Expression Omnibus page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181225). On the page, locate the file` GSE181225_LNCaP_p57_VO_and_p57_PIM1_RNA_Seq_normalizedCounts.txt.gz` in the `Supplementary files` section, which contains the normalized gene expression data. Download this `supplementary file` and save it to the directory specified by `baseDir`.

```{r}

library(GEOquery)

gse = getGEOSuppFiles("GSE181225", baseDir = getwd())

```
After downloading the file, you should find a new directory named `GSE181225`, which contains the gene expression file. To proceed, read the gene expression file into R. The file has the following structure:
*   The `Symbol` column contains the gene symbols.
*  The remaining four columns represent different samples.
 

```{r}
df = read.delim("GSE181225/GSE181225_LNCaP_p57_VO_and_p57_PIM1_RNA_Seq_normalizedCounts.txt.gz")
head(df)
```

Next, create the __expData__, __geneID__, and __annotation_db__.
```{r}
geneID = df[,1]

expData = df[, 2:ncol(df)]
rownames(expData) = geneID

library(org.Hs.eg.db)
```
To map gene symbols to Entrez identifiers using the __annotation_db__, you can use the `select` function from the `AnnotationDbi` package. Here’s how you can do it in R:

```{r}
library(AnnotationDbi)

genes = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(expData), 
                      columns=c("ENTREZID"), 
                      keytype="SYMBOL")
# initial dimension
print(dim(genes))
head(genes)
```
Remove genes with missing `SYMBOL` or `ENTREZID`.

```{r}
genes = genes[!is.na(genes$SYMBOL), ]
genes = genes[!is.na(genes$ENTREZID), ]

#dimension after dropping missing values
print(dim(genes))
head(genes)
```
 
Remove genes with duplicated `SYMBOL` or `ENTREZID`.
```{r}
genes = genes[!duplicated(genes$SYMBOL),]
genes = genes[!duplicated(genes$ENTREZID), ]
#dimension after dropping missing values
print(dim(genes))
print(head(genes))
```

Keep only rows in __expData__ that have corresponding gene identifiers present in `genes`.

```{r}
expData = data.frame(expData, SYMBOL = rownames(expData))
expData =  merge(expData, genes, by = "SYMBOL")
```

Produce __expData__.
```{r}
rownames(expData) = expData$ENTREZID
expData = expData[, c(2:6)]
print(head(expData))
```

Remove genes with zero variance from __expData__.

```{r}
# Dropping zero variance genes

vars = apply(expData, 1, var)
zeroInd = which(vars == 0)

if(length(zeroInd) != 0) {
  print(paste0("number of zero variance genes ", length(zeroInd)))
  expData = expData[-zeroInd, ]
  genes = genes[-zeroInd, ]
}

print(paste0("number of genes after dropping ", dim(genes)[1]))
```
Remove genes with no gene ontology mapping.

```{r}
## Remove genes with no GO mapping

xx = as.list(org.Hs.egGO[genes$ENTREZID])
haveGO  = sapply(xx,
                 function(x) {if (length(x) == 1 && is.na(x)) FALSE else TRUE })
numNoGO  = sum(!haveGO)
if(numNoGO != 0){
  print(paste0("number of genes with no GO mapping ", length(zeroInd)))
  expData = expData[haveGO, ]
  genes = genes[haveGO, ]
  
}
print(paste0("number of genes after dropping ", dim(genes)[1]))
```
Produce the final __expData__, __geneID__, __annotation_db__. Now, the input is ready for `SGCP`. Refer to
[SGCP Bioconductor page](https://bioconductor.org/packages/release/bioc/html/SGCP.html) in order to see how to use this input in `SGCP`. 

```{r}
expData = expData
print(head(expData))

geneID = genes$ENTREZID
print(head(geneID))

annotation_db = "org.Hs.eg.db" 
```
